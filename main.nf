#!/usr/bin/env nextflow
// Startup: source .startup.sh
// Command: nextflow run main.nf -profile singularity
// Test: SVABA[13m]; DELLY[1m]; LUMPY[3m]; MANTA[5m]; GRIDSSPL[36m]

nextflow.enable.dsl=2

include { SV_COHORT } from './modules/cohort.nf'
include { PREP_REFERENCE; PREP_REFERENCE_WITH_INDEX; PREP_GRIDSS_ASSETS; PREP_GRIDSS_ASSETS_WITH_BLACKLIST; PREP_GRIDSS_PON_FROM_FILES; PREP_GRIDSS_PON_FROM_BUNDLE } from './modules/reference_prep.nf'
include { MANTA } from './modules/manta.nf'
include { LUMPY } from './modules/lumpy.nf'
include { SVABA } from './modules/svaba.nf'
include { DELLY; DELLY_SPLIT } from './modules/delly.nf'
include { GRIDSS; GRIDSS_SOMATIC_FILTER_WITH_PON; GRIDSS_SOMATIC_FILTER_NO_PON; GRIDSS_TUMOUR_VCF; GRIDSS_NORMAL_VCF } from './modules/gridss.nf'

workflow {
    def requestedBuild = (params.genome_build ?: 'hg38').toString().toLowerCase()
    def buildAliases = params.genome_build_aliases ?: [:]
    def normalizedBuild = buildAliases.containsKey(requestedBuild) ? buildAliases[requestedBuild] : requestedBuild
    def buildDefaults = (params.reference_defaults ?: [:])[normalizedBuild] ?: [:]

    if (!['hg38', 'hg19', 'other'].contains(normalizedBuild)) {
        error "Unsupported --genome_build '${params.genome_build}'. Use hg38, hg19, or other."
    }

    def csvValue = params.csv
    def referenceFastaValue = params.reference_fasta ?: buildDefaults.reference_fasta
    def referenceIndexValue = params.reference_index

    // full: call, then filter and merge; call: raw caller VCFs only; cohort: filter and merge VCFs from a long-format sheet
    def entry = (params.entry ?: 'full').toString().toLowerCase()
    if (!['full', 'call', 'cohort'].contains(entry)) {
        error "Unsupported --entry '${params.entry}'. Use full, call, or cohort."
    }

    if (!csvValue) {
        error "Missing required parameter --csv. Provide a sample sheet or use -profile test."
    }

    if (entry != 'cohort' && !referenceFastaValue) {
        error "Missing required reference FASTA. Set --genome_build hg38/hg19 or provide --reference_fasta."
    }

    if (entry != 'cohort' && referenceFastaValue.toString().endsWith('.gz') && referenceIndexValue) {
        log.warn "Ignoring --reference_index for compressed FASTA inputs because the pipeline rebuilds the .fai after decompression."
    }

    if (!(csvValue.toString() ==~ /^(https?|ftp|s3):\/\/.+/) && !new File(csvValue.toString()).exists()) {
        error "Parameter --csv points to a missing local path: ${csvValue}"
    }

    if (entry != 'cohort' && !(referenceFastaValue.toString() ==~ /^(https?|ftp|s3):\/\/.+/) && !new File(referenceFastaValue.toString()).exists()) {
        error "Parameter --reference_fasta points to a missing local path: ${referenceFastaValue}"
    }

    if (entry != 'cohort' && referenceIndexValue && !(referenceIndexValue.toString() ==~ /^(https?|ftp|s3):\/\/.+/) && !new File(referenceIndexValue.toString()).exists()) {
        error "Parameter --reference_index points to a missing local path: ${referenceIndexValue}"
    }

    // Caller toggles + consensus threshold
    // .toString().toBoolean() guards against Nextflow parsing CLI overrides like
    // `--run_gridss false` as the (truthy) String "false" rather than Boolean false.
    def runManta = params.run_manta.toString().toBoolean()
    def runLumpy = params.run_lumpy.toString().toBoolean()
    def runSvaba = params.run_svaba.toString().toBoolean()
    def runDelly = params.run_delly.toString().toBoolean()
    def runGridss = params.run_gridss.toString().toBoolean()
    def enabledCallerCount = [runManta, runLumpy, runSvaba, runDelly, runGridss].count { it }

    if (entry == 'full' && enabledCallerCount < 2) {
        error "Consensus calling requires at least 2 enabled SV callers (got ${enabledCallerCount}). Enable more via --run_manta/--run_lumpy/--run_svaba/--run_delly/--run_gridss."
    }


    def rows = null

    if (entry != 'cohort') {
    def reference_fasta_source = Channel.value(file(referenceFastaValue))
    def prepared_reference
    def prepared_gridss_blacklist = null
    def prepared_gridss_properties = null
    def prepared_gridss_pondir = null
    def useGridssPon = false

    if (referenceIndexValue) {
        def reference_index_source = Channel.value(file(referenceIndexValue))
        PREP_REFERENCE_WITH_INDEX(reference_fasta_source, reference_index_source)
        prepared_reference = PREP_REFERENCE_WITH_INDEX.out.ref_dir
    } else {
        PREP_REFERENCE(reference_fasta_source)
        prepared_reference = PREP_REFERENCE.out.ref_dir
    }

    if (runGridss) {
        def usingDerivedBlacklist = params.gridss_blacklist == null && buildDefaults.gridss_blacklist != null
        def exactPonBuild = ['hg38', 'hg19'].contains(requestedBuild) ? requestedBuild : null

        def gridssBlacklistValue = params.gridss_blacklist ?: buildDefaults.gridss_blacklist
        def gridssPropertiesValue = params.gridss_properties ?: params.default_gridss_properties
        def gridssSvPonValue = params.gridss_sv_pon
        def gridssSglPonValue = params.gridss_sgl_pon
        def gridssResourceBundleDefault = exactPonBuild ? ((params.reference_defaults ?: [:])[exactPonBuild] ?: [:]).gridss_resource_bundle : null
        def gridssResourceBundleValue = params.gridss_resource_bundle ?: gridssResourceBundleDefault
        def hasExplicitGridssPon = gridssSvPonValue != null || gridssSglPonValue != null
        def useAutoGridssPon = !hasExplicitGridssPon && exactPonBuild != null
        useGridssPon = hasExplicitGridssPon || useAutoGridssPon

        if (gridssBlacklistValue && !(gridssBlacklistValue.toString() ==~ /^(https?|ftp|s3):\/\/.+/) && !new File(gridssBlacklistValue.toString()).exists()) {
            error "Parameter --gridss_blacklist points to a missing local path: ${gridssBlacklistValue}"
        }

        if (gridssPropertiesValue && !(gridssPropertiesValue.toString() ==~ /^(https?|ftp|s3):\/\/.+/) && !new File(gridssPropertiesValue.toString()).exists()) {
            error "Parameter --gridss_properties points to a missing local path: ${gridssPropertiesValue}"
        }

        if ((gridssSvPonValue == null) != (gridssSglPonValue == null)) {
            error "Provide both --gridss_sv_pon and --gridss_sgl_pon together, or omit both to use the build-derived GRIDSS resource bundle."
        }

        if (gridssSvPonValue && !(gridssSvPonValue.toString() ==~ /^(https?|ftp|s3):\/\/.+/) && !new File(gridssSvPonValue.toString()).exists()) {
            error "Parameter --gridss_sv_pon points to a missing local path: ${gridssSvPonValue}"
        }

        if (gridssSglPonValue && !(gridssSglPonValue.toString() ==~ /^(https?|ftp|s3):\/\/.+/) && !new File(gridssSglPonValue.toString()).exists()) {
            error "Parameter --gridss_sgl_pon points to a missing local path: ${gridssSglPonValue}"
        }

        if (useAutoGridssPon) {
            if (!gridssResourceBundleValue) {
                error "Missing GRIDSS resource bundle for --genome_build ${requestedBuild}. Provide --gridss_resource_bundle or both explicit GRIDSS PoN files."
            }

            if (!(gridssResourceBundleValue.toString() ==~ /^(https?|ftp|s3):\/\/.+/) && !new File(gridssResourceBundleValue.toString()).exists()) {
                error "Parameter --gridss_resource_bundle points to a missing local path: ${gridssResourceBundleValue}"
            }
        } else if (!hasExplicitGridssPon && normalizedBuild != 'other') {
            error "Automatic GRIDSS PoN extraction is only available for exact --genome_build hg38 and --genome_build hg19. Provide both --gridss_sv_pon and --gridss_sgl_pon for aliases or other builds."
        } else if (!useGridssPon) {
            log.warn "No GRIDSS PoN provided for --genome_build other. The workflow will run gridss_somatic_filter without --pondir."
        }

        if (params.reference_fasta && usingDerivedBlacklist) {
            log.warn "Using a build-derived GRIDSS blacklist with a custom reference FASTA. The workflow will remap simple contig aliases against the staged reference, but you should still override --gridss_blacklist if your FASTA uses a non-canonical naming scheme."
        }

        def gridss_properties_source = Channel.value(file(gridssPropertiesValue))

        if (gridssBlacklistValue) {
            def gridss_blacklist_source = Channel.value(file(gridssBlacklistValue))
            PREP_GRIDSS_ASSETS_WITH_BLACKLIST(prepared_reference, gridss_blacklist_source, gridss_properties_source)
            prepared_gridss_blacklist = PREP_GRIDSS_ASSETS_WITH_BLACKLIST.out.blacklist
            prepared_gridss_properties = PREP_GRIDSS_ASSETS_WITH_BLACKLIST.out.properties
        } else {
            PREP_GRIDSS_ASSETS(gridss_properties_source)
            prepared_gridss_blacklist = PREP_GRIDSS_ASSETS.out.blacklist
            prepared_gridss_properties = PREP_GRIDSS_ASSETS.out.properties
        }

        if (hasExplicitGridssPon) {
            def gridss_sv_pon_source = Channel.value(file(gridssSvPonValue))
            def gridss_sgl_pon_source = Channel.value(file(gridssSglPonValue))
            PREP_GRIDSS_PON_FROM_FILES(gridss_sv_pon_source, gridss_sgl_pon_source)
            prepared_gridss_pondir = PREP_GRIDSS_PON_FROM_FILES.out.pondir
        } else if (useAutoGridssPon) {
            def gridss_resource_bundle_source = Channel.value(file(gridssResourceBundleValue))
            PREP_GRIDSS_PON_FROM_BUNDLE(gridss_resource_bundle_source, requestedBuild)
            prepared_gridss_pondir = PREP_GRIDSS_PON_FROM_BUNDLE.out.pondir
        }
    } else {
        log.info "Skipping GRIDSS (--run_gridss=false): no GRIDSS asset validation or preparation will run."
    }

    // Load inputs from csv
    def bam_channel = Channel
            .fromPath(csvValue)
            .splitCsv(header: true, sep: ',', strip: true)
            .map { row -> tuple(
                row.id,
                file(row.tumor_bam),
                file(row.tumor_bai),
                file(row.normal_bam),
                file(row.normal_bai)
                ) }

    // The normal's identity decides how many distinct normals the panel sees. It defaults to the
    // normal BAM file name, so pairs that share a normal BAM share a normal.
    def pair_normals = Channel
            .fromPath(csvValue)
            .splitCsv(header: true, sep: ',', strip: true)
            .map { row -> tuple(row.id, (row.normal_id ?: file(row.normal_bam).name.replaceFirst(/\.(bam|cram)$/, '')).toString()) }

    // Run callers on BAM, each gated by its --run_<caller> toggle. Each caller adds its tumour and
    // normal VCFs to the rows passed to the cohort stage as (id, caller, role, vcf).
    def caller_rows = []

    if (runManta) {
        MANTA(bam_channel, prepared_reference)
        caller_rows << MANTA.out.vcf_t.map { id, vcf -> tuple(id, 'manta', 'tumor', vcf) }
        caller_rows << MANTA.out.vcf_n.map { id, vcf -> tuple(id, 'manta', 'normal', vcf) }
    } else {
        log.info "Skipping Manta (--run_manta=false)"
    }

    if (runLumpy) {
        LUMPY(bam_channel, prepared_reference)
        caller_rows << LUMPY.out.vcf_t.map { id, vcf -> tuple(id, 'lumpy', 'tumor', vcf) }
        caller_rows << LUMPY.out.vcf_n.map { id, vcf -> tuple(id, 'lumpy', 'normal', vcf) }
    } else {
        log.info "Skipping LUMPY (--run_lumpy=false)"
    }

    if (runSvaba) {
        SVABA(bam_channel, prepared_reference)
        caller_rows << SVABA.out.vcf_t.map { id, vcf -> tuple(id, 'svaba', 'tumor', vcf) }
        caller_rows << SVABA.out.vcf_n.map { id, vcf -> tuple(id, 'svaba', 'normal', vcf) }
    } else {
        log.info "Skipping SvABA (--run_svaba=false)"
    }

    if (runDelly) {
        DELLY(bam_channel, prepared_reference)
        DELLY_SPLIT(DELLY.out.id, DELLY.out.bcf)
        if (params.delly_filter.toString().toBoolean()) {
            // delly filter replaces the split tumour VCF; the cohort stage runs it on the joint BCF
            caller_rows << DELLY.out.bcf_t.map { id, bcf -> tuple(id, 'delly', 'bcf', bcf) }
        } else {
            caller_rows << DELLY_SPLIT.out.vcf_t.map { id, vcf -> tuple(id, 'delly', 'tumor', vcf) }
        }
        caller_rows << DELLY_SPLIT.out.vcf_n.map { id, vcf -> tuple(id, 'delly', 'normal', vcf) }
    } else {
        log.info "Skipping DELLY (--run_delly=false)"
    }

    if (runGridss) {
        GRIDSS(bam_channel, prepared_reference, prepared_gridss_blacklist, prepared_gridss_properties)
        def gridss_filtered_joint_vcf
        if (useGridssPon) {
            GRIDSS_SOMATIC_FILTER_WITH_PON(GRIDSS.out.joint_vcf, prepared_gridss_pondir)
            gridss_filtered_joint_vcf = GRIDSS_SOMATIC_FILTER_WITH_PON.out.joint_vcf
        } else {
            GRIDSS_SOMATIC_FILTER_NO_PON(GRIDSS.out.joint_vcf)
            gridss_filtered_joint_vcf = GRIDSS_SOMATIC_FILTER_NO_PON.out.joint_vcf
        }
        GRIDSS_NORMAL_VCF(GRIDSS.out.joint_vcf)
        GRIDSS_TUMOUR_VCF(gridss_filtered_joint_vcf)
        caller_rows << GRIDSS_TUMOUR_VCF.out.vcf_t.map { id, vcf -> tuple(id, 'gridss', 'tumor', vcf) }
        caller_rows << GRIDSS_NORMAL_VCF.out.vcf_n.map { id, vcf -> tuple(id, 'gridss', 'normal', vcf) }
    } else {
        log.info "Skipping GRIDSS (--run_gridss=false)"
    }

    if (entry == 'full') {
        rows = caller_rows
            .inject { a, b -> a.mix(b) }
            .combine(pair_normals, by: 0)
            .map { id, caller, role, vcf, nid -> tuple(id, nid, caller, role, vcf) }
    }
    }

    if (entry == 'cohort') {
        // Long-format sheet, one row per VCF: id, normal_id, caller, role, vcf
        def sheet = file(csvValue).splitCsv(header: true, sep: ',', strip: true)
        def required = ['id', 'caller', 'role', 'vcf']
        def missing = required.findAll { !(sheet ? sheet[0].containsKey(it) : true) }
        if (missing) {
            error "Cohort sheet ${csvValue} is missing required columns: ${missing.join(', ')}. Expected id, normal_id, caller, role, vcf."
        }
        if (!sheet) {
            error "Cohort sheet ${csvValue} has no rows."
        }

        def seen = [:]
        sheet.eachWithIndex { row, i ->
            def where = "row ${i + 2} of ${csvValue}"
            def role = (row.role ?: '').toString().toLowerCase()
            def caller = (row.caller ?: '').toString()
            if (!['tumor', 'normal', 'panel', 'bcf'].contains(role)) {
                error "Invalid role '${row.role}' in ${where}. Use tumor, normal, panel, or bcf."
            }
            if (!(caller ==~ /^[A-Za-z0-9-]+$/)) {
                error "Invalid caller '${row.caller}' in ${where}. Use letters, digits, and hyphens only."
            }
            if (!row.vcf) {
                error "Missing vcf in ${where}."
            }
            if (role == 'bcf' && caller != 'delly') {
                error "Role bcf is only for the DELLY joint BCF; got caller '${caller}' in ${where}."
            }
            if (role != 'panel' && !row.id) {
                error "Missing id in ${where}."
            }
            if (role == 'normal' && !row.normal_id) {
                error "Missing normal_id in ${where}. Normal rows need an identifier for the normal sample so a normal shared by several pairs is counted once."
            }
            def key = [role, row.id ?: '', caller]
            if (role != 'panel' && seen.containsKey(key)) {
                error "Duplicate ${role} row for id '${row.id}' and caller '${caller}' in ${where}."
            }
            if (role == 'panel' && seen.containsKey([role, '', caller])) {
                error "More than one panel row for caller '${caller}' in ${where}."
            }
            seen[role == 'panel' ? [role, '', caller] : key] = true
        }

        // With delly_filter on, a DELLY BCF row replaces that pair's DELLY tumour VCF row
        def useDellyFilter = params.delly_filter.toString().toBoolean()
        def bcfIds = sheet.findAll { it.role.toString().toLowerCase() == 'bcf' }.collect { it.id.toString() } as Set
        sheet = sheet.findAll { row ->
            def role = row.role.toString().toLowerCase()
            if (role == 'bcf') { return useDellyFilter }
            return !(useDellyFilter && role == 'tumor' && row.caller.toString() == 'delly' && bcfIds.contains(row.id.toString()))
        }

        def tumourCallers = sheet.findAll { it.role.toString().toLowerCase() in ['tumor', 'bcf'] }.collect { it.caller.toString() }.unique()
        if (tumourCallers.size() < 2) {
            error "Consensus calling needs tumor rows from at least 2 callers (found ${tumourCallers.size()}: ${tumourCallers})."
        }

        rows = Channel.fromList(sheet.collect { row ->
            tuple((row.id ?: '').toString(), (row.normal_id ?: '').toString(), row.caller.toString(), row.role.toString().toLowerCase(), file(row.vcf))
        })
    }

    if (entry != 'call') {
        SV_COHORT(rows)
    }
}
