#!/usr/bin/env nextflow
// Startup: source .startup.sh
// Command: nextflow run main.nf -profile singularity 
// Test: SVABA[13m]; DELLY[1m]; LUMPY[3m]; MANTA[5m]; GRIDSSPL[36m]

nextflow.enable.dsl=2

include { JASMINE_PANEL as MANTA_NORMAL_PANEL } from './modules/jasmine_panel.nf'
include { JASMINE_PANEL as LUMPY_NORMAL_PANEL } from './modules/jasmine_panel.nf'
include { JASMINE_PANEL as SVABA_NORMAL_PANEL } from './modules/jasmine_panel.nf'
include { JASMINE_PANEL as DELLY_NORMAL_PANEL } from './modules/jasmine_panel.nf'
include { JASMINE_PANEL as GRIDSS_NORMAL_PANEL } from './modules/jasmine_panel.nf'

include { TUMOUR_FILTER as MANTA_TUMOUR_FILTER } from './modules/tumour_filter.nf'
include { TUMOUR_FILTER as LUMPY_TUMOUR_FILTER } from './modules/tumour_filter.nf'
include { TUMOUR_FILTER as SVABA_TUMOUR_FILTER } from './modules/tumour_filter.nf'
include { TUMOUR_FILTER as DELLY_TUMOUR_FILTER } from './modules/tumour_filter.nf'
include { TUMOUR_FILTER as GRIDSS_TUMOUR_FILTER } from './modules/tumour_filter.nf'

include { TUMOUR_CONSENSUS_VCF } from './modules/tumour_consensus.nf'
include { TUMOUR_CONSENSUS_CALL } from './modules/tumour_consensus.nf'
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
    def usingDerivedBlacklist = params.gridss_blacklist == null && buildDefaults.gridss_blacklist != null
    def exactPonBuild = ['hg38', 'hg19'].contains(requestedBuild) ? requestedBuild : null

    if (!['hg38', 'hg19', 'other'].contains(normalizedBuild)) {
        error "Unsupported --genome_build '${params.genome_build}'. Use hg38, hg19, or other."
    }

    def csvValue = params.csv
    def referenceFastaValue = params.reference_fasta ?: buildDefaults.reference_fasta
    def referenceIndexValue = params.reference_index
    def gridssBlacklistValue = params.gridss_blacklist ?: buildDefaults.gridss_blacklist
    def gridssPropertiesValue = params.gridss_properties ?: params.default_gridss_properties
    def gridssSvPonValue = params.gridss_sv_pon
    def gridssSglPonValue = params.gridss_sgl_pon
    def gridssResourceBundleDefault = exactPonBuild ? ((params.reference_defaults ?: [:])[exactPonBuild] ?: [:]).gridss_resource_bundle : null
    def gridssResourceBundleValue = params.gridss_resource_bundle ?: gridssResourceBundleDefault
    def hasExplicitGridssPon = gridssSvPonValue != null || gridssSglPonValue != null
    def useAutoGridssPon = !hasExplicitGridssPon && exactPonBuild != null
    def useGridssPon = hasExplicitGridssPon || useAutoGridssPon

    if (!csvValue) {
        error "Missing required parameter --csv. Provide a sample sheet or use -profile test."
    }

    if (!referenceFastaValue) {
        error "Missing required reference FASTA. Set --genome_build hg38/hg19 or provide --reference_fasta."
    }

    if (referenceFastaValue.toString().endsWith('.gz') && referenceIndexValue) {
        log.warn "Ignoring --reference_index for compressed FASTA inputs because the pipeline rebuilds the .fai after decompression."
    }

    if (!(csvValue.toString() ==~ /^(https?|ftp):\/\/.+/) && !new File(csvValue.toString()).exists()) {
        error "Parameter --csv points to a missing local path: ${csvValue}"
    }

    if (!(referenceFastaValue.toString() ==~ /^(https?|ftp):\/\/.+/) && !new File(referenceFastaValue.toString()).exists()) {
        error "Parameter --reference_fasta points to a missing local path: ${referenceFastaValue}"
    }

    if (referenceIndexValue && !(referenceIndexValue.toString() ==~ /^(https?|ftp):\/\/.+/) && !new File(referenceIndexValue.toString()).exists()) {
        error "Parameter --reference_index points to a missing local path: ${referenceIndexValue}"
    }

    if (gridssBlacklistValue && !(gridssBlacklistValue.toString() ==~ /^(https?|ftp):\/\/.+/) && !new File(gridssBlacklistValue.toString()).exists()) {
        error "Parameter --gridss_blacklist points to a missing local path: ${gridssBlacklistValue}"
    }

    if (gridssPropertiesValue && !(gridssPropertiesValue.toString() ==~ /^(https?|ftp):\/\/.+/) && !new File(gridssPropertiesValue.toString()).exists()) {
        error "Parameter --gridss_properties points to a missing local path: ${gridssPropertiesValue}"
    }

    if ((gridssSvPonValue == null) != (gridssSglPonValue == null)) {
        error "Provide both --gridss_sv_pon and --gridss_sgl_pon together, or omit both to use the build-derived GRIDSS resource bundle."
    }

    if (gridssSvPonValue && !(gridssSvPonValue.toString() ==~ /^(https?|ftp):\/\/.+/) && !new File(gridssSvPonValue.toString()).exists()) {
        error "Parameter --gridss_sv_pon points to a missing local path: ${gridssSvPonValue}"
    }

    if (gridssSglPonValue && !(gridssSglPonValue.toString() ==~ /^(https?|ftp):\/\/.+/) && !new File(gridssSglPonValue.toString()).exists()) {
        error "Parameter --gridss_sgl_pon points to a missing local path: ${gridssSglPonValue}"
    }

    if (useAutoGridssPon) {
        if (!gridssResourceBundleValue) {
            error "Missing GRIDSS resource bundle for --genome_build ${requestedBuild}. Provide --gridss_resource_bundle or both explicit GRIDSS PoN files."
        }

        if (!(gridssResourceBundleValue.toString() ==~ /^(https?|ftp):\/\/.+/) && !new File(gridssResourceBundleValue.toString()).exists()) {
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

    def reference_fasta_source = Channel.value(file(referenceFastaValue))
    def gridss_properties_source = Channel.value(file(gridssPropertiesValue))
    def prepared_reference
    def prepared_gridss_blacklist
    def prepared_gridss_properties
    def prepared_gridss_pondir = null

    if (referenceIndexValue) {
        def reference_index_source = Channel.value(file(referenceIndexValue))
        PREP_REFERENCE_WITH_INDEX(reference_fasta_source, reference_index_source)
        prepared_reference = PREP_REFERENCE_WITH_INDEX.out.ref_dir
    } else {
        PREP_REFERENCE(reference_fasta_source)
        prepared_reference = PREP_REFERENCE.out.ref_dir
    }

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
    
    // Run callers on BAM
    MANTA(bam_channel, prepared_reference)
    LUMPY(bam_channel, prepared_reference)
    SVABA(bam_channel, prepared_reference)
    DELLY(bam_channel, prepared_reference)
    DELLY_SPLIT(DELLY.out.id, DELLY.out.bcf)
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

    // Combine VCFs from normal outputs
    manta_normals = MANTA.out.vcf_n.collect()
    manta_panel = MANTA_NORMAL_PANEL(manta_normals, 'manta')
    lumpy_normals = LUMPY.out.vcf_n.collect()
    lumpy_panel = LUMPY_NORMAL_PANEL(lumpy_normals, 'lumpy')
    svaba_normals = SVABA.out.vcf_n.collect()
    svaba_panel = SVABA_NORMAL_PANEL(svaba_normals, 'svaba')
    delly_normals = DELLY_SPLIT.out.vcf_n.collect()
    delly_panel = DELLY_NORMAL_PANEL(delly_normals, 'delly')
    gridss_normals = GRIDSS_NORMAL_VCF.out.vcf_n.collect()
    gridss_panel = GRIDSS_NORMAL_PANEL(gridss_normals, 'gridss')

    // Filter tumor VCFs using panels
    MANTA_TUMOUR_FILTER(MANTA.out.vcf_t, manta_panel, 'manta')
    LUMPY_TUMOUR_FILTER(LUMPY.out.vcf_t, lumpy_panel, 'lumpy')
    SVABA_TUMOUR_FILTER(SVABA.out.vcf_t, svaba_panel, 'svaba')
    DELLY_TUMOUR_FILTER(DELLY_SPLIT.out.vcf_t, delly_panel, 'delly')
    GRIDSS_TUMOUR_FILTER(GRIDSS_TUMOUR_VCF.out.vcf_t, gridss_panel, 'gridss')
    
    // Combine all tumor channels by Sample ID
    TUMOR_COMBINED = MANTA_TUMOUR_FILTER.out.filtered_vcf
        .combine(SVABA_TUMOUR_FILTER.out.filtered_vcf, by: 0)
        .combine(LUMPY_TUMOUR_FILTER.out.filtered_vcf, by: 0)
        .combine(DELLY_TUMOUR_FILTER.out.filtered_vcf, by: 0)
        .combine(GRIDSS_TUMOUR_FILTER.out.filtered_vcf, by: 0)
    
    // Aggregate VCF files for each sample
    //  Create channel with (val(id), [path(vcf), path(vcf), ...])
    TUMOR_FORMATTED = TUMOR_COMBINED.map { v -> tuple(v[0], v.subList(1, v.size())) } // `*vcfs` captures all elements except the first
    //  Create single VCF file for each sample
    TUMOUR_CONSENSUS_VCF(TUMOR_FORMATTED)

    // Call consensus SVs
    TUMOUR_CONSENSUS_CALL(TUMOUR_CONSENSUS_VCF.out.id, TUMOUR_CONSENSUS_VCF.out.consensus_vcf)

}
