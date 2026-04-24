process GRIDSS {

    tag { "${task.process}_${task.attempt}_${id}" }

    publishDir "${params.publishDir}", mode: 'copy'

    container 'gridss/gridss:2.13.2'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        path reference_dir
        path gridss_blacklist
        path gridss_properties

    output:
        path "gridsspl/${id}/", emit: outdir
        tuple val(id), path("gridsspl/${id}/gridss.vcf"), emit: joint_vcf

    script:
    def gridsspl_outdir = "gridsspl/${id}"
    def reference_fasta = "${reference_dir}/reference.fa"
    def jvmheap = "${(task.memory.toGiga().intValue() * 3) / 4}g"
    """
    mkdir -p ${gridsspl_outdir}

    blacklist_arg=""
    reference_fasta_abs="\$(pwd)/${reference_fasta}"
    if [[ -s ${gridss_blacklist} ]]; then
        blacklist_arg="-b ${gridss_blacklist}"
    fi

    properties_arg=""
    if [[ -s ${gridss_properties} ]]; then
        properties_arg="-c ${gridss_properties}"
    fi

    # Run GRIDSS jointly on the matched normal and tumour BAMs.
    set +e
    gridss \
            -o ${gridsspl_outdir}/gridss.vcf \
            -r "\$reference_fasta_abs" \
            --threads ${task.cpus} \
            --jvmheap ${jvmheap} \
            \$blacklist_arg \
            \$properties_arg \
            "${normal_bam}" "${tumor_bam}"
    gridss_exit=\$?
    set -e

    gridss_full_log=""
    shopt -s nullglob
    gridss_logs=(gridss.full.*.log)
    shopt -u nullglob
    if (( \${#gridss_logs[@]} > 0 )); then
        gridss_full_log="\${gridss_logs[0]}"
        cp "\$gridss_full_log" ${gridsspl_outdir}/gridss.full.log
    fi

    shopt -s nullglob
    diagnostic_files=(hs_err_pid*.log replay_pid*.log java_pid*.hprof)
    shopt -u nullglob
    for diagnostic_file in "\${diagnostic_files[@]}"; do
        cp "\$diagnostic_file" ${gridsspl_outdir}/
    done

    if [[ \$gridss_exit -ne 0 ]]; then
        echo "GRIDSS failed for sample ${id} with exit code \$gridss_exit" >&2

        if [[ -s ${gridsspl_outdir}/gridss.full.log ]]; then
            echo "==== tail: ${gridsspl_outdir}/gridss.full.log ====" >&2
            tail -n 200 ${gridsspl_outdir}/gridss.full.log >&2 || true
        fi

        shopt -s nullglob
        copied_diagnostics=(${gridsspl_outdir}/hs_err_pid*.log ${gridsspl_outdir}/replay_pid*.log ${gridsspl_outdir}/java_pid*.hprof)
        shopt -u nullglob
        for diagnostic_file in "\${copied_diagnostics[@]}"; do
            echo "==== tail: \$diagnostic_file ====" >&2
            if [[ "\$diagnostic_file" == *.hprof ]]; then
                ls -lh "\$diagnostic_file" >&2 || true
            else
                tail -n 200 "\$diagnostic_file" >&2 || true
            fi
        done

        if [[ ! -s ${gridsspl_outdir}/gridss.full.log && \${#copied_diagnostics[@]} -eq 0 ]]; then
            echo "No GRIDSS/JVM diagnostic files found in the task working directory." >&2
        fi

        exit \$gridss_exit
    fi
    """
}

process GRIDSS_SOMATIC_FILTER_WITH_PON {

    tag { "${task.process}_${task.attempt}_${id}" }

    publishDir "${params.publishDir}", mode: 'copy'

    container 'gridss/gridss:2.13.2'

    input:
        tuple val(id), path(gridss_joint_vcf)
        path gridss_pon_dir

    output:
        tuple val(id), path("gridsspl/${id}/${id}_somatic_filtered.vcf"), emit: joint_vcf

    script:
    def gridsspl_outdir = "gridsspl/${id}"
    def filtered_vcf = "${gridsspl_outdir}/${id}_somatic_filtered.vcf"
    """
    set -euo pipefail

    mkdir -p ${gridsspl_outdir}

    gridss_somatic_filter \
        --input ${gridss_joint_vcf} \
        --output ${filtered_vcf} \
        --pondir ${gridss_pon_dir} \
        --normalordinal 1 \
        --tumourordinal 2
    """
}

process GRIDSS_SOMATIC_FILTER_NO_PON {

    tag { "${task.process}_${task.attempt}_${id}" }

    publishDir "${params.publishDir}", mode: 'copy'

    container 'gridss/gridss:2.13.2'

    input:
        tuple val(id), path(gridss_joint_vcf)

    output:
        tuple val(id), path("gridsspl/${id}/${id}_somatic_filtered.vcf"), emit: joint_vcf

    script:
    def gridsspl_outdir = "gridsspl/${id}"
    def filtered_vcf = "${gridsspl_outdir}/${id}_somatic_filtered.vcf"
    """
    set -euo pipefail

    mkdir -p ${gridsspl_outdir}

    gridss_somatic_filter \
        --input ${gridss_joint_vcf} \
        --output ${filtered_vcf} \
        --normalordinal 1 \
        --tumourordinal 2
    """
}

process GRIDSS_TUMOUR_VCF {

    tag { "${task.process}_${task.attempt}_${id}" }

    label 'utility'

    publishDir "${params.publishDir}", mode: 'copy'

    container 'gridss/gridss:2.13.2'

    input:
        tuple val(id), path(gridss_filtered_joint_vcf)

    output:
        tuple val(id), path("gridsspl/${id}/${id}_tumor.vcf"), emit: vcf_t

    script:
    def gridsspl_outdir = "gridsspl/${id}"
    def tumour_vcf = "${gridsspl_outdir}/${id}_tumor.vcf"
    """
    set -euo pipefail

    mkdir -p ${gridsspl_outdir}

    GRIDSS_INPUT="${gridss_filtered_joint_vcf}" GRIDSS_OUTPUT="${tumour_vcf}" Rscript - <<'RSCRIPT'
    suppressPackageStartupMessages(library(VariantAnnotation))

    input_path <- Sys.getenv("GRIDSS_INPUT")
    output_path <- Sys.getenv("GRIDSS_OUTPUT")

    vcf <- readVcf(input_path)
    sample_names <- samples(header(vcf))

    if (length(sample_names) < 2) {
        stop(sprintf("Expected a two-sample GRIDSS VCF but found %d sample(s) in %s", length(sample_names), input_path))
    }

    tumour_vcf <- vcf[, 2, drop = FALSE]
    writeVcf(tumour_vcf, output_path)
RSCRIPT
    """
}

process GRIDSS_NORMAL_VCF {

    tag { "${task.process}_${task.attempt}_${id}" }

    label 'utility'

    publishDir "${params.publishDir}", mode: 'copy'

    container 'gridss/gridss:2.13.2'

    input:
        tuple val(id), path(gridss_joint_vcf)

    output:
        path "gridsspl/${id}/${id}_normal.vcf", emit: vcf_n

    script:
    def gridsspl_outdir = "gridsspl/${id}"
    def normal_vcf = "${gridsspl_outdir}/${id}_normal.vcf"
    """
    set -euo pipefail

    mkdir -p ${gridsspl_outdir}

    GRIDSS_INPUT="${gridss_joint_vcf}" GRIDSS_OUTPUT="${normal_vcf}" Rscript - <<'RSCRIPT'
    suppressPackageStartupMessages(library(VariantAnnotation))
    suppressPackageStartupMessages(library(S4Vectors))

    input_path <- Sys.getenv("GRIDSS_INPUT")
    output_path <- Sys.getenv("GRIDSS_OUTPUT")

    vcf <- readVcf(input_path)
    sample_names <- samples(header(vcf))

    if (length(sample_names) < 2) {
        stop(sprintf("Expected a two-sample GRIDSS VCF but found %d sample(s) in %s", length(sample_names), input_path))
    }

    genotype <- geno(vcf)
    required_fields <- c("REF", "REFPAIR", "VF", "BVF")
    missing_fields <- setdiff(required_fields, names(genotype))

    if (length(missing_fields) > 0) {
        stop(sprintf(
            "Missing GRIDSS FORMAT field(s) required for normal VCF extraction: %s",
            paste(missing_fields, collapse = ", ")
        ))
    }

    extract_numeric_field <- function(field_name, sample_index) {
        values <- genotype[[field_name]][, sample_index, drop = TRUE]
        values <- as.numeric(values)
        values[is.na(values)] <- 0
        values
    }

    safe_af <- function(variant_support, ref_support, refpair_support) {
        depth <- variant_support + ref_support + refpair_support
        af <- ifelse(depth > 0, variant_support / depth, 0)
        af[is.na(af)] <- 0
        af
    }

    mateid <- info(vcf)$MATEID
    if (is.null(mateid)) {
        stop("Missing INFO/MATEID in GRIDSS VCF; cannot distinguish breakpoint and single-breakend records.")
    }

    if (inherits(mateid, "CharacterList")) {
        is_breakpoint <- elementNROWS(mateid) > 0
    } else {
        mateid_values <- as.character(mateid)
        is_breakpoint <- !is.na(mateid_values) & nzchar(mateid_values)
    }

    is_single_breakend <- !is_breakpoint
    is_pass <- as.character(fixed(vcf)$FILTER) == "PASS"

    normal_ref <- extract_numeric_field("REF", 1)
    normal_refpair <- extract_numeric_field("REFPAIR", 1)
    tumour_ref <- extract_numeric_field("REF", 2)
    tumour_refpair <- extract_numeric_field("REFPAIR", 2)

    normal_vf <- extract_numeric_field("VF", 1)
    tumour_vf <- extract_numeric_field("VF", 2)
    normal_bvf <- extract_numeric_field("BVF", 1)
    tumour_bvf <- extract_numeric_field("BVF", 2)

    normal_bp_af <- safe_af(normal_vf, normal_ref, normal_refpair)
    tumour_bp_af <- safe_af(tumour_vf, tumour_ref, tumour_refpair)
    normal_be_af <- safe_af(normal_bvf, normal_ref, normal_refpair)
    tumour_be_af <- safe_af(tumour_bvf, tumour_ref, tumour_refpair)

    keep_breakpoints <- is_breakpoint &
        is_pass &
        normal_bp_af >= 0.10 &
        normal_vf >= 4 &
        tumour_bp_af <= (3 * normal_bp_af)

    keep_single_breakends <- is_single_breakend &
        is_pass &
        normal_be_af >= 0.10 &
        normal_bvf >= 4 &
        tumour_be_af <= (3 * normal_be_af)

    keep_in_normal <- keep_breakpoints | keep_single_breakends
    normal_vcf <- vcf[keep_in_normal, 1, drop = FALSE]

    writeVcf(normal_vcf, output_path)
RSCRIPT
    """
}
