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

    # The GRIDSS 2.13.2 wrapper can carry a CRLF Rscript shebang in some
    # container builds, so invoke the R script explicitly.
    gridss_somatic_filter_script="\$(command -v gridss_somatic_filter)"
    gridss_somatic_filter_lib="\$(find /opt/gridss /usr/local -name libgridss.R -print -quit 2>/dev/null || true)"
    if [[ -z "\$gridss_somatic_filter_lib" ]]; then
        echo "Could not find libgridss.R required by gridss_somatic_filter" >&2
        exit 1
    fi
    gridss_somatic_filter_scriptdir="\$(dirname "\$gridss_somatic_filter_lib")"

    echo "gridss_somatic_filter: \$gridss_somatic_filter_script" >&2
    echo "gridss_somatic_filter scriptdir: \$gridss_somatic_filter_scriptdir" >&2
    echo "GRIDSS input VCF:" >&2
    ls -lh ${gridss_joint_vcf} >&2
    echo "GRIDSS PoN directory:" >&2
    ls -lh ${gridss_pon_dir} >&2

    Rscript "\$gridss_somatic_filter_script" \
        --scriptdir "\$gridss_somatic_filter_scriptdir" \
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

    # The GRIDSS 2.13.2 wrapper can carry a CRLF Rscript shebang in some
    # container builds, so invoke the R script explicitly.
    gridss_somatic_filter_script="\$(command -v gridss_somatic_filter)"
    gridss_somatic_filter_lib="\$(find /opt/gridss /usr/local -name libgridss.R -print -quit 2>/dev/null || true)"
    if [[ -z "\$gridss_somatic_filter_lib" ]]; then
        echo "Could not find libgridss.R required by gridss_somatic_filter" >&2
        exit 1
    fi
    gridss_somatic_filter_scriptdir="\$(dirname "\$gridss_somatic_filter_lib")"

    echo "gridss_somatic_filter: \$gridss_somatic_filter_script" >&2
    echo "gridss_somatic_filter scriptdir: \$gridss_somatic_filter_scriptdir" >&2
    echo "GRIDSS input VCF:" >&2
    ls -lh ${gridss_joint_vcf} >&2

    Rscript "\$gridss_somatic_filter_script" \
        --scriptdir "\$gridss_somatic_filter_scriptdir" \
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

    awk '
    BEGIN {
        FS = OFS = "\\t"
    }
    function die(message) {
        print message > "/dev/stderr"
        exit 1
    }
    /^##/ {
        print
        next
    }
    /^#CHROM/ {
        if (NF < 11) {
            die("Expected at least two sample columns in GRIDSS VCF header")
        }
        print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$11
        saw_header = 1
        next
    }
    {
        if (!saw_header) {
            die("Missing #CHROM header in GRIDSS VCF")
        }
        if (NF < 11) {
            die("Expected at least two sample columns in GRIDSS VCF record")
        }
        print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$11
    }
    END {
        if (!saw_header) {
            die("Missing #CHROM header in GRIDSS VCF")
        }
    }' ${gridss_filtered_joint_vcf} > ${tumour_vcf}
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

    awk '
    BEGIN {
        FS = OFS = "\\t"
    }
    function die(message) {
        print message > "/dev/stderr"
        exit 1
    }
    function split_format(format_string, index_by_name,   count, i, fields) {
        delete index_by_name
        count = split(format_string, fields, ":")
        for (i = 1; i <= count; i++) {
            index_by_name[fields[i]] = i
        }
    }
    function require_format(index_by_name, field_name) {
        if (!(field_name in index_by_name)) {
            die("Missing GRIDSS FORMAT field required for normal VCF extraction: " field_name)
        }
    }
    function numeric_sample_field(sample_string, index_by_name, field_name,   fields, value) {
        split(sample_string, fields, ":")
        value = fields[index_by_name[field_name]]
        if (value == "" || value == "." || value ~ /,/) {
            return 0
        }
        return value + 0
    }
    function has_mateid(info_string,   count, i, fields, value) {
        if (info_string == "" || info_string == ".") {
            return 0
        }
        count = split(info_string, fields, ";")
        for (i = 1; i <= count; i++) {
            if (fields[i] == "MATEID") {
                return 1
            }
            if (fields[i] ~ /^MATEID=/) {
                value = substr(fields[i], 8)
                return value != "" && value != "."
            }
        }
        return 0
    }
    function safe_af(variant_support, ref_support, refpair_support,   depth) {
        depth = variant_support + ref_support + refpair_support
        if (depth > 0) {
            return variant_support / depth
        }
        return 0
    }
    /^##INFO=<ID=MATEID[,>]/ {
        has_mateid_header = 1
        print
        next
    }
    /^##/ {
        print
        next
    }
    /^#CHROM/ {
        if (NF < 11) {
            die("Expected at least two sample columns in GRIDSS VCF header")
        }
        print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10
        saw_header = 1
        next
    }
    {
        if (!saw_header) {
            die("Missing #CHROM header in GRIDSS VCF")
        }
        if (!has_mateid_header) {
            die("Missing INFO/MATEID header in GRIDSS VCF; cannot distinguish breakpoint and single-breakend records.")
        }
        if (NF < 11) {
            die("Expected at least two sample columns in GRIDSS VCF record")
        }

        split_format(\$9, format_index)
        require_format(format_index, "REF")
        require_format(format_index, "REFPAIR")
        require_format(format_index, "VF")
        require_format(format_index, "BVF")

        normal_ref = numeric_sample_field(\$10, format_index, "REF")
        normal_refpair = numeric_sample_field(\$10, format_index, "REFPAIR")
        tumour_ref = numeric_sample_field(\$11, format_index, "REF")
        tumour_refpair = numeric_sample_field(\$11, format_index, "REFPAIR")

        normal_vf = numeric_sample_field(\$10, format_index, "VF")
        tumour_vf = numeric_sample_field(\$11, format_index, "VF")
        normal_bvf = numeric_sample_field(\$10, format_index, "BVF")
        tumour_bvf = numeric_sample_field(\$11, format_index, "BVF")

        normal_bp_af = safe_af(normal_vf, normal_ref, normal_refpair)
        tumour_bp_af = safe_af(tumour_vf, tumour_ref, tumour_refpair)
        normal_be_af = safe_af(normal_bvf, normal_ref, normal_refpair)
        tumour_be_af = safe_af(tumour_bvf, tumour_ref, tumour_refpair)

        is_breakpoint = has_mateid(\$8)
        is_single_breakend = !is_breakpoint
        is_pass = \$7 == "PASS"

        keep_breakpoint = is_breakpoint && is_pass && normal_bp_af >= 0.10 && normal_vf >= 4 && tumour_bp_af <= (3 * normal_bp_af)
        keep_single_breakend = is_single_breakend && is_pass && normal_be_af >= 0.10 && normal_bvf >= 4 && tumour_be_af <= (3 * normal_be_af)

        if (keep_breakpoint || keep_single_breakend) {
            print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10
        }
    }
    END {
        if (!saw_header) {
            die("Missing #CHROM header in GRIDSS VCF")
        }
        if (!has_mateid_header) {
            die("Missing INFO/MATEID header in GRIDSS VCF; cannot distinguish breakpoint and single-breakend records.")
        }
    }' ${gridss_joint_vcf} > ${normal_vcf}
    """
}
