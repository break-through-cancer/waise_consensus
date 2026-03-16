process GRIDSS {

    publishDir "${params.publishDir}", mode: 'copy'

    container 'gridss/gridss:2.13.2'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        path reference_dir
        path gridss_blacklist
        path gridss_properties

    output:
        val "${id}", emit: id
        path "gridsspl/${id}/", emit: outdir
        path "gridsspl/${id}/${id}_normal.vcf", emit: vcf_n
        tuple val(id), path("gridsspl/${id}/${id}_tumor.vcf"), emit: vcf_t

    script:
    def gridsspl_outdir = "gridsspl/${id}"
    def vcf_pre = "gridsspl/${id}/${id}"
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

    # Run gridss
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

    # Split BAMs
    sids=`bcftools query -l ${gridsspl_outdir}/gridss.vcf`
    normal_id=`echo \$sids | cut -f1 -d' '` # ID order is as given to GRIDSSPL
    tumor_id=`echo \$sids | cut -f2 -d' '`

    # Split VCF into tumor and normal
    echo "Splitting joint output VCF"
    bcftools view -c1 -Ov -s \$normal_id -o ${vcf_pre}_normal.vcf ${gridsspl_outdir}/gridss.vcf
    bcftools view -c1 -Ov -s \$tumor_id -o ${vcf_pre}_tumor.vcf ${gridsspl_outdir}/gridss.vcf
    """
}
