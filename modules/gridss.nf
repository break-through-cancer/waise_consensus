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
    def gridsspl_outdir = "gridsspl/${id}/"
    def vcf_pre = "gridsspl/${id}/${id}"
    def reference_fasta = "${reference_dir}/reference.fa"
    """
    blacklist_arg=""
    if [[ -s ${gridss_blacklist} ]]; then
        blacklist_arg="-b ${gridss_blacklist}"
    fi

    properties_arg=""
    if [[ -s ${gridss_properties} ]]; then
        properties_arg="-c ${gridss_properties}"
    fi

    # Run gridss
    gridss \
            -o ${gridsspl_outdir}/gridss.vcf \
            -r ${reference_fasta} \
            --threads $params.threads \
            --jvmheap $params.jvmheap \
            \$blacklist_arg \
            \$properties_arg \
            $normal_bam $tumor_bam

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
