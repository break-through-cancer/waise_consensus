process LUMPY {

    publishDir "${params.publishDir}", mode: 'copy'
    
    container 'quay.io/biocontainers/smoove:0.2.8--h9ee0642_1'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        path reference_dir

    output:
        path "lumpy/${id}/", emit: outdir
        path "lumpy/${id}/*_normal.vcf", emit: vcf_n
        tuple val(id), path("lumpy/${id}/*_tumor.vcf"), emit: vcf_t

    script:
    def lumpy_outdir = "lumpy/${id}/"
    def reference_fasta = "${reference_dir}/reference.fa"

    """
    reference_fasta_abs="\$(pwd)/${reference_fasta}"
    smoove call \
        --name ${id} \
        --outdir ${lumpy_outdir} \
        -f "\$reference_fasta_abs" \
        -processes ${task.cpus} \
        --removepr \
        --support 3 \
        "${normal_bam}" \
        "${tumor_bam}"

    # Get sample IDs from VCF header
    echo "Getting sample IDs from joint VCF"
    lumpy_file="${lumpy_outdir}/${id}-smoove.vcf"
    gunzip -k \${lumpy_file}.gz
    sids=`vcfutils.pl listsam \$lumpy_file`
    normal_id=`echo \$sids | cut -f1 -d' '`
    tumor_id=`echo \$sids | cut -f2 -d' '`

    # Split VCF into tumor and normal
    echo "Splitting joint output VCF"
    vcfutils.pl subsam \$lumpy_file \$normal_id > ${lumpy_outdir}/${id}_normal.vcf
    vcfutils.pl subsam \$lumpy_file \$tumor_id > ${lumpy_outdir}/${id}_tumor.vcf
    """

}
