process DELLY {

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/delly:1.0.3--h358d541_4'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        path reference_dir

    output:
        val "${id}", emit: id
        path "delly/${id}/", emit: outdir
        path "delly/${id}/*_delly.bcf", emit: bcf

    script:
    def delly_outdir = "delly/${id}/"
    def reference_fasta = "${reference_dir}/reference.fa"
    
    """
    mkdir -p ${delly_outdir}
    reference_fasta_abs="\$(pwd)/${reference_fasta}"
    delly call \
    -g "\$reference_fasta_abs" \
    -o ${delly_outdir}/${id}_delly.bcf \
    "${tumor_bam}" \
    "${normal_bam}"
    """
}

process DELLY_SPLIT {

    label 'utility'

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
        val id
        path delly_bcf

    output:
        path "delly/${id}/*_normal.vcf", emit: vcf_n
        tuple val(id), path("delly/${id}/*_tumor.vcf"), emit: vcf_t

    script:
    def delly_outdir = "delly/${id}/"

    """
    mkdir -p ${delly_outdir} # For publishing
    sids=`bcftools query -l ${delly_bcf}`
    tumor_id=`echo \$sids | cut -f1 -d' '` # ID order is as given to DELLY_CALL
    normal_id=`echo \$sids | cut -f2 -d' '`

    # Split VCF into tumor and normal
    echo "Splitting joint output VCF"
    bcftools view -c1 -Ov -s \$normal_id -o ${delly_outdir}/${id}_normal.vcf ${delly_bcf}
    bcftools view -c1 -Ov -s \$tumor_id -o ${delly_outdir}/${id}_tumor.vcf ${delly_bcf}
    """
}
