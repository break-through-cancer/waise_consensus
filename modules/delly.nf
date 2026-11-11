process DELLY {

    tag { "${task.process}_${task.attempt}_${id}" }

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/delly:1.0.3--h358d541_4'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        path reference_dir

    output:
        val "${id}", emit: id
        path "delly/${id}/", emit: outdir
        path "delly/${id}/*_delly.bcf", emit: bcf
        tuple val(id), path("delly/${id}/*_delly.bcf"), emit: bcf_t

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

    tag { "${task.process}_${task.attempt}_${id}" }

    label 'utility'

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
        val id
        path delly_bcf

    output:
        tuple val(id), path("delly/${id}/*_normal.vcf"), emit: vcf_n
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

process DELLY_SOMATIC_PREP {
    // Index the joint BCF and name the tumour and control samples for delly filter.
    // The tumour is the first sample in the BCF, as in DELLY_SPLIT.

    tag { id }

    label 'utility'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
        tuple val(id), path(delly_bcf, stageAs: 'input.bcf')

    output:
        tuple val(id), path('input.bcf'), path('input.bcf.csi'), path('samples.tsv'), emit: prepared

    script:
    """
    bcftools index -f input.bcf
    sids=`bcftools query -l input.bcf`
    tumor_id=`echo \$sids | cut -f1 -d' '`
    normal_id=`echo \$sids | cut -f2 -d' '`
    printf '%s\\ttumor\\n%s\\tcontrol\\n' "\$tumor_id" "\$normal_id" > samples.tsv
    """
}

process DELLY_FILTER {
    // Somatic filter from DELLY itself: alt fraction and coverage in the tumour, and no alt support in
    // the control. Extra options (for example -a or -p) go in --delly_filter_args.

    tag { id }

    label 'utility'

    container 'quay.io/biocontainers/delly:1.0.3--h358d541_4'

    input:
        tuple val(id), path(delly_bcf), path(delly_csi), path(samples)

    output:
        tuple val(id), path("${id}_delly_somatic.bcf"), emit: somatic

    script:
    """
    delly filter -f somatic ${params.delly_filter_args ?: ''} -o ${id}_delly_somatic.bcf -s ${samples} ${delly_bcf}
    """
}

process DELLY_TO_VCF {
    tag { id }

    label 'utility'

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.19--h8b25389_0'

    input:
        tuple val(id), path(somatic_bcf)

    output:
        tuple val(id), val('delly'), path("delly_filter/${id}/${id}_delly_somatic.vcf"), emit: vcf

    script:
    """
    mkdir -p delly_filter/${id}
    # Keep the tumour sample only; Jasmine needs every input VCF to have the same samples
    tumor_id=`bcftools query -l ${somatic_bcf} | head -n 1`
    bcftools view -s "\$tumor_id" ${somatic_bcf} -o delly_filter/${id}/${id}_delly_somatic.vcf
    """
}
