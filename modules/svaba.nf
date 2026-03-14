process SVABA {

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/svaba:1.1.0--h468198e_3'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        path reference_dir

    output:
        path "svaba/${id}/", emit: outdir
        path "svaba/${id}/*svaba.germline.sv.vcf", emit: vcf_n
        tuple val(id), path("svaba/${id}/*svaba.unfiltered.somatic.sv.vcf"), emit: vcf_t

    script:
    def svaba_outdir = "svaba/${id}/"
    """
    mkdir -p ${svaba_outdir}
    tbam=\$(pwd)/${tumor_bam}
    nbam=\$(pwd)/${normal_bam}
    ref=\$(pwd)/${reference_dir}/reference.fa

    cd ${svaba_outdir}
    svaba run \
        -t \$tbam \
        -n \$nbam \
        -a ${id} \
        --threads ${task.cpus} \
        --reference-genome \$ref
    """
}
