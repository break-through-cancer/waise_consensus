process MANTA {
    
    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/manta:1.6.0--py27_0'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        path reference_dir

    output:
        path "manta/${id}/", emit: outdir
        path "manta/${id}/results/variants/${id}*diploidSV.vcf", emit: vcf_n
        tuple val(id), path("manta/${id}/results/variants/${id}*somaticSV.vcf"), emit: vcf_t

    script:
    def manta_outdir = "manta/${id}/"
    def reference_fasta = "${reference_dir}/reference.fa"
    def manta_memory_gb = task.memory.toGiga().intValue()
    """
    configManta.py \
    --tumorBam ${tumor_bam} \
    --normalBam ${normal_bam} \
    --reference ${reference_fasta} \
    --runDir ${manta_outdir}

    # Run manta
    ${manta_outdir}/runWorkflow.py -m local -j ${task.cpus} -g ${manta_memory_gb}

    # Unzip for jasmine compatibility
    gunzip -c ${manta_outdir}/results/variants/somaticSV.vcf.gz > ${manta_outdir}/results/variants/${id}_somaticSV.vcf
    gunzip -c ${manta_outdir}/results/variants/diploidSV.vcf.gz > ${manta_outdir}/results/variants/${id}_diploidSV.vcf
    """

}
