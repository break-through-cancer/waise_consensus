process LUMPY {

    tag { "${task.process}_${task.attempt}_${id}" }

    publishDir "${params.publishDir}", mode: 'copy'
    
    container 'quay.io/biocontainers/smoove:0.2.8--h9ee0642_1'

    input:
        tuple val(id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        path reference_dir

    output:
        path "lumpy/${id}/", emit: outdir
        tuple val(id), path("lumpy/${id}/*_normal.vcf"), emit: vcf_n
        tuple val(id), path("lumpy/${id}/*_tumor.vcf"), emit: vcf_t

    script:
    def lumpy_outdir = "lumpy/${id}/"
    def reference_fasta = "${reference_dir}/reference.fa"

    """
    # The image's bcftools links libopenblas.so.0, which is absent; the image's libblas.so.3 satisfies it.
    # Smoove runs bcftools for --genotype. The 0.2.8--h9ee0642_0 build has the opposite fault (lumpy_filter
    # needs libcrypto.so.1.0.0), so this build stays.
    mkdir -p blas_shim
    ln -s /usr/local/lib/libblas.so.3 blas_shim/libopenblas.so.0
    export LD_LIBRARY_PATH="\$(pwd)/blas_shim"
    bcftools --version | head -1

    reference_fasta_abs="\$(pwd)/${reference_fasta}"
    smoove call \
        --name ${id} \
        --outdir ${lumpy_outdir} \
        -f "\$reference_fasta_abs" \
        -processes ${task.cpus} \
        --genotype \
        --removepr \
        --support 3 \
        "${normal_bam}" \
        "${tumor_bam}"

    # Get sample IDs from VCF header. With --genotype smoove writes the genotyped file.
    echo "Getting sample IDs from joint VCF"
    lumpy_file="${lumpy_outdir}/${id}-smoove.genotyped.vcf"
    gunzip -k \${lumpy_file}.gz
    sids=`vcfutils.pl listsam \$lumpy_file`
    normal_id=`echo \$sids | cut -f1 -d' '`
    tumor_id=`echo \$sids | cut -f2 -d' '`

    # Split the joint VCF into tumor and normal, keeping the sites where that sample has an alt
    # genotype (the same rule as bcftools view -c1 for DELLY and GRIDSS), read here with awk.
    echo "Splitting joint output VCF"
    for pair in "normal:\$normal_id" "tumor:\$tumor_id"; do
        role=\${pair%%:*}
        sample=\${pair#*:}
        vcfutils.pl subsam \$lumpy_file \$sample | awk 'BEGIN { FS = OFS = "\\t" }
            /^#/ { print; next }
            {
                n = split(\$9, fmt, ":")
                for (i = 1; i <= n; i++) if (fmt[i] == "GT") gt = i
                split(\$10, val, ":")
                if (val[gt] ~ /1/) print
            }' > ${lumpy_outdir}/${id}_\${role}.vcf
    done
    """

}
