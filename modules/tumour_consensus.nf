
process TUMOUR_CONSENSUS_VCF {
    // Merge VCF files into a single VCF for downstream processing

    publishDir "${params.publishDir}/somatic", mode: 'copy'

    container 'quay.io/biocontainers/jasminesv:1.1.5--hdfd78af_0'

    input:
        tuple val(id), path(vcf_files)

    output:
        val id, emit: id
        path "${id}_consensus_jasmine_out.vcf", emit: consensus_vcf
        path "consensus_filelist.txt", emit: filelist

    script:
    """
    # Create a file list for Jasmine
    ls -1 ${vcf_files} > consensus_filelist.txt

    # Run Jasmine
    jasmine file_list=consensus_filelist.txt out_file=${id}_consensus_jasmine_out.vcf
    """
}

process TUMOUR_CONSENSUS_CALL {
    // Merge VCF files into a single VCF for downstream processing

    label 'utility'

    publishDir "${params.publishDir}/somatic", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

    input:
        val id
        path vcf
        val n_callers

    output:
        path "${id}_tumor_cons_min*.vcf", emit: level_vcfs
        path "${id}_consensus_levels.tsv", emit: level_summary
        val id, emit: id

    script:

    """
    awk '/^##INFO=<ID=SUPP/ { sub("Type=String", "Type=Integer"); } { print }' ${vcf} > tmp.vcf0
    printf 'min_support\\tn_records\\n' > ${id}_consensus_levels.tsv
    for k in \$(seq 2 ${n_callers}); do
        bcftools filter -i "INFO/SUPP >= \$k" tmp.vcf0 > ${id}_tumor_cons_min\${k}.vcf
        printf '%s\\t%s\\n' "\$k" "\$(grep -vc '^#' ${id}_tumor_cons_min\${k}.vcf || true)" >> ${id}_consensus_levels.tsv
    done
    """
}
