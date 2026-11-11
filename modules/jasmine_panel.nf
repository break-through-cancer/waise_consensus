// Process VCF files from normal samples with Jasmine.
// For each caller, PANEL_MERGE and PANEL_FILTER create a VCF with SVs observed in at least
// `min_support` distinct normals. Each distinct normal is passed in once.


process PANEL_MERGE {
    // Merge one VCF per distinct normal into a single VCF

    tag { tool }

    publishDir "${params.publishDir}/jasmine", mode: 'copy'

    container 'quay.io/biocontainers/jasminesv:1.1.5--hdfd78af_0'

    input:
        tuple val(tool), path(vcf_files), val(min_support)

    output:
        tuple val(tool), path("${tool}_jasmine_out.vcf"), val(min_support), emit: merged
        path "${tool}_jasmine_filelist.txt", emit: filelist

    script:
    """
    printf '%s\\n' ${vcf_files} > ${tool}_jasmine_filelist.txt

    # Fewer normals than the support threshold cannot produce a panel record, so skip the merge
    if [ "\$(wc -l < ${tool}_jasmine_filelist.txt)" -lt ${min_support} ]; then
        grep '^#' "\$(head -n 1 ${tool}_jasmine_filelist.txt)" > ${tool}_jasmine_out.vcf
    else
        jasmine file_list=${tool}_jasmine_filelist.txt out_file=${tool}_jasmine_out.vcf
    fi
    """
}

process PANEL_FILTER {
    // Keep the merged SVs supported by at least min_support normals

    tag { tool }

    label 'utility'

    publishDir "${params.publishDir}/panel", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

    input:
        tuple val(tool), path(vcf), val(min_support)

    output:
        tuple val(tool), path("${tool}_jasmine_cons.vcf"), emit: panel

    script:
    """
    if grep -q '^##INFO=<ID=SUPP' ${vcf}; then
        awk '/^##INFO=<ID=SUPP/ { sub("Type=String", "Type=Integer"); } { print }' ${vcf} > tmp.vcf0
        bcftools filter -i "INFO/SUPP >= ${min_support}" tmp.vcf0 > ${tool}_jasmine_cons.vcf
    else
        cp ${vcf} ${tool}_jasmine_cons.vcf
    fi
    """
}
