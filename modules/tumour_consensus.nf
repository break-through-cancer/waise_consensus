process NORMALIZE_VCF {
    // Jasmine groups records by chromosome, type, and strand; LUMPY is the only caller that writes
    // STRANDS on non-breakend records, so without this its calls never merge with the other callers.
    // The field is removed from every caller's copy that enters the consensus merge.

    tag { "${id}_${caller}" }

    label 'utility'

    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

    input:
        tuple val(id), val(caller), path(vcf, stageAs: 'input.vcf')

    output:
        tuple val(id), val(caller), path("${id}_${caller}_consensus_input.vcf"), emit: vcf

    script:
    """
    bcftools annotate -x INFO/STRANDS input.vcf -o ${id}_${caller}_consensus_input.vcf
    """
}

process TUMOUR_CONSENSUS_VCF {
    // Merge VCF files into a single VCF for downstream processing. The files are listed in the
    // order of `callers`, which is the order of the positions in SUPP_VEC.

    tag { id }

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/jasminesv:1.1.5--hdfd78af_0'

    input:
        tuple val(id), val(callers), path(vcf_files)

    output:
        val id, emit: id
        val "${callers.size()}", emit: n_callers
        path "somatic/${id}/${id}_consensus_jasmine_out.vcf", emit: consensus_vcf
        path "somatic/${id}/consensus_filelist.txt", emit: filelist
        path "somatic/${id}/${id}_consensus_callers.tsv", emit: caller_order

    script:
    """
    mkdir -p somatic/${id}

    # Create a file list for Jasmine, in caller order
    printf '%s\\n' ${vcf_files} > somatic/${id}/consensus_filelist.txt

    # Record which caller each position of SUPP_VEC refers to
    printf 'position\\tcaller\\n' > somatic/${id}/${id}_consensus_callers.tsv
    i=1
    for c in ${callers.join(' ')}; do
        printf '%s\\t%s\\n' "\$i" "\$c" >> somatic/${id}/${id}_consensus_callers.tsv
        i=\$((i + 1))
    done

    # Run Jasmine
    jasmine file_list=somatic/${id}/consensus_filelist.txt out_file=somatic/${id}/${id}_consensus_jasmine_out.vcf
    """
}

process TUMOUR_CONSENSUS_CALL {
    // Write one consensus VCF per support level for downstream processing

    tag { id }

    label 'utility'

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

    input:
        val id
        path vcf
        val n_callers

    output:
        path "somatic/${id}/${id}_tumor_cons_min*.vcf", emit: level_vcfs
        path "somatic/${id}/${id}_consensus_levels.tsv", emit: level_summary
        val id, emit: id

    script:

    """
    mkdir -p somatic/${id}
    awk '/^##INFO=<ID=SUPP/ { sub("Type=String", "Type=Integer"); } { print }' ${vcf} > tmp.vcf0
    printf 'min_support\\tn_records\\n' > somatic/${id}/${id}_consensus_levels.tsv
    for k in \$(seq 2 ${n_callers}); do
        bcftools filter -i "INFO/SUPP >= \$k" tmp.vcf0 > somatic/${id}/${id}_tumor_cons_min\${k}.vcf
        printf '%s\\t%s\\n' "\$k" "\$(grep -vc '^#' somatic/${id}/${id}_tumor_cons_min\${k}.vcf || true)" >> somatic/${id}/${id}_consensus_levels.tsv
    done
    """
}
