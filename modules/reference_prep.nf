process PREP_REFERENCE {

    container 'gridss/gridss:2.13.2'

    input:
        path reference_fasta_source

    output:
        path 'reference', emit: ref_dir

    script:
    """
    mkdir -p reference

    if [[ "${reference_fasta_source}" == *.gz ]]; then
        gzip -dc ${reference_fasta_source} > reference/reference.fa
    else
        cp ${reference_fasta_source} reference/reference.fa
    fi

    samtools faidx reference/reference.fa
    bwa index reference/reference.fa
    """
}

process PREP_REFERENCE_WITH_INDEX {

    container 'gridss/gridss:2.13.2'

    input:
        path reference_fasta_source
        path reference_index_source

    output:
        path 'reference', emit: ref_dir

    script:
    """
    mkdir -p reference

    if [[ "${reference_fasta_source}" == *.gz ]]; then
        gzip -dc ${reference_fasta_source} > reference/reference.fa
        samtools faidx reference/reference.fa
    else
        cp ${reference_fasta_source} reference/reference.fa
        cp ${reference_index_source} reference/reference.fa.fai
    fi

    bwa index reference/reference.fa
    """
}

process PREP_GRIDSS_ASSETS {

    label 'utility'

    container 'gridss/gridss:2.13.2'

    input:
        path gridss_properties_source, name: 'gridss.properties.input'

    output:
        path 'gridss_blacklist.bed', emit: blacklist
        path 'gridss.properties', emit: properties

    script:
    """
    : > gridss_blacklist.bed
    cp gridss.properties.input gridss.properties
    """
}

process PREP_GRIDSS_ASSETS_WITH_BLACKLIST {

    label 'utility'

    container 'gridss/gridss:2.13.2'

    input:
        path gridss_blacklist_source
        path gridss_properties_source, name: 'gridss.properties.input'
        val normalize_default_blacklist

    output:
        path 'gridss_blacklist.bed', emit: blacklist
        path 'gridss.properties', emit: properties

    script:
    """
    if [[ "${gridss_blacklist_source}" == *.gz ]]; then
        gzip -dc ${gridss_blacklist_source} > raw_gridss_blacklist.bed
    else
        cp ${gridss_blacklist_source} raw_gridss_blacklist.bed
    fi

    if [[ "${normalize_default_blacklist}" == "true" && -s raw_gridss_blacklist.bed ]]; then
        awk 'BEGIN { OFS="\\t" } { sub(/^chr/, "", \$1); if (\$1 == "M") \$1 = "MT"; print }' raw_gridss_blacklist.bed > gridss_blacklist.bed
    else
        mv raw_gridss_blacklist.bed gridss_blacklist.bed
    fi

    cp gridss.properties.input gridss.properties
    """
}
