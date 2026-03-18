// Both prep paths normalize the staged FASTA name to reference/reference.fa.
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
        path reference_dir
        path gridss_blacklist_source
        path gridss_properties_source, name: 'gridss.properties.input'

    output:
        path 'gridss_blacklist.bed', emit: blacklist
        path 'gridss.properties', emit: properties

    script:
    """
    reference_index="${reference_dir}/reference.fa.fai"

    if [[ "${gridss_blacklist_source}" == *.gz ]]; then
        gzip -dc ${gridss_blacklist_source} > raw_gridss_blacklist.bed
    else
        cp ${gridss_blacklist_source} raw_gridss_blacklist.bed
    fi

    if [[ ! -s "\$reference_index" ]]; then
        echo "Missing staged reference index at \$reference_index" >&2
        exit 1
    fi

    if [[ -s raw_gridss_blacklist.bed ]]; then
        : > missing_blacklist_contigs.txt
        : > missing_blacklist_count.txt

        # GRIDSS expects blacklist contigs to match the staged reference dictionary exactly.
        # Human reference bundles are inconsistent about using chr-prefixed versus bare
        # chromosome names, so remap only simple aliases against the staged .fai first.
        # Any remaining non-matching contigs are skipped so custom blacklist files can
        # still contribute intervals on contigs present in the reference.
        awk 'BEGIN { OFS="\\t" }
        FNR == NR {
            ref[\$1] = 1
            next
        }
        function is_header() {
            return \$0 ~ /^[[:space:]]*$/ || \$0 ~ /^#/ || \$0 ~ /^track([[:space:]]|$)/ || \$0 ~ /^browser([[:space:]]|$)/
        }
        function map_contig(contig, core, candidate) {
            core = contig
            sub(/^chr/, "", core)

            if (contig in ref) return contig
            if (core in ref) return core

            candidate = "chr" core
            if (candidate in ref) return candidate

            if (core == "M" || core == "MT") {
                if ("chrM" in ref) return "chrM"
                if ("MT" in ref) return "MT"
                if ("M" in ref) return "M"
            }

            return contig
        }
        is_header() {
            print
            next
        }
        {
            \$1 = map_contig(\$1)
            if (!(\$1 in ref)) {
                missing[\$1] = 1
                missing_count += 1
                next
            }
            print
        }
        END {
            for (contig in missing) {
                print contig > "missing_blacklist_contigs.txt"
            }
            print missing_count > "missing_blacklist_count.txt"
        }' "\$reference_index" raw_gridss_blacklist.bed > gridss_blacklist.bed

        # Surface dropped contigs clearly so users can still spot custom blacklist and
        # reference mismatches without blocking GRIDSS on partial overlap cases.
        LC_ALL=C sort -o missing_blacklist_contigs.txt missing_blacklist_contigs.txt

        if [[ -s missing_blacklist_contigs.txt ]]; then
            missing_count="\$(cat missing_blacklist_count.txt)"
            echo "Skipping \$missing_count GRIDSS blacklist intervals whose contigs do not exist in the staged reference after normalization." >&2
            echo "Skipped contigs:" >&2
            head -n 10 missing_blacklist_contigs.txt >&2
            echo "Provide --gridss_blacklist with a BED built for this reference if you need those regions retained." >&2
        fi
    else
        : > gridss_blacklist.bed
    fi

    cp gridss.properties.input gridss.properties
    """
}
