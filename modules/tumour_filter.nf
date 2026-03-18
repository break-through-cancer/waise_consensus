
process TUMOUR_FILTER {
    // Filter tumour VCFs against caller-specific normal panels after dropping malformed primary coordinates.

    label 'utility'

    publishDir "${params.publishDir}/filtered", mode: 'copy'

    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'

    input:
        tuple val(id), path(vcf_file)
        path normal_panel
        val tool

    output:
        tuple val(id), path("${id}_${tool}_filtered.vcf"), emit: filtered_vcf

    script:
    '''
    set -euo pipefail

    tumour_valid_vcf="tumour_primary_valid.vcf"
    panel_valid_vcf="panel_primary_valid.vcf"

    sanitize_vcf() {
        local input_vcf="$1"
        local output_vcf="$2"

        : > "${output_vcf}"

        # Keep the VCF header, but drop records whose primary POS falls outside the contig bounds.
        VALID_VCF="${output_vcf}" awk '
        BEGIN {
            FS = OFS = "\\t"
            valid_vcf = ENVIRON["VALID_VCF"]
        }

        function store_contig_length(line, payload, parts) {
            payload = line
            sub(/^##contig=<ID=/, "", payload)
            split(payload, parts, /,length=|>/)
            if (parts[1] != "" && parts[2] ~ /^[0-9]+$/) {
                contig_length[parts[1]] = parts[2] + 0
            }
        }

        function primary_invalid(chrom, pos) {
            return (pos <= 0) || ((chrom in contig_length) && pos > contig_length[chrom])
        }

        /^##contig=<ID=/ {
            store_contig_length($0)
            print > valid_vcf
            next
        }

        /^#/ {
            print > valid_vcf
            next
        }

        {
            pos = $2 + 0

            # Drop invalid primary positions before bedtools sees the VCF body.
            if (!primary_invalid($1, pos)) {
                print > valid_vcf
            }
        }
        ' "${input_vcf}"
    }

    sanitize_vcf !{vcf_file} "${tumour_valid_vcf}"
    sanitize_vcf !{normal_panel} "${panel_valid_vcf}"

    # Rebuild the final VCF explicitly because bedtools subtract does not document header preservation.
    grep '^#' "${tumour_valid_vcf}" > !{id}_!{tool}_filtered.vcf

    if grep -qv '^#' "${panel_valid_vcf}"; then
        bedtools subtract -a "${tumour_valid_vcf}" -b "${panel_valid_vcf}" | awk '!/^#/' >> !{id}_!{tool}_filtered.vcf
    else
        awk '!/^#/' "${tumour_valid_vcf}" >> !{id}_!{tool}_filtered.vcf
    fi
    '''
}
