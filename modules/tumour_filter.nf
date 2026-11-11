process TUMOUR_FILTER {
    // Filter tumour VCFs against a caller-specific normal panel and the pair's own normal calls
    // after dropping malformed primary coordinates.

    tag { "${id}_${tool}" }

    label 'utility'

    publishDir "${params.publishDir}", mode: 'copy'

    container 'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2'

    input:
        tuple val(id), val(tool), path(vcf_file, stageAs: 'tumour_input.vcf'), path(normal_panel, stageAs: 'panel_input.vcf'), val(use_panel), path(matched_normal, stageAs: 'matched_normal_input.vcf'), val(use_matched)

    output:
        tuple val(id), val(tool), path("filtered/${id}/${id}_${tool}_filtered.vcf"), emit: filtered_vcf

    script:
    def subtract_flag = params.panel_subtract_any_overlap.toString().toBoolean() ? '-A' : ''
    """
    set -euo pipefail

    mkdir -p filtered/${id}

    tumour_valid_vcf="tumour_primary_valid.vcf"
    filtered_vcf="filtered/${id}/${id}_${tool}_filtered.vcf"

    sanitize_vcf() {
        local input_vcf="\$1"
        local output_vcf="\$2"

        : > "\${output_vcf}"

        # Keep the VCF header, but drop records whose primary POS falls outside the contig bounds.
        VALID_VCF="\${output_vcf}" awk '
        BEGIN {
            FS = OFS = "\\t"
            valid_vcf = ENVIRON["VALID_VCF"]
        }

        function store_contig_length(line, payload, parts) {
            payload = line
            sub(/^##contig=<ID=/, "", payload)
            split(payload, parts, /,length=|>/)
            if (parts[1] != "" && parts[2] ~ /^[0-9]+\$/) {
                contig_length[parts[1]] = parts[2] + 0
            }
        }

        function primary_invalid(chrom, pos) {
            return (pos <= 0) || ((chrom in contig_length) && pos > contig_length[chrom])
        }

        /^##contig=<ID=/ {
            store_contig_length(\$0)
            print > valid_vcf
            next
        }

        /^#/ {
            print > valid_vcf
            next
        }

        {
            pos = \$2 + 0

            # Drop invalid primary positions before bedtools sees the VCF body.
            if (!primary_invalid(\$1, pos)) {
                print > valid_vcf
            }
        }
        ' "\${input_vcf}"
    }

    sanitize_vcf tumour_input.vcf "\${tumour_valid_vcf}"

    # Records to subtract: the normal panel plus the pair's own normal calls
    : > exclude_body.txt

    if [ "${use_panel}" = "true" ]; then
        sanitize_vcf panel_input.vcf panel_primary_valid.vcf
        awk '!/^#/' panel_primary_valid.vcf >> exclude_body.txt
    fi

    if [ "${use_matched}" = "true" ]; then
        sanitize_vcf matched_normal_input.vcf matched_primary_valid.vcf
        awk '!/^#/' matched_primary_valid.vcf >> exclude_body.txt
    fi

    # Rebuild the final VCF explicitly because bedtools subtract does not document header preservation.
    grep '^#' "\${tumour_valid_vcf}" > "\${filtered_vcf}"

    if [ -s exclude_body.txt ]; then
        grep '^#' "\${tumour_valid_vcf}" > exclude.vcf
        cat exclude_body.txt >> exclude.vcf

        # bedtools subtract writes a record once per fragment left after removing an overlap, so
        # identical lines are collapsed to keep one copy of each record.
        bedtools subtract ${subtract_flag} -a "\${tumour_valid_vcf}" -b exclude.vcf | awk '!/^#/ && !seen[\$0]++' >> "\${filtered_vcf}"
    else
        awk '!/^#/' "\${tumour_valid_vcf}" >> "\${filtered_vcf}"
    fi
    """
}

process NORMAL_ALT_ONLY {
    // Keep the normal records where the normal carries the SV, so the panel and the matched-normal
    // subtraction only see calls with evidence in the normal. Callers differ in how they record this:
    //  - with an expression for the caller (--normal_filters), that bcftools expression is used;
    //  - otherwise records with an alt genotype (GT[0]) are kept when the normal VCF has genotype calls;
    //  - when the VCF has no usable genotype (LUMPY and GRIDSS leave it empty), the file passes through
    //    unchanged and is marked unusable, so its calls are not subtracted from the tumour calls.

    tag { "${id}_${caller}" }

    label 'utility'

    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

    input:
        tuple val(id), val(nid), val(caller), path(vcf, stageAs: 'input.vcf'), val(expr)

    output:
        tuple val(id), val(nid), val(caller), path("${id}_${caller}_normal_alt.vcf"), path("${id}_${caller}_normal_usable.txt"), emit: vcf

    script:
    """
    out=${id}_${caller}_normal_alt.vcf
    flag=${id}_${caller}_normal_usable.txt

    if [ -n '${expr}' ]; then
        bcftools view -i '${expr}' input.vcf -o \$out
        echo true > \$flag
    elif bcftools view -h input.vcf | grep -q '^##FORMAT=<ID=GT,'; then
        n_alt=\$(bcftools view -H -i 'GT[0]="alt"' input.vcf | wc -l)
        n_all=\$(bcftools view -H input.vcf | wc -l)
        if [ "\$n_alt" -gt 0 ] || [ "\$n_all" -eq 0 ]; then
            bcftools view -i 'GT[0]="alt"' input.vcf -o \$out
            echo true > \$flag
        else
            echo "No genotype calls in the ${caller} normal VCF for ${id}; matched-normal subtraction is skipped for it"
            cp input.vcf \$out
            echo false > \$flag
        fi
    else
        echo "No GT field in the ${caller} normal VCF for ${id}; matched-normal subtraction is skipped for it"
        cp input.vcf \$out
        echo false > \$flag
    fi
    """
}
