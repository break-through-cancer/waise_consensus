# Consensus Structural Variant Caller (Waise et al., 2025)

This repository contains a Nextflow implementation of the consensus structural variant (SV) pipeline presented by Waise et al. (2025) _bioRxiv_.

![img: Consensus SV pipeline](assets/NASA_pipeline.png)

## Installation

Install Nextflow from [https://www.nextflow.io/](https://www.nextflow.io/).

Install [Docker](https://docs.docker.com/get-started/get-docker/) or [Singularity](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html) on your system.

Clone the pipeline from github to your local working directory:

```bash
git clone https://github.com/sarawaise/NASA-SV.git
cd NASA-SV
```

## Quickstart

The pipeline requires a `samplesheet.csv` file listing tumor-normal BAM file pairs and their index files:

```
id,tumor_bam,tumor_bai,normal_bam,normal_bai
nasa1,t1.bam,t1.bam.bai,n1.bam,n1.bam.bai
nasa2,t2.bam,t1.bam.bai,n2.bam,n1.bam.bai
nasa3,t3.bam,t1.bam.bai,n3.bam,n1.bam.bai
```

> [!NOTE]
> BAM index files are required for Nextflow to include them in the working directory.

## Entry Modes

`--entry` selects how much of the pipeline runs:

| Entry | What runs | Input sheet |
| --- | --- | --- |
| `full` (default) | The callers, then the cohort stage | BAM sheet (above) |
| `call` | The callers only, publishing each caller's raw tumor and normal VCFs | BAM sheet |
| `cohort` | The cohort stage only, on VCFs that already exist | Long-format sheet (below) |

Calling is independent for each pair: a pair's caller output does not depend on which other pairs are in the run. The normal panel is the only step that looks across pairs, so it is the part that can be rerun on its own with `--entry cohort`.

The BAM sheet accepts an optional `normal_id` column naming the normal sample. It defaults to the normal BAM file name, so pairs that share a normal BAM share a normal.

## Cohort Stage

The cohort stage builds a normal panel for each caller, filters each pair's tumor calls, and merges the callers into a consensus. It takes a long-format sheet with one row per VCF:

```
id,normal_id,caller,role,vcf
pair1,patientA,manta,tumor,pair1_somaticSV.vcf
pair1,patientA,manta,normal,pair1_diploidSV.vcf
pair1,patientA,delly,tumor,pair1_delly_tumor.vcf
pair1,patientA,delly,normal,pair1_delly_normal.vcf
pair2,patientA,manta,tumor,pair2_somaticSV.vcf
pair2,patientA,manta,normal,pair2_diploidSV.vcf
,,manta,panel,frozen_manta_panel.vcf
```

- `caller` is any label (letters, digits, and hyphens). The five built-in callers are not special; VCFs from other callers can be added as more `caller` values, and the consensus uses every caller present. The caller order in `SUPP_VEC` is alphabetical and is written to `somatic/<id>/<id>_consensus_callers.tsv`.
- `role` is `tumor` (the pair's somatic calls), `normal` (the pair's normal calls), `panel` (a frozen panel for one caller, used instead of building one), or `bcf` (the DELLY joint BCF, caller `delly` only; see below).
- `normal_id` is required on normal rows. A normal that appears in several pairs is counted once in the panel.
- A caller without normal rows and without a frozen panel is not filtered against a panel.
- At least two callers need tumor rows.

Inputs must use the same reference build and contig names. Each caller's VCFs should come from the same caller settings, since a panel built from one caller's normals only makes sense for that caller.

## Normal Panel and Matched Normal

Each caller's panel keeps the SVs found in at least `--panel_min_support` distinct normals (default 2). The panel is merged with Jasmine from one VCF per distinct normal, so several tumors from one patient do not make that patient's normal calls look recurrent. With fewer distinct normals than the threshold the panel is empty.

Normal records enter the panel and the subtraction only when the normal carries the SV (`--normal_alt_only true`). By default that means an alt genotype (`GT[0]`). The calling step runs `smoove call --genotype` for LUMPY and writes only the sites where each sample has an alt genotype, as it does for DELLY and GRIDSS, so the LUMPY normal VCF holds the SVs the normal carries. LUMPY VCFs from earlier runs have no genotype and list every site LUMPY genotyped, so subtracting them would remove every tumor call. For such a VCF the file is used as it is in the panel, matched-normal subtraction is skipped for that caller, and the task log says so. To decide with another measure, set `normal_filters` in a params file, for example:

```yaml
normal_filters:
  lumpy: 'FORMAT/SU[0]>=3'
```

The expression is a `bcftools view -i` expression and cannot contain single quotes. Normal VCFs should hold a single sample.

Two options adjust how callers take part:

- `--normal_filtered_callers` (default `[gridss]`) lists callers whose tumor calls the caller itself already filtered against the normal. GRIDSS runs `gridss_somatic_filter` with a panel of normals, and neither the cohort panel nor the subtraction removed any GRIDSS record in testing. These callers skip both steps, and no panel is built for them, which avoids merging the large GRIDSS normal VCFs.
- `--matched_normal_by_id` (default `[lumpy]`) lists callers whose tumor and normal VCFs come from one joint run and share record ids. For them the matched-normal subtraction removes the tumor record with the same id as an alt normal record, instead of removing every tumor record that overlaps a normal interval, which removed hundreds more.

By default each pair's own normal calls are also subtracted from its tumor calls (`--matched_normal_subtract true`). Normal calls that enter this step come from the joint tumor-normal run for that pair.

Tumor records are removed when they are fully covered by the panel and matched-normal intervals, as `bedtools subtract` does. `bedtools subtract` writes a record once for each fragment that remains, so identical lines are collapsed to a single copy. `--panel_subtract_any_overlap true` removes a record on any overlap instead (`bedtools subtract -A`).

## DELLY Somatic Filter

Waise et al. identified somatic DELLY calls with `delly filter` instead of subtracting the germline calls. With `--delly_filter true` (default) the cohort stage runs `delly filter -f somatic` on the DELLY joint BCF, a `bcf` row in the sheet, in place of that pair's DELLY tumor VCF. The tumor is the first sample in the BCF, as in the calling step. DELLY's defaults apply: a minimum alt fraction of 0.2 and coverage of 10 in the tumor, and no alt support in the control. Other options, such as `-a 0.1` or `-p` for PASS only, go in `--delly_filter_args`. The matched-normal subtraction is not applied to these calls, because the filter already tests the control. The panel still is. The BCF needs no index; one is built. With `--delly_filter false` the `bcf` row is ignored and the `tumor` row is used.

On one tumor-normal pair the DELLY filter kept 595 of 67,265 records, where the split tumor VCF had 27,910 and the matched-normal subtraction left 615. Only 66 records were in both sets, so the two filters are not interchangeable.

## Merge Inputs

Jasmine writes inconsistent records when its input VCFs have different numbers of samples, and `bcftools` then stops reading the merged VCF part way through. The somatic Manta and SvABA VCFs have two samples, the others one. Every VCF is therefore reduced to its first sample before the merge, and the consensus step stops with an error if `bcftools` cannot parse every record of the merged VCF.

## Tumour Filter Safeguards

The tumour filter drops malformed VCF records whose primary `POS` is outside the allowed contig bounds before running `bedtools`. This covers both `POS <= 0` and records whose `POS` exceeds a declared `##contig` length, preventing caller-specific edge cases from aborting tumour filtering.

## Reference Selection

The pipeline no longer assumes reference assets are stored under `projectDir`.

For the common human builds, set `--genome_build` and let the workflow derive the canonical reference FASTA, GRIDSS blacklist, GRIDSS properties, and the HMF GRIDSS panel-of-normals bundle used by `gridss_somatic_filter`:

```bash
nextflow run main.nf -profile singularity --csv samplesheet.csv --genome_build hg38
nextflow run main.nf -profile singularity --csv samplesheet.csv --genome_build hg19
```

Automatic GRIDSS PoN bundle extraction is only enabled for the exact `--genome_build hg38` and `--genome_build hg19` values.

If you use a non-default or nonhuman reference, set `--genome_build other` and provide the FASTA explicitly. The workflow will build the `.fai` and BWA index files if they are missing. GRIDSS will still run `gridss_somatic_filter` in this mode; if you do not provide a PoN, it falls back to matched-normal filtering without `--pondir`:

```bash
nextflow run main.nf -profile singularity \
  --csv samplesheet.csv \
  --genome_build other \
  --reference_fasta /path/to/genome.fa
```

If you have species-specific or custom GRIDSS PoN files, provide both explicitly:

```bash
nextflow run main.nf -profile singularity \
  --csv samplesheet.csv \
  --genome_build other \
  --reference_fasta /path/to/genome.fa \
  --gridss_sv_pon /path/to/gridss_pon_breakpoint.bedpe \
  --gridss_sgl_pon /path/to/gridss_pon_single_breakend.bed
```

You can override any derived human defaults with custom files:

```bash
nextflow run main.nf -profile singularity \
  --csv samplesheet.csv \
  --genome_build hg38 \
  --reference_fasta /path/to/custom_hg38.fa \
  --gridss_blacklist /path/to/custom_blacklist.bed \
  --gridss_sv_pon /path/to/custom_gridss_pon_breakpoint.bedpe \
  --gridss_sgl_pon /path/to/custom_gridss_pon_single_breakend.bed
```

> [!NOTE]
> The FASTA must still match the BAM contig names exactly. If your alignments were produced against a custom or differently named human reference, override the default FASTA and blacklist.
> For the built-in human defaults, the workflow also normalizes the ENCODE blacklist contig names to match the 1000 Genomes reference naming.
> Exact `--genome_build hg38` and `--genome_build hg19` still require a GRIDSS PoN and will auto-stage the HMF bundle unless you override it. `--genome_build other` is the only mode that permits running GRIDSS somatic filtering without a PoN.

If you want to override the build-derived HMF bundle itself instead of providing extracted PoN files, set `--gridss_resource_bundle` to a local or remote `.tar.gz`/gzip-compressed tar archive containing the expected HMF resource layout.

An archive file (10GB) with test BAM files and an indexed hg38 reference is [available on Zenodo (ID:15226469)](https://zenodo.org/records/15226469). To test the pipeline with the default online hg38 assets, extract the data archive into the repository root and run:

```bash
nextflow run main.nf -profile singularity,test
```

To reuse the original repository-relative reference defaults from the extracted bundle, add the compatibility profile:

```bash
nextflow run main.nf -profile singularity,test,bundle_paths
```

## Credits

The Genomics England implementation used for Waise et al., 2025 was written by Alex Cornish and Sara Waise. The Nextflow implementation provided in this repository was written by Nana Mensah and Sara Waise.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
