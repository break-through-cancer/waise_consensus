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

## Tumour Filter Safeguards

The tumour-panel subtraction step now drops malformed VCF records whose primary `POS` is outside the allowed contig bounds before running `bedtools`. This covers both `POS <= 0` and records whose `POS` exceeds a declared `##contig` length, preventing caller-specific edge cases from aborting tumour filtering.

## Reference Selection

The pipeline no longer assumes reference assets are stored under `projectDir`.

For the common human builds, set `--genome_build` and let the workflow derive the canonical reference FASTA, GRIDSS blacklist, and GRIDSS properties:

```bash
nextflow run main.nf -profile singularity --csv samplesheet.csv --genome_build hg38
nextflow run main.nf -profile singularity --csv samplesheet.csv --genome_build hg19
```

If you use a non-default or nonhuman reference, set `--genome_build other` and provide the FASTA explicitly. The workflow will build the `.fai` and BWA index files if they are missing:

```bash
nextflow run main.nf -profile singularity \
  --csv samplesheet.csv \
  --genome_build other \
  --reference_fasta /path/to/genome.fa
```

You can override any derived human defaults with custom files:

```bash
nextflow run main.nf -profile singularity \
  --csv samplesheet.csv \
  --genome_build hg38 \
  --reference_fasta /path/to/custom_hg38.fa \
  --gridss_blacklist /path/to/custom_blacklist.bed
```

> [!NOTE]
> The FASTA must still match the BAM contig names exactly. If your alignments were produced against a custom or differently named human reference, override the default FASTA and blacklist.
> For the built-in human defaults, the workflow also normalizes the ENCODE blacklist contig names to match the 1000 Genomes reference naming.

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
