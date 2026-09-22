#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd

# 1. Get parameters from cirro pipeline call
ds = PreprocessDataset.from_running()
ds.logger.info("List of starting params")
ds.logger.info(ds.params)

# 1b. Validate SV caller selection / consensus threshold before doing any other
# work, so a bad combination fails fast here rather than after the Nextflow
# workflow has already started. Mirrors the equivalent checks in main.nf.
_caller_defaults = {
    "run_manta": True,
    "run_lumpy": True,
    "run_svaba": True,
    "run_delly": True,
    "run_gridss": True,
}
enabled_callers = {name: bool(ds.params.get(name, default)) for name, default in _caller_defaults.items()}
enabled_count = sum(enabled_callers.values())
ds.logger.info(f"Enabled callers: {sorted(n for n, on in enabled_callers.items() if on)} ({enabled_count} of {len(enabled_callers)})")

if enabled_count < 2:
    raise ValueError(
        f"Consensus calling requires at least 2 enabled SV callers, got {enabled_count} "
        f"({enabled_callers}). Enable more callers before launching this analysis."
    )

consensus_min_support = ds.params.get("consensus_min_support")
if consensus_min_support is not None:
    consensus_min_support = int(consensus_min_support)
    if not (2 <= consensus_min_support <= enabled_count):
        raise ValueError(
            f"consensus_min_support must be between 2 and the number of enabled callers "
            f"({enabled_count}); got {consensus_min_support}."
        )

ds.logger.info('checking ds.files')
files = ds.files
ds.logger.info(files.head())
ds.logger.info(files.columns)

# 2. Add samplesheet parameter and set equal to ds.samplesheet
ds.logger.info("Checking samplesheet parameter")
ds.logger.info(ds.samplesheet)
samples = ds.samplesheet.set_index("sample") if "sample" in ds.samplesheet.columns else pd.DataFrame()

# create samplesheet from ds.files
df = ds.files.copy()

# Some sources (e.g. sarek align output) provide both a 'recalibrated' and a
# 'markduplicates' BAM/BAI per sample. Without this, the bam/bai pivot below
# groups by filetype only and can pair a BAM with a BAI from the other
# bamType (e.g. recal.bam with md.bam.bai) since neither is unique per sample.
if "bamType" in df.columns:
    has_recal = df.groupby("sample")["bamType"].transform(lambda s: (s == "recalibrated").any())
    df = df[~((df["bamType"] == "markduplicates") & has_recal)]

# file type from path
df["is_bai"] = df["file"].str.endswith(".bai")
df["filetype"] = df["is_bai"].map({True: "bai", False: "bam"})

# case_id: prefer an explicit 'patient' column from the Cirro sample sheet so
# grouping doesn't depend on any particular sample-naming scheme. Fall back to
# the "<case>.S<n>..."/"<case>.PBMC" convention for datasets that rely on
# naming instead of sample metadata, then to the sample name itself (each
# sample is its own case) if neither yields a value.
name_case_id = df["sample"].str.extract(r"^(.*?)(?:\.S\d+.*|\.PBMC)$", expand=False)
meta_case_id = df["sample"].map(samples["patient"]) if "patient" in samples.columns else pd.Series(pd.NA, index=df.index, dtype=object)
df["case_id"] = meta_case_id.fillna(name_case_id).fillna(df["sample"]).astype(str)

# role: prefer an explicit 'status' column (tumor/normal, case-insensitive, or
# sarek-style 1/0). Fall back to the ".PBMC" naming convention, then default
# to "tumor" -- this pipeline's primary use case is SV calling on a subject BAM.
_status_map = {"tumor": "tumor", "1": "tumor", "normal": "normal", "0": "normal"}
meta_role = (
    df["sample"].map(samples["status"]).astype(str).str.lower().map(_status_map)
    if "status" in samples.columns else pd.Series(pd.NA, index=df.index, dtype=object)
)
name_role = df["sample"].str.contains(r"\.PBMC$", regex=True).map({True: "normal", False: "tumor"})
df["role"] = meta_role.fillna(name_role)

# one row per sample with bam/bai split out
per_sample = (
    df.pivot_table(
        index=["case_id", "sample", "role"],
        columns="filetype",
        values="file",
        aggfunc="first"
    )
    .reset_index()
)

tumor = (
    per_sample[per_sample["role"] == "tumor"]
    .rename(columns={
        "sample": "id",
        "bam": "tumor_bam",
        "bai": "tumor_bai",
    })[["case_id", "id", "tumor_bam", "tumor_bai"]]
)

normal = (
    per_sample[per_sample["role"] == "normal"]
    .rename(columns={
        "bam": "normal_bam",
        "bai": "normal_bai",
    })[["case_id", "normal_bam", "normal_bai"]]
)

samplesheet = (
    tumor.merge(normal, on="case_id", how="left")
         [["id", "tumor_bam", "tumor_bai", "normal_bam", "normal_bai"]]
         .sort_values("id")
         .reset_index(drop=True)
)

pd.set_option('display.max_columns', None)
print(samplesheet)


#%%
samplesheet.head(1).to_csv("samplesheet.csv", index=False)
ds.add_param("csv", "samplesheet.csv")

ds.logger.info(ds.params)
