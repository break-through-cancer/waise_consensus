#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import numpy as np
import re

# 1. Get parameters from cirro pipeline call
ds = PreprocessDataset.from_running()
ds.logger.info("List of starting params")
ds.logger.info(ds.params)

ds.logger.info('checking ds.files')
files = ds.files
ds.logger.info(files.head())
ds.logger.info(files.columns)

# 2. Add samplesheet parameter and set equal to ds.samplesheet
ds.logger.info("Checking samplesheet parameter")
ds.logger.info(ds.samplesheet)

# create samplesheet from ds.files
df = ds.files.copy()

# file type from path
df["is_bai"] = df["file"].str.endswith(".bai")
df["filetype"] = df["is_bai"].map({True: "bai", False: "bam"})

# case id = everything before .S<number>... or .PBMC
df["case_id"] = df["sample"].str.extract(r"^(.*?)(?:\.S\d+.*|\.PBMC)$", expand=False)

# role
df["role"] = df["sample"].str.contains(r"\.PBMC$", regex=True).map({True: "normal", False: "tumor"})

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