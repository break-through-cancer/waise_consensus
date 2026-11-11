#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd

# 1. Get parameters from cirro pipeline call
ds = PreprocessDataset.from_running()
ds.logger.info("List of starting params")
ds.logger.info(ds.params)

entry = str(ds.params.get("entry", "full")).lower()
if entry not in ("full", "call", "cohort"):
    raise ValueError(f"Unsupported entry '{entry}'. Use full, call, or cohort.")

# 1b. Validate SV caller selection / consensus threshold before doing any other
# work, so a bad combination fails fast here rather than after the Nextflow
# workflow has already started. Mirrors the equivalent checks in main.nf.
import re


# Caller and role of each VCF this pipeline publishes, by path. The cohort entry takes these
# straight from the datasets the callers wrote, so no sheet has to be written by hand.
_VCF_PATTERNS = [
    ("manta", "tumor", re.compile(r"/manta/(?P<id>[^/]+)/results/variants/[^/]*_somaticSV\.vcf$")),
    ("manta", "normal", re.compile(r"/manta/(?P<id>[^/]+)/results/variants/[^/]*_diploidSV\.vcf$")),
    ("svaba", "tumor", re.compile(r"/svaba/(?P<id>[^/]+)/[^/]*svaba\.unfiltered\.somatic\.sv\.vcf$")),
    ("svaba", "normal", re.compile(r"/svaba/(?P<id>[^/]+)/[^/]*svaba\.germline\.sv\.vcf$")),
    ("delly", "tumor", re.compile(r"/delly/(?P<id>[^/]+)/[^/]*_tumor\.vcf$")),
    ("delly", "normal", re.compile(r"/delly/(?P<id>[^/]+)/[^/]*_normal\.vcf$")),
    ("lumpy", "tumor", re.compile(r"/lumpy/(?P<id>[^/]+)/[^/]*_tumor\.vcf$")),
    ("lumpy", "normal", re.compile(r"/lumpy/(?P<id>[^/]+)/[^/]*_normal\.vcf$")),
    ("gridss", "tumor", re.compile(r"/gridsspl/(?P<id>[^/]+)/[^/]*_tumor\.vcf$")),
    ("gridss", "normal", re.compile(r"/gridsspl/(?P<id>[^/]+)/[^/]*_normal\.vcf$")),
]


def vcf_rows(paths, normal_ids):
    """Long-format rows (id, normal_id, caller, role, vcf) from published caller VCF paths.

    normal_ids maps a pair id to the identifier of its normal sample. Paths that match no
    pattern are ignored.
    """
    rows = []
    for path in paths:
        for caller, role, pattern in _VCF_PATTERNS:
            match = pattern.search(path)
            if match:
                pair_id = match.group("id")
                rows.append({
                    "id": pair_id,
                    "normal_id": normal_ids.get(pair_id, pair_id),
                    "caller": caller,
                    "role": role,
                    "vcf": path,
                })
                break
    return rows


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

if entry == "full" and enabled_count < 2:
    raise ValueError(
        f"Consensus calling requires at least 2 enabled SV callers, got {enabled_count} "
        f"({enabled_callers}). Enable more callers before launching this analysis."
    )

ds.logger.info('checking ds.files')
files = ds.files
ds.logger.info(files.head())
ds.logger.info(files.columns)

# 2. Add samplesheet parameter and set equal to ds.samplesheet
ds.logger.info("Checking samplesheet parameter")
ds.logger.info(ds.samplesheet)
samples = ds.samplesheet.set_index("sample") if "sample" in ds.samplesheet.columns else pd.DataFrame()

if entry == "cohort":
    # One row per VCF. The normal sample is identified by the patient, taken from the sample
    # metadata when present, so pairs from the same patient share one normal in the panel.
    df = ds.files.copy()
    patient = samples["patient"].astype(str).to_dict() if "patient" in samples.columns else {}
    rows = vcf_rows([str(p) for p in df["file"]], patient)
    if not rows:
        raise ValueError("No caller VCFs recognized among the input files. Select datasets that contain this pipeline's caller outputs.")
    missing = sorted({r["id"] for r in rows if r["id"] not in patient})
    if missing:
        ds.logger.warning(
            "No patient in the sample metadata for: %s. Their normal is identified by the pair id, "
            "so a normal shared with another pair would be counted once per pair." % ", ".join(missing)
        )
    sheet = pd.DataFrame(rows, columns=["id", "normal_id", "caller", "role", "vcf"]).sort_values(["id", "caller", "role"]).reset_index(drop=True)
    print(sheet)
    sheet.to_csv("samplesheet.csv", index=False)
    ds.add_param("csv", "samplesheet.csv")
else:
    # create samplesheet from ds.files
    df = ds.files.copy()

    # Some sources (e.g. sarek align output) provide both a 'recalibrated' and a
    # 'markduplicates' BAM/BAI per sample. Without this, the bam/bai pivot below
    # groups by filetype only and can pair a BAM with a BAI from the other
    # bamType (e.g. recal.bam with md.bam.bai) since neither is unique per sample.
    if "bamType" in df.columns:
        has_recal = df.groupby("sample")["bamType"].transform(lambda s: (s == "recalibrated").any())
        df = df[~((df["bamType"] == "markduplicates") & has_recal)]

    # file type from path: BAM+BAI or CRAM+CRAI. The "bam"/"bai" labels name the alignment and
    # its index; the callers read both formats given the reference FASTA.
    df["is_bai"] = df["file"].str.contains(r"\.(?:bai|crai)$", regex=True)
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
             .assign(normal_id=lambda d: d["case_id"])
             [["id", "tumor_bam", "tumor_bai", "normal_bam", "normal_bai", "normal_id"]]
             .sort_values("id")
             .reset_index(drop=True)
    )

    pd.set_option('display.max_columns', None)
    print(samplesheet)


    #%%
    samplesheet.to_csv("samplesheet.csv", index=False)
    ds.add_param("csv", "samplesheet.csv")

ds.logger.info(ds.params)
