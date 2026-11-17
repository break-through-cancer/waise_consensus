#!/usr/bin/env python3
"""Preprocess for the SV cohort process: reads the caller VCFs published by the SV calling process
in one or more datasets and builds the long-format cohort sheet (--entry cohort)."""

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import re

# 1. Get parameters from cirro pipeline call
ds = PreprocessDataset.from_running()
ds.logger.info("List of starting params")
ds.logger.info(ds.params)

ds.add_param("entry", "cohort")

ds.logger.info("Checking samplesheet parameter")
ds.logger.info(ds.samplesheet)
samples = ds.samplesheet.set_index("sample") if "sample" in ds.samplesheet.columns else pd.DataFrame()

# Caller and role of each VCF this pipeline publishes, by path. The cohort entry takes these
# straight from the datasets the callers wrote, so no sheet has to be written by hand.
_VCF_PATTERNS = [
    ("manta", "tumor", re.compile(r"/manta/(?P<id>[^/]+)/results/variants/[^/]*_somaticSV\.vcf$")),
    ("manta", "normal", re.compile(r"/manta/(?P<id>[^/]+)/results/variants/[^/]*_diploidSV\.vcf$")),
    ("svaba", "tumor", re.compile(r"/svaba/(?P<id>[^/]+)/[^/]*svaba\.unfiltered\.somatic\.sv\.vcf$")),
    ("svaba", "normal", re.compile(r"/svaba/(?P<id>[^/]+)/[^/]*svaba\.germline\.sv\.vcf$")),
    ("delly", "bcf", re.compile(r"/delly/(?P<id>[^/]+)/[^/]*_delly\.bcf$")),
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



# One row per VCF. The normal sample is identified by the patient, taken from the sample
# metadata when present, so pairs from the same patient share one normal in the panel.
df = ds.files.copy()
patient = samples["patient"].astype(str).to_dict() if "patient" in samples.columns else {}
rows = vcf_rows([str(p) for p in df["file"]], patient)
if not rows:
    raise ValueError("No caller VCFs recognized among the input files. Select datasets that contain the SV calling process outputs.")
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

ds.logger.info(ds.params)
