import csv
import logging
import os
import runpy
import sys
import tempfile
import types
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


class FakeDataset:
    files = None
    samplesheet = None
    params = {}
    added = {}
    logger = logging.getLogger("fake")

    @classmethod
    def from_running(cls):
        return cls()

    def add_param(self, key, value):
        FakeDataset.added[key] = value


def run_preprocess(folder, files, samplesheet, params):
    FakeDataset.files = pd.DataFrame(files)
    FakeDataset.samplesheet = pd.DataFrame(samplesheet)
    FakeDataset.params = params
    FakeDataset.added = {}
    helpers = types.ModuleType("cirro.helpers.preprocess_dataset")
    helpers.PreprocessDataset = FakeDataset
    saved = {k: sys.modules.get(k) for k in ("cirro", "cirro.helpers", "cirro.helpers.preprocess_dataset")}
    sys.modules["cirro"] = types.ModuleType("cirro")
    sys.modules["cirro.helpers"] = types.ModuleType("cirro.helpers")
    sys.modules["cirro.helpers.preprocess_dataset"] = helpers
    cwd = os.getcwd()
    try:
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)
            runpy.run_path(str(ROOT / ".cirro" / folder / "preprocess.py"), run_name="__main__")
            with open("samplesheet.csv") as handle:
                return list(csv.DictReader(handle)), dict(FakeDataset.added)
    finally:
        os.chdir(cwd)
        for key, value in saved.items():
            if value is None:
                sys.modules.pop(key, None)
            else:
                sys.modules[key] = value


class CallPreprocess(unittest.TestCase):
    def test_pairs_and_fixed_entry(self):
        files = {
            "sample": ["T1", "T1", "T2", "T2", "N1", "N1"],
            "file": ["s3://b/T1.bam", "s3://b/T1.bam.bai", "s3://b/T2.cram", "s3://b/T2.cram.crai", "s3://b/N1.bam", "s3://b/N1.bam.bai"],
        }
        sheet = {"sample": ["T1", "T2", "N1"], "patient": ["P1", "P1", "P1"], "status": ["tumor", "tumor", "normal"]}
        files["sample"] = ["T1", "N1"] + files["sample"]
        files["file"] = ["s3://b/preprocessing/recal_table/T1/T1.table", "s3://b/preprocessing/recal_table/N1/N1.table"] + files["file"]
        rows, added = run_preprocess("call", files, sheet, {})
        self.assertEqual(added["entry"], "call")
        self.assertTrue(all(r["normal_bam"].endswith(".bam") and r["tumor_bam"].endswith((".bam", ".cram")) for r in rows))
        self.assertEqual([r["id"] for r in rows], ["T1", "T2"])
        self.assertEqual({r["normal_id"] for r in rows}, {"P1"})
        self.assertEqual(rows[1]["tumor_bai"], "s3://b/T2.cram.crai")


class CohortPreprocess(unittest.TestCase):
    def test_vcf_rows_and_fixed_entry(self):
        base = "s3://b/datasets/x/data/"
        files = {"file": [
            base + "manta/T1/results/variants/T1_somaticSV.vcf",
            base + "manta/T1/results/variants/T1_diploidSV.vcf",
            base + "delly/T1/T1_delly.bcf",
            base + "lumpy/T1/T1_tumor.vcf",
            base + "gridsspl/T1/T1_normal.vcf",
            base + "other/readme.txt",
        ]}
        sheet = {"sample": ["T1"], "patient": ["P1"]}
        rows, added = run_preprocess("cohort", files, sheet, {})
        self.assertEqual(added["entry"], "cohort")
        self.assertEqual(len(rows), 5)
        self.assertEqual({r["normal_id"] for r in rows}, {"P1"})
        self.assertIn(("delly", "bcf"), {(r["caller"], r["role"]) for r in rows})

    def test_no_vcfs_fails(self):
        with self.assertRaises(ValueError):
            run_preprocess("cohort", {"file": ["s3://b/x.txt"]}, {"sample": ["T1"]}, {})


if __name__ == "__main__":
    unittest.main()
