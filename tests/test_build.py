"""Offline tests for scripts/build.py and scripts/validate.py.

No network: fetch_cpic / fetch_pharmgkb / check_alleles are replaced with fixtures that mirror the
shapes of the real CPIC API responses and the PharmGKB clinicalAnnotations.zip download.
Run:  python -I -m unittest discover -s tests
"""
from __future__ import annotations

import io
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import build  # noqa: E402
import validate  # noqa: E402

FX = ROOT / "tests" / "fixtures"


def fixture_cpic() -> dict:
    pairs = json.loads((FX / "cpic_pairs.json").read_text())
    drugs = {d["drugid"]: d for d in json.loads((FX / "cpic_drugs.json").read_text())}
    guidelines = {g["id"]: g for g in json.loads((FX / "cpic_guidelines.json").read_text())}
    recs_by_guideline = json.loads((FX / "cpic_recommendations.json").read_text())
    recommendations = build.expand_recommendations([r for rows in recs_by_guideline.values() for r in rows])
    alleles = json.loads((FX / "cpic_alleles.json").read_text())
    return {"pairs": pairs, "drugs": drugs, "guidelines": guidelines, "recommendations": recommendations, "alleles": alleles}


def fixture_pharmgkb_zip() -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("clinical_annotations.tsv", (FX / "pharmgkb_clinical_annotations.tsv").read_text())
        zf.writestr("clinical_ann_alleles.tsv", (FX / "pharmgkb_clinical_ann_alleles.tsv").read_text())
        zf.writestr("clinical_ann_evidence.tsv", (FX / "pharmgkb_clinical_ann_evidence.tsv").read_text())
        zf.writestr("LICENSE.txt", (FX / "pharmgkb_LICENSE.txt").read_text())
        zf.writestr("CREATED_2026-01-15.txt", "")
    return buf.getvalue()


def fixture_hlaverify(alleles: list[str]) -> dict:
    rows = {a: {"current_name": a, "two_field": a, "flags": ""} for a in alleles}
    return {"url": build.HLAVERIFY_NORMALIZE, "ok": True, "ipd_release": "3.65.0", "alleles_checked": len(rows),
            "unresolvable": [], "rows": rows, "error": ""}


class ParseTests(unittest.TestCase):
    def test_lookup_key(self):
        # the real CPIC shape for ALLELE_STATUS genes
        self.assertEqual(build.parse_lookup_key({"HLA-B": "*57:01 positive"}), ("HLA-B", "HLA-B*57:01", "positive"))
        self.assertEqual(build.parse_lookup_keys({"CYP2C9": "*1/*1", "HLA-B": "*15:02 negative"}), [("HLA-B", "HLA-B*15:02", "negative")])
        self.assertEqual(build.parse_lookup_keys({"HLA-A": "*31:01 negative", "HLA-B": "*15:02 positive"}),
                         [("HLA-A", "HLA-A*31:01", "negative"), ("HLA-B", "HLA-B*15:02", "positive")])
        # tolerated alternative shape
        self.assertEqual(build.parse_lookup_key({"HLA-B*57:01": "positive"}), ("HLA-B", "HLA-B*57:01", "positive"))
        self.assertEqual(build.parse_lookup_key('{"HLA-A*31:01": "negative"}'), ("HLA-A", "HLA-A*31:01", "negative"))
        self.assertEqual(build.parse_lookup_key({"CYP2C9": "Normal Metabolizer"}), (None, None, None))
        gene, allele, status = build.parse_lookup_key({"CYP2C9": "Normal Metabolizer", "HLA-B*15:02": "positive"})
        self.assertEqual((gene, allele, status), ("HLA-B", "HLA-B*15:02", "positive"))
        self.assertEqual(build.parse_lookup_keys({"HLA-A*31:01": "negative", "HLA-B*15:02": "positive"}),
                         [("HLA-A", "HLA-A*31:01", "negative"), ("HLA-B", "HLA-B*15:02", "positive")])
        self.assertEqual(len(build.expand_recommendations([{"lookupkey": {"HLA-A*31:01": "negative", "HLA-B*15:02": "positive"}},
                                                           {"lookupkey": {"CYP2C9": "Normal Metabolizer"}}])), 2)

    def test_two_field(self):
        self.assertEqual(build.two_field("57:01"), "57:01")
        self.assertEqual(build.two_field("15:02:01"), "15:02")
        self.assertEqual(build.two_field("24:09N"), "24:09N")
        self.assertEqual(build.two_field("02:01:01:02L"), "02:01")

    def test_pharmgkb_allele_status(self):
        self.assertEqual(build.pharmgkb_allele_status("*57:01", "HLA-B"), ("HLA-B*57:01", "positive"))
        self.assertEqual(build.pharmgkb_allele_status("HLA-B*57:01 positive", "HLA-B"), ("HLA-B*57:01", "positive"))
        self.assertEqual(build.pharmgkb_allele_status("*57:01 negative", "HLA-B"), ("HLA-B*57:01", "negative"))
        self.assertEqual(build.pharmgkb_allele_status("no *15:02", "HLA-B"), ("HLA-B*15:02", "negative"))
        self.assertEqual(build.pharmgkb_allele_status("other", "HLA-B"), ("", ""))
        self.assertEqual(build.pharmgkb_allele_status("", "HLA-B"), ("", ""))
        self.assertEqual(build.pharmgkb_allele_status("*07:01", "HLA-DQA1", ["HLA-DRB1*07:01", "HLA-DQA1*02:01"]), ("HLA-DRB1*07:01", "positive"))
        self.assertEqual(build.pharmgkb_allele_status("*02:01", "HLA-DQA1", ["HLA-DRB1*07:01", "HLA-DQA1*02:01"]), ("HLA-DQA1*02:01", "positive"))

    def test_misc(self):
        self.assertTrue(build.is_hla_gene_field("HLA-DQA1, HLA-DRB1"))
        self.assertFalse(build.is_hla_gene_field("CYP2C9"))
        self.assertEqual(build.text_of({"HLA-B": "x"}), "x")
        self.assertEqual(build.text_of({"CYP2C9": "a", "HLA-B": "b"}), "CYP2C9: a; HLA-B: b")
        self.assertEqual(build.pmids(["1", "2"]), "1;2")
        self.assertEqual(build.pmids(None), "")
        self.assertTrue(build.ALLELE_RE.search("HLA-DRB1*07:01, HLA-DQA1*02:01"))


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name) / "data"

    def tearDown(self):
        self.tmp.cleanup()

    def run_build(self):
        with mock.patch.object(build, "fetch_cpic", fixture_cpic), \
             mock.patch.object(build, "fetch_pharmgkb", lambda: (fixture_pharmgkb_zip(), "fixture://pharmgkb")), \
             mock.patch.object(build, "check_alleles", fixture_hlaverify):
            rc = build.main(["--out", str(self.out)])
        self.assertEqual(rc, 0)

    def test_end_to_end_and_validate(self):
        self.run_build()
        for name in ("hla_pgx.csv", "hla_pgx.json", "cpic_hla_pairs.csv", "cpic_hla_recommendations.csv",
                     "cpic_hla_allele_frequencies.csv", "pharmgkb_hla_clinical_annotations.csv", "manifest.json",
                     "sources/pharmgkb_LICENSE.txt", "sources/cpic_pull.json"):
            self.assertTrue((self.out / name).exists(), name)
        rc = validate.main([str(self.out)])
        self.assertEqual(rc, 0, "validator must pass on the fixture build")

    def test_flattened_content(self):
        self.run_build()
        flat = json.loads((self.out / "hla_pgx.json").read_text())
        # the CYP2C9-only phenytoin row is not an HLA row and must be skipped
        self.assertFalse(any("no HLA key" in r["recommendation"] for r in flat))
        # phenytoin's HLA row keeps the full lookup key as its label
        pht = [r for r in flat if r["source"] == "cpic" and r["drug"] == "phenytoin"]
        self.assertTrue(pht and pht[0]["source_allele_label"].startswith("{"), pht)
        self.assertEqual(pht[0]["allele"], "HLA-B*15:02")
        # Level C pairs appear as gene–drug rows with a blank allele
        dap = [r for r in flat if r["source"] == "cpic" and r["drug"] == "dapsone"]
        self.assertEqual(len(dap), 1)
        self.assertEqual((dap[0]["allele"], dap[0]["evidence_level"], dap[0]["clinpgx_level"]), ("", "CPIC C", "2A"))
        # PharmGKB: CYP2C9 annotation filtered out; multi-drug annotation fans out; 'other' row is negative
        self.assertFalse(any(r["source"] == "pharmgkb" and r["gene"] == "CYP2C9" for r in flat))
        aba = [r for r in flat if r["source"] == "pharmgkb" and r["source_record_id"] == "1183621186"]
        self.assertEqual(sorted({r["drug"] for r in aba}), ["abacavir", "lamivudine"])
        other = [r for r in aba if r["source_allele_label"] == "other"]
        self.assertTrue(other and all(r["allele_status"] == "negative" and r["allele"] == "HLA-B*57:01" for r in other))
        # multi-gene annotation keeps the first HLA gene and parses the first allele
        lap = [r for r in flat if r["source"] == "pharmgkb" and r["drug"] == "lapatinib"]
        self.assertTrue(lap, "lapatinib rows present")
        pos = [r for r in lap if r["source_allele_label"] == "*07:01"]
        self.assertEqual((pos[0]["gene"], pos[0]["allele"], pos[0]["allele_status"]), ("HLA-DRB1", "HLA-DRB1*07:01", "positive"))
        oth = [r for r in lap if r["source_allele_label"] == "other"]
        self.assertEqual(oth[0]["allele_status"], "")  # two distinct alleles in the annotation: 'other' is not a clean negative
        # HLA-Verify columns filled from the (fixture) check
        b5701 = [r for r in flat if r["allele"] == "HLA-B*57:01"]
        self.assertTrue(b5701 and all(r["ipd_release"] == "3.65.0" for r in b5701))
        # every row has exactly the documented columns
        self.assertTrue(all(list(r.keys()) == build.FLAT_COLUMNS for r in flat))

    def test_manifest(self):
        self.run_build()
        m = json.loads((self.out / "manifest.json").read_text())
        self.assertEqual(m["license"], "CC-BY-SA-4.0")
        self.assertEqual(m["sources"]["cpic"]["pairs"], 16)
        self.assertEqual(m["sources"]["pharmgkb"]["annotations"], 15)
        self.assertEqual(m["sources"]["pharmgkb"]["created_stamp"], "CREATED_2026-01-15.txt")
        self.assertIn("hla_pgx.csv", m["files"])
        self.assertEqual(m["files"]["cpic_hla_pairs.csv"]["rows"], 16)

    def test_cpic_only_build(self):
        with mock.patch.object(build, "fetch_cpic", fixture_cpic), mock.patch.object(build, "check_alleles", fixture_hlaverify):
            rc = build.main(["--out", str(self.out), "--skip-pharmgkb"])
        self.assertEqual(rc, 0)
        self.assertFalse((self.out / "pharmgkb_hla_clinical_annotations.csv").exists())
        self.assertEqual(validate.main([str(self.out)]), 1)  # PharmGKB layer required by default
        self.assertEqual(validate.main([str(self.out), "--allow-no-pharmgkb"]), 0)


if __name__ == "__main__":
    unittest.main()
