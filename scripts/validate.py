#!/usr/bin/env python3
"""Check the invariants a usable hla-pgx-open snapshot must satisfy.

Usage:  python -I scripts/validate.py data [--allow-no-pharmgkb]

Exit 0 when every check passes, 1 otherwise. Checks are deliberately about facts that have been
stable for a decade (abacavir/HLA-B*57:01 is Level A and "Testing Required"); if one of them
ever fails the upstream changed shape or content and a human should look before publishing.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

TWO_FIELD_RE = re.compile(r"^HLA-[A-Z]+[0-9]*\*\d{2,4}:\d{2,4}[A-Z]?$")
CORE_ALLELES = {"HLA-B*57:01", "HLA-B*15:02", "HLA-B*58:01", "HLA-A*31:01"}


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("data", type=Path)
    ap.add_argument("--allow-no-pharmgkb", action="store_true")
    args = ap.parse_args(argv)
    d: Path = args.data
    problems: list[str] = []
    notes: list[str] = []

    def check(cond: bool, msg: str) -> None:
        (notes if cond else problems).append(("ok   " if cond else "FAIL ") + msg)

    manifest = json.loads((d / "manifest.json").read_text())
    pairs = read_csv(d / "cpic_hla_pairs.csv")
    recs = read_csv(d / "cpic_hla_recommendations.csv")
    freqs = read_csv(d / "cpic_hla_allele_frequencies.csv")
    flat = read_csv(d / "hla_pgx.csv")
    flat_json = json.loads((d / "hla_pgx.json").read_text())

    # --- CPIC pairs
    check(len(pairs) >= 10, f"CPIC has at least 10 HLA pairs (got {len(pairs)})")
    aba = [p for p in pairs if p["gene"] == "HLA-B" and p["drug"].lower() == "abacavir"]
    check(bool(aba), "abacavir/HLA-B pair present")
    if aba:
        check(aba[0]["cpic_level"] == "A", f"abacavir/HLA-B is CPIC Level A (got {aba[0]['cpic_level']!r})")
        check(aba[0]["pgx_testing"] == "Testing Required", f"abacavir/HLA-B is 'Testing Required' (got {aba[0]['pgx_testing']!r})")
        check(bool(aba[0]["guideline_url"]), "abacavir/HLA-B carries a guideline URL")
        check(bool(aba[0]["citations_pmid"]), "abacavir/HLA-B carries PMIDs")
    for gene, drug in (("HLA-B", "carbamazepine"), ("HLA-A", "carbamazepine"), ("HLA-B", "allopurinol"), ("HLA-B", "phenytoin")):
        check(any(p["gene"] == gene and p["drug"].lower() == drug and p["cpic_level"] == "A" for p in pairs),
              f"{drug}/{gene} present at CPIC Level A")

    # --- CPIC recommendations
    check(len(recs) >= 8, f"at least 8 HLA allele-status recommendation rows (got {len(recs)})")
    aba_pos = [r for r in recs if r["allele"] == "HLA-B*57:01" and r["allele_status"] == "positive" and r["drug"].lower() == "abacavir"]
    check(bool(aba_pos) and bool(aba_pos[0]["recommendation"]), "HLA-B*57:01 positive → abacavir recommendation text present")
    check(all(TWO_FIELD_RE.match(r["allele"]) for r in recs), "every recommendation allele is a two-field HLA name")
    level_a_pairs = {(p["gene"], p["drug_id"]) for p in pairs if p["cpic_level"] == "A" and p["used_for_recommendation"] in ("True", "true", "1")}
    rec_pairs = {(r["gene"], r["drug_id"]) for r in recs}
    missing = sorted(level_a_pairs - rec_pairs)
    check(not missing, f"every CPIC Level A pair used for recommendation has at least one recommendation row (missing: {missing})")

    # --- CPIC frequencies
    freq_alleles = {f["allele"] for f in freqs}
    check(CORE_ALLELES <= freq_alleles, f"frequencies cover the four core alleles (have {sorted(freq_alleles)})")
    groups = {f["population_group"] for f in freqs}
    check(len(groups) >= 5, f"at least 5 biogeographic groups in the frequency table (got {len(groups)})")
    check(all(0.0 <= float(f["allele_frequency"]) <= 1.0 for f in freqs if f["allele_frequency"] != ""), "frequencies are within [0,1]")

    # --- PharmGKB
    pg_path = d / "pharmgkb_hla_clinical_annotations.csv"
    if pg_path.exists():
        pg = read_csv(pg_path)
        ids = {r["clinical_annotation_id"] for r in pg}
        check(len(ids) >= 15, f"at least 15 PharmGKB HLA clinical annotations (got {len(ids)})")
        check(any("57:01" in r["variant_haplotypes"] and "abacavir" in r["drugs"].lower() and r["level_of_evidence"].startswith("1")
                  for r in pg), "PharmGKB has a Level 1 HLA-B*57:01/abacavir annotation")
        check(all(r["url"].startswith("http") for r in pg), "every PharmGKB row has a URL")
        check(sum(1 for r in pg if r["annotation_text"]) >= 0.8 * len(pg), "most PharmGKB rows carry allele annotation text")
        check(bool(manifest["sources"]["pharmgkb"].get("license_file_in_zip")), "PharmGKB license file found and copied")
    else:
        check(args.allow_no_pharmgkb, "PharmGKB layer present (pass --allow-no-pharmgkb to build CPIC-only)")

    # --- flattened table
    check(len(flat) == len(flat_json), f"hla_pgx.csv and hla_pgx.json have the same row count ({len(flat)} vs {len(flat_json)})")
    check(len(flat) >= len(recs), "flattened table has at least one row per CPIC recommendation")
    bad = [r["allele"] for r in flat if r["allele"] and not TWO_FIELD_RE.match(r["allele"])]
    check(not bad, f"every non-empty allele in hla_pgx.csv is a two-field HLA name (bad: {bad[:5]})")
    check(all(r["source"] in ("cpic", "pharmgkb") for r in flat), "source column is cpic|pharmgkb")
    check(all(r["allele_status"] in ("", "positive", "negative", "carrier") for r in flat), "allele_status vocabulary")
    check(any(r["source"] == "cpic" and r["allele"] == "HLA-B*57:01" and r["allele_status"] == "positive" and r["drug"] == "abacavir" for r in flat),
          "flattened table has the CPIC HLA-B*57:01 positive / abacavir row")

    # --- HLA-Verify check (soft)
    hv = manifest["sources"]["hlaverify"]
    if hv.get("ok"):
        checked = hv.get("alleles_checked", 0)
        unres = hv.get("unresolvable", [])
        check(checked > 0, f"HLA-Verify checked {checked} alleles against IPD-IMGT/HLA {hv.get('ipd_release')}")
        check(len(unres) <= max(1, 0.2 * checked), f"≤20% of alleles unresolvable in IPD-IMGT/HLA (unresolvable: {unres})")
        if unres:
            notes.append(f"note HLA-Verify could not resolve: {unres} — inspect; may be a deprecated or renamed name upstream")
    else:
        notes.append(f"note HLA-Verify check not run ({hv.get('error', 'skipped')}); hlaverify_* columns are empty")

    # --- manifest hashes match files
    for rel, entry in manifest["files"].items():
        p = d / rel
        if not p.exists():
            problems.append(f"FAIL manifest lists missing file {rel}")
            continue
        import hashlib
        if hashlib.sha256(p.read_bytes()).hexdigest() != entry["sha256"]:
            problems.append(f"FAIL sha256 mismatch for {rel}")
    notes.append(f"ok   manifest covers {len(manifest['files'])} files")

    for line in notes + problems:
        print(line)
    print(f"\n{len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
