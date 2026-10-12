#!/usr/bin/env python3
"""Build the hla-pgx-open dataset from CPIC and PharmGKB, with IPD-IMGT/HLA allele checks.

Standard library only. Run:  python -I scripts/build.py --out data

Sources
-------
CPIC      https://api.cpicpgx.org/v1/  (PostgREST; tables pair, recommendation, guideline, drug, allele)
PharmGKB  clinicalAnnotations.zip from the PharmGKB / ClinPGx download API
HLA-Verify  POST https://api.hlaverify.com/v1/normalize  (anonymous tier; non-fatal if unreachable)

Everything written under --out is CC BY-SA 4.0 (see LICENSE-DATA); this code is Apache-2.0.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

BUILDER_VERSION = "0.1.0"
USER_AGENT = f"hla-pgx-open/{BUILDER_VERSION} (+https://github.com/jasonbrelsford/hla-pgx-open)"

CPIC_BASE = "https://api.cpicpgx.org/v1"
HLA_GENES_CPIC = ("HLA-A", "HLA-B")
PHARMGKB_URLS = (
    "https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip",
    "https://api.clinpgx.org/v1/download/file/data/clinicalAnnotations.zip",
)
HLAVERIFY_NORMALIZE = "https://api.hlaverify.com/v1/normalize"
HLAVERIFY_BATCH = 200  # under the anonymous tier's per-call cap

# HLA-B*57:01, HLA-DRB1*07:01, HLA-B*15:02 — gene prefix, two fields, optional expression suffix.
ALLELE_RE = re.compile(r"\bHLA-([A-Z]+[0-9]*)\*(\d{2,4}:\d{2,4}(?::\d{2,4}){0,2}[A-Z]?)\b")
TWO_FIELD_RE = re.compile(r"^HLA-[A-Z]+[0-9]*\*\d{2,4}:\d{2,4}[A-Z]?$")

FLAT_COLUMNS = [
    "source", "source_record_id", "gene", "allele", "allele_status", "source_allele_label",
    "drug", "drug_id", "evidence_level", "clinpgx_level", "pgx_testing", "phenotype_category",
    "implication", "recommendation", "classification", "population", "score",
    "guideline_name", "guideline_url", "annotation_url", "citations_pmid",
    "hlaverify_current_name", "hlaverify_2field", "hlaverify_flags", "ipd_release", "source_version",
]


def log(msg: str) -> None:
    print(f"[build] {msg}", file=sys.stderr, flush=True)


# ----------------------------------------------------------------------------- HTTP helpers
def http_get(url: str, *, accept: str = "application/json", retries: int = 3, timeout: int = 60) -> bytes:
    last: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": accept})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:  # noqa: PERF203
            last = e
            log(f"GET {url} failed (attempt {attempt}/{retries}): {e}")
            time.sleep(2 * attempt)
    raise RuntimeError(f"GET {url} failed after {retries} attempts: {last}")


def http_post_json(url: str, payload: dict, *, timeout: int = 60) -> dict:
    body = json.dumps(payload).encode()
    req = urllib.request.Request(url, data=body, method="POST",
                                 headers={"User-Agent": USER_AGENT, "Accept": "application/json",
                                          "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def cpic_get(table: str, **params: str) -> list[dict]:
    qs = urllib.parse.urlencode(params, safe="().,:*")  # PostgREST operators stay readable
    url = f"{CPIC_BASE}/{table}" + (f"?{qs}" if qs else "")
    data = json.loads(http_get(url).decode())
    if not isinstance(data, list):
        raise RuntimeError(f"CPIC {table}: expected a JSON array, got {type(data).__name__}: {str(data)[:200]}")
    return data


# ----------------------------------------------------------------------------- CPIC
FIELDS_RE = re.compile(r"\*(\d{2,4}:\d{2,4}(?::\d{2,4}){0,2}[A-Z]?)")


def parse_lookup_keys(lookupkey) -> list[tuple[str, str, str | None]]:
    """CPIC lookup keys for ALLELE_STATUS genes come as {"HLA-B": "*57:01 positive"}; the allele-keyed
    form {"HLA-B*57:01": "positive"} is accepted too. Returns [(gene, allele, status)] for every HLA
    entry; non-HLA keys (CYP2C9 in the phenytoin guideline) are skipped."""
    if isinstance(lookupkey, str):
        try:
            lookupkey = json.loads(lookupkey)
        except json.JSONDecodeError:
            return []
    if not isinstance(lookupkey, dict):
        return []
    out = []
    for key, value in lookupkey.items():
        key_s, val_s = str(key).strip(), ("" if value is None else str(value).strip())
        m = ALLELE_RE.search(key_s)
        if m:  # allele in the key, status in the value
            gene = f"HLA-{m.group(1)}"
            out.append((gene, f"{gene}*{two_field(m.group(2))}", val_s.lower() or None))
            continue
        if not key_s.upper().startswith("HLA-"):
            continue
        gene = key_s.upper()
        mv = ALLELE_RE.search(val_s) or FIELDS_RE.search(val_s)
        if not mv:
            continue
        fields = mv.group(2) if mv.re is ALLELE_RE else mv.group(1)
        if mv.re is ALLELE_RE:
            gene = f"HLA-{mv.group(1)}"
        low = val_s.lower()
        status = "negative" if "negative" in low else ("positive" if "positive" in low else None)
        out.append((gene, f"{gene}*{two_field(fields)}", status))
    return out


def parse_lookup_key(lookupkey) -> tuple[str | None, str | None, str | None]:
    """First HLA allele in a lookup key, or (None, None, None)."""
    found = parse_lookup_keys(lookupkey)
    return found[0] if found else (None, None, None)


def expand_recommendations(rows: list[dict]) -> list[dict]:
    """One entry per HLA allele named in each recommendation's lookup key; rows with no HLA key are dropped
    (a CYP2C9-only row in the phenytoin guideline, for example)."""
    out = []
    for r in rows:
        for gene, allele, status in parse_lookup_keys(r.get("lookupkey")):
            e = dict(r)
            e["_gene"], e["_allele"], e["_status"] = gene, allele, status
            out.append(e)
    return out


def two_field(fields: str) -> str:
    parts = fields.split(":")
    suffix = ""
    if parts[-1][-1:].isalpha():
        suffix = parts[-1][-1]
        parts[-1] = parts[-1][:-1]
    return ":".join(parts[:2]) + (suffix if len(parts) <= 2 else "")


def pmids(citations) -> str:
    if not citations:
        return ""
    if isinstance(citations, str):
        try:
            citations = json.loads(citations)
        except json.JSONDecodeError:
            citations = [c.strip() for c in citations.split(",")]
    return ";".join(str(c) for c in citations if str(c).strip())


def fetch_cpic() -> dict:
    log("CPIC: pairs")
    pairs: list[dict] = []
    for gene in HLA_GENES_CPIC:
        pairs.extend(cpic_get("pair", genesymbol=f"eq.{gene}"))
    pairs = [p for p in pairs if not p.get("removed")]
    if not pairs:
        raise RuntimeError("CPIC returned no HLA pairs")

    log("CPIC: drugs")
    drugs = {d["drugid"]: d for d in cpic_get("drug", select="drugid,name")}

    guideline_ids = sorted({p["guidelineid"] for p in pairs if p.get("guidelineid")})
    log(f"CPIC: guidelines {guideline_ids}")
    guidelines: dict[int, dict] = {}
    for gid in guideline_ids:
        rows = cpic_get("guideline", id=f"eq.{gid}")
        if rows:
            guidelines[gid] = rows[0]

    log("CPIC: recommendations")
    recommendations: list[dict] = []
    for gid in guideline_ids:
        recommendations.extend(expand_recommendations(cpic_get("recommendation", guidelineid=f"eq.{gid}")))

    log("CPIC: alleles (frequencies)")
    alleles: list[dict] = []
    for gene in HLA_GENES_CPIC:
        alleles.extend(cpic_get("allele", genesymbol=f"eq.{gene}"))

    return {"pairs": pairs, "drugs": drugs, "guidelines": guidelines,
            "recommendations": recommendations, "alleles": alleles}


def cpic_pair_rows(c: dict) -> list[dict]:
    rows = []
    for p in c["pairs"]:
        g = c["guidelines"].get(p.get("guidelineid")) or {}
        d = c["drugs"].get(p["drugid"]) or {}
        rows.append({
            "pairid": p["pairid"], "gene": p["genesymbol"], "drug_id": p["drugid"], "drug": d.get("name", ""),
            "cpic_level": p.get("cpiclevel") or "", "clinpgx_level": p.get("clinpgxlevel") or "",
            "pgx_testing": p.get("pgxtesting") or "", "used_for_recommendation": p.get("usedforrecommendation"),
            "guideline_id": p.get("guidelineid") or "", "guideline_name": g.get("name", ""),
            "guideline_url": g.get("url", ""), "guideline_clinpgx_id": g.get("clinpgxid") or g.get("pharmgkbid") or "",
            "citations_pmid": pmids(p.get("citations")), "pair_version": p.get("version", ""),
            "removed": p.get("removed", False), "removed_date": p.get("removeddate") or "",
            "removed_reason": p.get("removedreason") or "",
        })
    rows.sort(key=lambda r: (r["gene"], str(r["cpic_level"]), r["drug"], r["pairid"]))
    return rows


def cpic_recommendation_rows(c: dict) -> list[dict]:
    rows = []
    for r in c["recommendations"]:
        g = c["guidelines"].get(r.get("guidelineid")) or {}
        d = c["drugs"].get(r.get("drugid")) or {}
        rows.append({
            "recommendation_id": r.get("id"), "guideline_id": r.get("guidelineid"), "guideline_name": g.get("name", ""),
            "drug_id": r.get("drugid"), "drug": d.get("name", ""), "gene": r["_gene"], "allele": r["_allele"],
            "allele_status": r["_status"] or "", "implication": text_of(r.get("implications")),
            "recommendation": r.get("drugrecommendation") or "", "classification": r.get("classification") or "",
            "population": r.get("population") or "", "comments": r.get("comments") or "",
            "dosing_information": r.get("dosinginformation"), "alternate_drug_available": r.get("alternatedrugavailable"),
            "other_prescribing_guidance": r.get("otherprescribingguidance"),
            "lookup_key_json": json.dumps(r.get("lookupkey"), sort_keys=True) if r.get("lookupkey") is not None else "",
            "recommendation_version": r.get("version", ""),
        })
    rows.sort(key=lambda r: (r["gene"], r["allele"], r["drug"], r["allele_status"], r["population"], str(r["recommendation_id"])))
    return rows


def text_of(v) -> str:
    """CPIC 'implications' is a {gene: text} object; flatten to 'GENE: text; GENE: text' or the plain string."""
    if v is None:
        return ""
    if isinstance(v, str):
        try:
            v = json.loads(v)
        except json.JSONDecodeError:
            return v
    if isinstance(v, dict):
        if len(v) == 1:
            return str(next(iter(v.values())))
        return "; ".join(f"{k}: {val}" for k, val in sorted(v.items()))
    return str(v)


def cpic_frequency_rows(c: dict) -> list[dict]:
    rows = []
    for a in c["alleles"]:
        name = a.get("name") or ""
        gene = a.get("genesymbol") or ""
        allele = name if name.startswith("HLA-") else (f"{gene}{name}" if name.startswith("*") else f"{gene}*{name}")
        freq = a.get("frequency") or {}
        if isinstance(freq, str):
            try:
                freq = json.loads(freq)
            except json.JSONDecodeError:
                freq = {}
        for group, value in sorted((freq or {}).items()):
            rows.append({"gene": gene, "allele": allele, "population_group": group, "allele_frequency": value,
                         "allele_record_version": a.get("version", ""), "inferred_frequency": a.get("inferredfrequency")})
    rows.sort(key=lambda r: (r["gene"], r["allele"], r["population_group"]))
    return rows


# ----------------------------------------------------------------------------- PharmGKB
def fetch_pharmgkb() -> tuple[bytes, str]:
    last = None
    for url in PHARMGKB_URLS:
        try:
            log(f"PharmGKB: downloading {url}")
            data = http_get(url, accept="application/zip, application/octet-stream, */*", retries=2, timeout=180)
            zipfile.ZipFile(io.BytesIO(data)).testzip()
            return data, url
        except Exception as e:  # noqa: BLE001
            last = e
            log(f"PharmGKB: {url} unusable: {e}")
    raise RuntimeError(f"PharmGKB clinicalAnnotations.zip unavailable from all known URLs: {last}")


def read_tsv(zf: zipfile.ZipFile, suffix: str) -> list[dict]:
    names = [n for n in zf.namelist() if n.lower().endswith(suffix.lower())]
    if not names:
        raise RuntimeError(f"PharmGKB zip has no file ending in {suffix}; members: {zf.namelist()[:20]}")
    text = zf.read(names[0]).decode("utf-8-sig")
    return list(csv.DictReader(io.StringIO(text), delimiter="\t", quoting=csv.QUOTE_NONE))


def is_hla_gene_field(gene_field: str) -> bool:
    return any(part.strip().upper().startswith("HLA-") for part in (gene_field or "").split(","))


def pharmgkb_allele_status(label: str, gene: str, candidates: list[str] | None = None) -> tuple[str, str]:
    """Return (allele, status) for a PharmGKB Genotype/Allele label like '*57:01', 'HLA-B*57:01 positive',
    '*57:01 negative', 'other', 'no *57:01'. A bare '*07:01' is attributed to the gene among `candidates`
    (the annotation's Variant/Haplotypes alleles) that carries those fields, else to `gene`.
    Unmappable labels return ('', '')."""
    raw = (label or "").strip()
    low = raw.lower()
    m = ALLELE_RE.search(raw)
    if m:
        g, fields = f"HLA-{m.group(1)}", m.group(2)
    else:
        m2 = re.search(r"\*(\d{2,4}:\d{2,4}(?::\d{2,4}){0,2}[A-Z]?)", raw)
        if not m2:
            return "", ""
        fields = m2.group(1)
        g = gene
        tf = two_field(fields)
        owners = sorted({c.split("*")[0] for c in (candidates or []) if c.split("*", 1)[1:] == [tf]})
        if len(owners) == 1:
            g = owners[0]
    allele = f"{g}*{two_field(fields)}"
    negative = any(w in low for w in ("negative", "non-carrier", "noncarrier", "absence", "absent", "without", " not ", "no *", "homozygous for other"))
    if low.startswith("not ") or low.startswith("no "):
        negative = True
    return allele, ("negative" if negative else "positive")


def pharmgkb_rows(zip_bytes: bytes) -> tuple[list[dict], dict]:
    zf = zipfile.ZipFile(io.BytesIO(zip_bytes))
    created = next((n for n in zf.namelist() if Path(n).name.upper().startswith("CREATED")), "")
    license_name = next((n for n in zf.namelist() if "license" in Path(n).name.lower()), "")
    license_text = zf.read(license_name).decode("utf-8", "replace") if license_name else ""

    ann = read_tsv(zf, "clinical_annotations.tsv")
    alleles = read_tsv(zf, "clinical_ann_alleles.tsv")
    evidence = read_tsv(zf, "clinical_ann_evidence.tsv")
    headers = {"clinical_annotations.tsv": list(ann[0].keys()) if ann else [],
               "clinical_ann_alleles.tsv": list(alleles[0].keys()) if alleles else [],
               "clinical_ann_evidence.tsv": list(evidence[0].keys()) if evidence else []}

    ev_pmids: dict[str, set[str]] = defaultdict(set)
    for e in evidence:
        pm = (e.get("PMID") or "").strip()
        if pm:
            ev_pmids[e.get("Clinical Annotation ID", "").strip()].add(pm)

    allele_rows: dict[str, list[dict]] = defaultdict(list)
    for a in alleles:
        allele_rows[a.get("Clinical Annotation ID", "").strip()].append(a)

    out = []
    kept = 0
    for r in ann:
        gene_field = r.get("Gene", "")
        variants = r.get("Variant/Haplotypes", "")
        if not (is_hla_gene_field(gene_field) or "HLA-" in variants.upper()):
            continue
        kept += 1
        cid = r.get("Clinical Annotation ID", "").strip()
        base = {
            "clinical_annotation_id": cid, "gene": gene_field, "variant_haplotypes": variants,
            "level_of_evidence": r.get("Level of Evidence", ""), "level_override": r.get("Level Override", ""),
            "level_modifiers": r.get("Level Modifiers", ""), "score": r.get("Score", ""),
            "phenotype_category": r.get("Phenotype Category", ""), "pmid_count": r.get("PMID Count", ""),
            "evidence_count": r.get("Evidence Count", ""), "drugs": r.get("Drug(s)", ""),
            "phenotypes": r.get("Phenotype(s)", ""), "latest_history_date": r.get("Latest History Date (YYYY-MM-DD)", r.get("Latest History Date", "")),
            "url": r.get("URL", ""), "specialty_population": r.get("Specialty Population", ""),
            "evidence_pmids": ";".join(sorted(ev_pmids.get(cid, ()), key=lambda s: (len(s), s))),
        }
        rows = allele_rows.get(cid) or [{}]
        for a in rows:
            row = dict(base)
            row["genotype_allele"] = a.get("Genotype/Allele", "")
            row["annotation_text"] = a.get("Annotation Text", "")
            row["allele_function"] = a.get("Allele Function", "")
            out.append(row)
    out.sort(key=lambda r: (r["gene"], r["variant_haplotypes"], r["drugs"], r["clinical_annotation_id"], r["genotype_allele"]))
    meta = {"created_stamp": Path(created).name if created else "", "license_file": license_name,
            "license_text": license_text, "annotations": kept, "allele_rows": len(out), "members": zf.namelist(),
            "headers": headers}
    return out, meta


# ----------------------------------------------------------------------------- HLA-Verify
def check_alleles(alleles: list[str], *, max_calls: int = 60) -> dict:
    """Run every distinct allele through HLA-Verify's /v1/normalize. Non-fatal: a failing batch is
    bisected so one bad name (or a server error on one input) does not blank every column."""
    result = {"url": HLAVERIFY_NORMALIZE, "ok": False, "ipd_release": "", "alleles_checked": 0,
              "unresolvable": [], "rows": {}, "error": "", "failed_inputs": [], "errors": [], "calls": 0}
    names = sorted(set(a for a in alleles if a))
    if not names:
        return result

    def call(batch: list[str]) -> None:
        if result["calls"] >= max_calls:
            result["failed_inputs"].extend(batch)
            return
        result["calls"] += 1
        try:
            resp = http_post_json(HLAVERIFY_NORMALIZE, {"typings": batch})
        except urllib.error.HTTPError as e:
            body = ""
            try:
                body = e.read().decode("utf-8", "replace")[:300]
            except Exception:  # noqa: BLE001
                pass
            result["errors"].append({"status": e.code, "batch_size": len(batch), "body": body})
            if len(batch) == 1:
                result["failed_inputs"].extend(batch)
                return
            mid = len(batch) // 2
            call(batch[:mid])
            call(batch[mid:])
            return
        except Exception as e:  # noqa: BLE001
            result["errors"].append({"status": None, "batch_size": len(batch), "body": str(e)[:300]})
            result["failed_inputs"].extend(batch)
            return
        result["ipd_release"] = resp.get("release", result["ipd_release"])
        for row in resp.get("rows", []):
            reported = row.get("reported")
            result["rows"][reported] = {"current_name": row.get("current_name", ""),
                                       "two_field": row.get("allele_2field", ""),
                                       "flags": ";".join(row.get("flags") or [])}
            if row.get("allele_2field") == "UNRESOLVABLE":
                result["unresolvable"].append(reported)

    for i in range(0, len(names), HLAVERIFY_BATCH):
        call(names[i:i + HLAVERIFY_BATCH])
    result["alleles_checked"] = len(result["rows"])
    result["ok"] = result["alleles_checked"] > 0
    if not result["ok"]:
        result["error"] = (result["errors"][0]["body"] or f"HTTP {result['errors'][0]['status']}") if result["errors"] else "no rows returned"
    log(f"HLA-Verify: checked {result['alleles_checked']}/{len(names)} alleles against IPD-IMGT/HLA {result['ipd_release'] or '?'} "
        f"in {result['calls']} call(s); unresolvable: {result['unresolvable'] or 'none'}; failed inputs: {result['failed_inputs'] or 'none'}")
    return result


# ----------------------------------------------------------------------------- flatten
def flatten(c: dict, pair_rows: list[dict], rec_rows: list[dict], pg_rows: list[dict], hv: dict, pg_meta: dict) -> list[dict]:
    pairs_by_key = {(p["gene"], p["drug_id"]): p for p in pair_rows}
    flat: list[dict] = []

    def hv_cols(allele: str) -> dict:
        r = hv["rows"].get(allele) or {}
        return {"hlaverify_current_name": r.get("current_name", ""), "hlaverify_2field": r.get("two_field", ""),
                "hlaverify_flags": r.get("flags", ""), "ipd_release": hv["ipd_release"] if r else ""}

    # CPIC rows with recommendations: one per recommendation row.
    covered_pairs = set()
    for r in rec_rows:
        p = pairs_by_key.get((r["gene"], r["drug_id"])) or {}
        covered_pairs.add((r["gene"], r["drug_id"]))
        lk = json.loads(r["lookup_key_json"]) if r["lookup_key_json"] else {}
        label = f"{r['allele']} {r['allele_status']}".strip() if len(lk) <= 1 else r["lookup_key_json"]
        flat.append({
            "source": "cpic", "source_record_id": f"{p.get('pairid', '')}/{r['recommendation_id']}",
            "gene": r["gene"], "allele": r["allele"], "allele_status": r["allele_status"], "source_allele_label": label,
            "drug": r["drug"], "drug_id": r["drug_id"], "evidence_level": f"CPIC {p.get('cpic_level', '')}".strip(),
            "clinpgx_level": p.get("clinpgx_level", ""), "pgx_testing": p.get("pgx_testing", ""), "phenotype_category": "",
            "implication": r["implication"], "recommendation": r["recommendation"], "classification": r["classification"],
            "population": r["population"], "score": "", "guideline_name": r["guideline_name"],
            "guideline_url": p.get("guideline_url", ""), "annotation_url": "", "citations_pmid": p.get("citations_pmid", ""),
            **hv_cols(r["allele"]),
            "source_version": f"pair:{p.get('pair_version', '')};recommendation:{r['recommendation_version']}",
        })
    # CPIC pairs with no allele-level recommendation (Level B/C associations): gene–drug row, allele blank.
    for p in pair_rows:
        if (p["gene"], p["drug_id"]) in covered_pairs:
            continue
        flat.append({
            "source": "cpic", "source_record_id": f"{p['pairid']}/", "gene": p["gene"], "allele": "", "allele_status": "",
            "source_allele_label": "", "drug": p["drug"], "drug_id": p["drug_id"],
            "evidence_level": f"CPIC {p['cpic_level']}".strip(), "clinpgx_level": p["clinpgx_level"],
            "pgx_testing": p["pgx_testing"], "phenotype_category": "", "implication": "", "recommendation": "",
            "classification": "", "population": "", "score": "", "guideline_name": p["guideline_name"],
            "guideline_url": p["guideline_url"], "annotation_url": "", "citations_pmid": p["citations_pmid"],
            "hlaverify_current_name": "", "hlaverify_2field": "", "hlaverify_flags": "", "ipd_release": "",
            "source_version": f"pair:{p['pair_version']}",
        })
    # PharmGKB rows: one per (annotation allele row, drug).
    for r in pg_rows:
        gene_first = next((g.strip() for g in r["gene"].split(",") if g.strip().upper().startswith("HLA-")), r["gene"].strip())
        found = [f"HLA-{m.group(1)}*{two_field(m.group(2))}" for m in ALLELE_RE.finditer(r["variant_haplotypes"])]
        allele, status = pharmgkb_allele_status(r["genotype_allele"], gene_first, found)
        if not allele:
            allele = found[0] if found else ""
            if len(set(found)) == 1 and r["genotype_allele"].strip().lower().startswith("other"):
                status = "negative"
        row_gene = allele.split("*")[0] if allele else gene_first
        drugs = [d.strip() for d in r["drugs"].split(",") if d.strip()] or [""]
        for drug in drugs:
            flat.append({
                "source": "pharmgkb", "source_record_id": r["clinical_annotation_id"], "gene": row_gene,
                "allele": allele, "allele_status": status, "source_allele_label": r["genotype_allele"],
                "drug": drug, "drug_id": "", "evidence_level": f"PharmGKB {r['level_of_evidence']}".strip(),
                "clinpgx_level": r["level_of_evidence"], "pgx_testing": "", "phenotype_category": r["phenotype_category"],
                "implication": r["annotation_text"], "recommendation": "", "classification": "",
                "population": r["specialty_population"], "score": r["score"], "guideline_name": "",
                "guideline_url": "", "annotation_url": r["url"], "citations_pmid": r["evidence_pmids"],
                **hv_cols(allele),
                "source_version": pg_meta.get("created_stamp", ""),
            })
    flat.sort(key=lambda r: (r["gene"], r["allele"], r["source"], r["drug"], r["allele_status"], r["source_record_id"]))
    return flat


# ----------------------------------------------------------------------------- output
def write_csv(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    columns = columns or (list(rows[0].keys()) if rows else [])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in columns})


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:  # noqa: BLE001
        return None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=Path("data"))
    ap.add_argument("--skip-hlaverify", action="store_true", help="do not call api.hlaverify.com")
    ap.add_argument("--skip-pharmgkb", action="store_true", help="build the CPIC layer only (manifest records it)")
    args = ap.parse_args(argv)
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    (out / "sources").mkdir(exist_ok=True)

    built_at = dt.datetime.now(dt.timezone.utc).replace(microsecond=0)
    snapshot_date = built_at.date().isoformat()

    cpic = fetch_cpic()
    pair_rows = cpic_pair_rows(cpic)
    rec_rows = cpic_recommendation_rows(cpic)
    freq_rows = cpic_frequency_rows(cpic)
    log(f"CPIC: {len(pair_rows)} pairs, {len(rec_rows)} HLA recommendations, {len(freq_rows)} frequency rows")

    pg_rows: list[dict] = []
    pg_meta: dict = {"ok": False}
    pg_url = ""
    if not args.skip_pharmgkb:
        zip_bytes, pg_url = fetch_pharmgkb()
        pg_rows, pg_meta = pharmgkb_rows(zip_bytes)
        pg_meta["ok"] = True
        if pg_meta.get("license_text"):
            (out / "sources" / "pharmgkb_LICENSE.txt").write_text(pg_meta["license_text"], encoding="utf-8")
        log(f"PharmGKB: {pg_meta['annotations']} HLA annotations, {pg_meta['allele_rows']} allele rows ({pg_meta['created_stamp']})")

    alleles = sorted({r["allele"] for r in rec_rows} | {r["allele"] for r in freq_rows}
                     | {pharmgkb_allele_status(r["genotype_allele"], r["gene"].split(",")[0].strip(),
                                               [f"HLA-{m.group(1)}*{two_field(m.group(2))}" for m in ALLELE_RE.finditer(r["variant_haplotypes"])])[0] for r in pg_rows}
                     | {f"HLA-{m.group(1)}*{two_field(m.group(2))}" for r in pg_rows for m in ALLELE_RE.finditer(r["variant_haplotypes"])})
    alleles = [a for a in alleles if a]
    hv = {"url": HLAVERIFY_NORMALIZE, "ok": False, "ipd_release": "", "alleles_checked": 0, "unresolvable": [], "rows": {}, "error": "skipped"}
    if not args.skip_hlaverify:
        hv = check_alleles(alleles)

    flat = flatten(cpic, pair_rows, rec_rows, pg_rows, hv, pg_meta)

    write_csv(out / "cpic_hla_pairs.csv", pair_rows)
    write_csv(out / "cpic_hla_recommendations.csv", rec_rows)
    write_csv(out / "cpic_hla_allele_frequencies.csv", freq_rows)
    if pg_rows:
        write_csv(out / "pharmgkb_hla_clinical_annotations.csv", pg_rows)
    write_csv(out / "hla_pgx.csv", flat, FLAT_COLUMNS)
    (out / "hla_pgx.json").write_text(json.dumps(flat, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")

    (out / "sources" / "cpic_pull.json").write_text(json.dumps({
        "base_url": CPIC_BASE, "pulled_at": built_at.isoformat(), "genes": HLA_GENES_CPIC,
        "guideline_ids": sorted(cpic["guidelines"]), "pair_ids": sorted(p["pairid"] for p in cpic["pairs"]),
        "allele_names": sorted(a.get("name", "") for a in cpic["alleles"]),
    }, indent=1) + "\n", encoding="utf-8")

    files = {}
    for p in sorted(out.rglob("*")):
        if p.is_file() and p.name != "manifest.json":
            rel = p.relative_to(out).as_posix()
            entry = {"sha256": sha256_file(p), "bytes": p.stat().st_size}
            if p.suffix == ".csv":
                with p.open(encoding="utf-8") as f:
                    entry["rows"] = sum(1 for _ in f) - 1
            files[rel] = entry

    manifest = {
        "name": "hla-pgx-open", "snapshot_date": snapshot_date, "built_at": built_at.isoformat(), "git_commit": git_commit(),
        "license": "CC-BY-SA-4.0", "builder": f"scripts/build.py@{BUILDER_VERSION}",
        "sources": {
            "cpic": {"base_url": CPIC_BASE, "ok": True, "pairs": len(pair_rows), "recommendations": len(rec_rows),
                     "allele_frequency_rows": len(freq_rows), "guideline_ids": sorted(cpic["guidelines"])},
            "pharmgkb": {"url_used": pg_url, "ok": bool(pg_meta.get("ok")), "created_stamp": pg_meta.get("created_stamp", ""),
                         "annotations": pg_meta.get("annotations", 0), "allele_rows": pg_meta.get("allele_rows", 0),
                         "license_file_in_zip": pg_meta.get("license_file", ""), "headers": pg_meta.get("headers", {})},
            "hlaverify": {k: v for k, v in hv.items() if k != "rows"},
        },
        "files": files,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    log(f"wrote {len(files)} files + manifest.json to {out}/ ; flattened rows: {len(flat)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
