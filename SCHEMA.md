# Schema

All files are UTF-8, comma-separated, header row first, one record per line, RFC 4180 quoting.
Empty string means "not provided by the source", never zero or false. Allele names are written
with the gene prefix and two fields (`HLA-B*57:01`); `hlaverify_*` columns hold what the
IPD-IMGT/HLA check returned for that name.

## `hla_pgx.csv` — the flattened table

One row per (source, source record, allele, allele status, drug). A PharmGKB annotation that
names three alleles and two drugs becomes six rows; a CPIC recommendation that applies to
"HLA-B\*57:01 positive" for abacavir is one row.

| column | type | meaning |
|---|---|---|
| `source` | `cpic` \| `pharmgkb` | Which upstream the row comes from |
| `source_record_id` | string | CPIC `pairid` + recommendation `id` (`110078/8479702`), or the PharmGKB clinical annotation id |
| `gene` | string | `HLA-A`, `HLA-B`, … |
| `allele` | string | Two-field allele with gene prefix, as written by the source |
| `allele_status` | `positive` \| `negative` \| `carrier` \| `` | CPIC: the allele status the recommendation applies to. PharmGKB: the genotype/allele label of the annotation row (`*57:01`, `*57:01 negative`, …) mapped to positive/negative when unambiguous, else copied verbatim into `source_allele_label` and left blank here |
| `source_allele_label` | string | The allele/genotype label exactly as the source wrote it |
| `drug` | string | Drug name as the source writes it (lower case for CPIC) |
| `drug_id` | string | CPIC: `RxNorm:…` id. PharmGKB: empty (the download carries names only) |
| `evidence_level` | string | CPIC: `CPIC A`, `CPIC B/C`, `CPIC C`. PharmGKB: `PharmGKB 1A` … `PharmGKB 4` |
| `clinpgx_level` | string | CPIC pair's `clinpgxlevel` (`1A`, `2A`, `3`, …), the PharmGKB level CPIC records for the pair |
| `pgx_testing` | string | CPIC pair's testing recommendation (`Testing Required`, `Testing Recommended`, `Actionable PGx`, `Informative PGx`) |
| `phenotype_category` | string | PharmGKB `Phenotype Category` (`Toxicity`, `Efficacy`, `Dosage`, `Metabolism/PK`, `Other`); empty for CPIC |
| `implication` | string | CPIC `implications` text for the allele status, or the PharmGKB per-allele annotation text |
| `recommendation` | string | CPIC `drugrecommendation` text; empty for PharmGKB (annotations describe associations, they do not prescribe) |
| `classification` | string | CPIC recommendation strength (`Strong`, `Moderate`, `Optional`, `No Recommendation`) |
| `population` | string | CPIC recommendation `population` (`general`, …) or PharmGKB `Specialty Population` |
| `score` | number | PharmGKB clinical annotation score; empty for CPIC |
| `guideline_name` | string | CPIC guideline title |
| `guideline_url` | url | CPIC guideline page (currently on clinpgx.org) |
| `annotation_url` | url | PharmGKB clinical annotation page |
| `citations_pmid` | string | Semicolon-separated PMIDs (CPIC: the pair's guideline citations; PharmGKB: PMIDs of the annotation's evidence rows) |
| `hlaverify_current_name` | string | The current IPD-IMGT/HLA name the allele resolves to, or `UNRESOLVABLE` |
| `hlaverify_2field` | string | The comparable two-field name |
| `hlaverify_flags` | string | Semicolon-separated flags from the check (`deprecated_name`, `null_allele`, …); empty when none |
| `ipd_release` | string | IPD-IMGT/HLA release the check ran against (e.g. `3.65.0`); empty if the check was skipped |
| `source_version` | string | CPIC: the `version` integers of the pair and recommendation rows; PharmGKB: the `CREATED_*` stamp of the download |
| `snapshot_date` | date | UTC date of the build |

## `cpic_hla_pairs.csv`

Direct projection of the CPIC `pair` table for `HLA-A` and `HLA-B`, one row per pair, joined
to `drug` and `guideline`:
`pairid, gene, drug_id, drug, cpic_level, clinpgx_level, pgx_testing, used_for_recommendation, guideline_id, guideline_name, guideline_url, guideline_clinpgx_id, citations_pmid, pair_version, removed, removed_date, removed_reason`.

## `cpic_hla_recommendations.csv`

One row per CPIC `recommendation` whose `lookupkey` names an HLA allele:
`recommendation_id, guideline_id, guideline_name, drug_id, drug, gene, allele, allele_status, implication, recommendation, classification, population, comments, dosing_information, alternate_drug_available, other_prescribing_guidance, recommendation_version`.

`allele` and `allele_status` are parsed from `lookupkey` (`{"HLA-B*57:01": "positive"}`). A
recommendation whose lookup key involves another gene as well (phenytoin also keys on CYP2C9)
keeps the full key in `lookup_key_json` so nothing is lost in the projection.

## `cpic_hla_allele_frequencies.csv`

One row per (allele, biogeographic group) from the CPIC `allele` table's frequency object:
`gene, allele, population_group, allele_frequency, allele_record_version, inferred_frequency`.
These are **allele** frequencies, not carrier frequencies; carrier frequency for a
Hardy–Weinberg population is `1 − (1 − p)²`.

## `pharmgkb_hla_clinical_annotations.csv`

One row per (clinical annotation, allele row) for annotations whose `Gene` column is an HLA
gene, carrying the annotation-level fields on every row:
`clinical_annotation_id, gene, variant_haplotypes, level_of_evidence, level_override, level_modifiers, score, phenotype_category, pmid_count, evidence_count, drugs, phenotypes, latest_history_date, url, specialty_population, genotype_allele, annotation_text, allele_function, evidence_pmids`.

## `manifest.json`

```
{
  "snapshot_date": "YYYY-MM-DD",
  "built_at": "ISO-8601 UTC",
  "git_commit": "sha or null",
  "sources": {
    "cpic": {"base_url": "...", "ok": true, "pairs": 16, "recommendations": N, "alleles": N},
    "pharmgkb": {"url_used": "...", "ok": true, "created_stamp": "...", "annotations": N, "allele_rows": N},
    "hlaverify": {"url": "...", "ok": true, "ipd_release": "3.65.0", "alleles_checked": N, "unresolvable": [...]}
  },
  "files": {"hla_pgx.csv": {"rows": N, "sha256": "..."}, ...},
  "license": "CC-BY-SA-4.0",
  "builder": "scripts/build.py@<version>"
}
```
