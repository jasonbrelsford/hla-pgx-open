# hla-pgx-open

**An open, machine-readable table of HLA–drug associations (pharmacogenomics), built from CPIC and PharmGKB/ClinPGx, with every allele name checked against the current IPD-IMGT/HLA release.**

Status: **v0 — first snapshots.** Data is rebuilt by CI every Monday and on demand; each build is a dated snapshot with a manifest. Treat it as a well-labelled mirror of the upstream sources, not as a new source of truth.

Companion dataset to [HLA-Verify](https://hlaverify.com) (deterministic HLA nomenclature verification). The API endpoints that consume this table live there; the data lives here so that it stays open, citable and separately licensed.

## What is in `data/`

| File | Rows are | Source |
|---|---|---|
| `hla_pgx.csv` / `hla_pgx.json` | **The flattened table**: one row per (source, allele, allele status, drug) | CPIC + PharmGKB |
| `cpic_hla_pairs.csv` | CPIC gene–drug pairs for HLA-A and HLA-B with CPIC level, ClinPGx level, testing recommendation, guideline, PMIDs | CPIC API `pair` |
| `cpic_hla_recommendations.csv` | CPIC allele-status recommendations (implication, recommendation text, strength, population) | CPIC API `recommendation` |
| `cpic_hla_allele_frequencies.csv` | CPIC allele frequencies by biogeographic group for the actionable HLA alleles | CPIC API `allele` |
| `pharmgkb_hla_clinical_annotations.csv` | PharmGKB/ClinPGx clinical annotations whose gene is an HLA gene, with level of evidence, score, phenotype category and per-allele annotation text | PharmGKB `clinicalAnnotations.zip` |
| `manifest.json` | Snapshot date, which upstream URLs answered, row counts, SHA-256 of every file, the IPD-IMGT/HLA release used for allele checks | build |
| `sources/` | Upstream license files and raw-pull notes copied verbatim at build time | upstream |

Column definitions: [`SCHEMA.md`](SCHEMA.md).

## What it is for

- Giving an agent, a LIMS, or a notebook one file to join a **two-field HLA allele** to its **known drug associations**, with the evidence level and the citation on the same row.
- Pairing with an HLA nomenclature normalizer. Lab reports and registry records arrive as G groups, MAC codes, legacy spellings and four-field names; PharmCAT, the reference PGx tool, [requires an outside two-field HLA call](https://pharmcat.clinpgx.org/using/Calling-HLA/) and will not impute one. Normalize first, then look up here.
- Teaching: the whole actionable HLA layer of pharmacogenomics is small enough to read in one sitting (CPIC lists 16 HLA gene–drug pairs, 7 of them Level A with guideline recommendations, at the 2026-10 snapshot).

## What it is not

**Not clinical decision support and not a medical device.** This is a lookup table over published guidelines. Any clinical use goes back to the cited CPIC guideline or PharmGKB annotation, and to a clinician. Nothing here is validated for patient care, and nothing here should receive patient data.

In particular, after an allogeneic haematopoietic cell transplant the recipient's blood carries the **donor's** genome, and which genome should drive an HLA-mediated drug alert is an open research question (see the NCI natural-history study [NCT06856226](https://clinicaltrials.gov/study/NCT06856226) and [Langman et al., J Mol Diagn 2016](https://mdanderson.elsevierpure.com/en/publications/challenges-in-determining-genotypes-for-pharmacogenetics-in-allog)). This table says nothing about that question; it only makes the allele-level facts easy to look up for either party.

## Using it

Raw files, always the latest snapshot on `main`:

```
https://raw.githubusercontent.com/jasonbrelsford/hla-pgx-open/main/data/hla_pgx.csv
https://raw.githubusercontent.com/jasonbrelsford/hla-pgx-open/main/data/hla_pgx.json
https://raw.githubusercontent.com/jasonbrelsford/hla-pgx-open/main/data/manifest.json
```

Pin a snapshot by commit SHA (shown in `manifest.json`) or by a release tag once releases start.

```python
import csv, urllib.request
url = "https://raw.githubusercontent.com/jasonbrelsford/hla-pgx-open/main/data/hla_pgx.csv"
rows = list(csv.DictReader(urllib.request.urlopen(url).read().decode().splitlines()))
hits = [r for r in rows if r["allele"] == "HLA-B*57:01" and r["allele_status"] == "positive"]
for r in hits:
    print(r["source"], r["drug"], r["evidence_level"], r["recommendation"], r["guideline_url"] or r["annotation_url"])
```

## How it is built

`scripts/build.py` (Python standard library only) does, in order:

1. **CPIC** — reads `pair`, `recommendation`, `guideline`, `drug` and `allele` from the public CPIC PostgREST API (`https://api.cpicpgx.org/v1/`) for `HLA-A` and `HLA-B`.
2. **PharmGKB / ClinPGx** — downloads `clinicalAnnotations.zip`, keeps annotations whose gene is `HLA-*`, and joins the per-allele annotation text and the evidence rows. The license file shipped inside the zip is copied to `data/sources/`.
3. **Allele check** — sends every distinct allele name through HLA-Verify's public `POST /v1/normalize` (anonymous tier) and records the current name, the comparable two-field name, any flags, and the IPD-IMGT/HLA release they were checked against. If the API is unreachable the columns are left empty and the manifest says so; the build does not fail on it.
4. **Flatten + validate** — writes the files above, then `scripts/validate.py` checks the invariants that must hold for the snapshot to be usable (for example: abacavir/HLA-B\*57:01 is present, Level A, "Testing Required"; every CPIC Level A pair has at least one recommendation row; every allele parses as a two-field name).

Run it yourself:

```
python -I scripts/build.py --out data
python -I scripts/validate.py data
python -I -m unittest discover -s tests
```

The CI workflow (`.github/workflows/build.yml`) runs the same three commands weekly and commits the result when anything changed, with the snapshot date in the commit message.

## Licenses

- **Data** (`data/`): [Creative Commons Attribution-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/) — see [`LICENSE-DATA`](LICENSE-DATA) and [`ATTRIBUTION.md`](ATTRIBUTION.md). This is inherited: CPIC and PharmGKB publish their content under CC BY-SA 4.0, so a derived table must carry the same license and credit them. Reuse it freely, keep the attribution, share derivatives alike.
- **Code** (`scripts/`, `tests/`, workflow): [Apache License 2.0](LICENSE-CODE).

## Citing

If you use the table, cite the upstream sources (they did the work) and this repository as the packaging. `CITATION.cff` has the repository entry; `ATTRIBUTION.md` has the upstream citations, including the PMIDs of each CPIC guideline present in the snapshot.

## Roadmap

- [ ] First CI snapshot on `main` (this repo's first job).
- [ ] Confirm the exact license text on both upstream sites against `LICENSE-DATA` before tagging **v1.0.0** ([CPIC licensing](https://cpicpgx.org/?p=1245) · [PharmGKB data usage policy](https://www.pharmgkb.org/page/dataUsagePolicy)).
- [ ] Zenodo DOI per release (enable the GitHub integration at <https://zenodo.org/account/settings/github/>, then tag).
- [ ] HLA-Verify: `GET /v1/pgx/table`, `POST /v1/pgx/screen`, `POST /v1/pgx/pair` and an MCP `pgx_lookup` tool that read this table ([tracking issue](https://github.com/jasonbrelsford/verifiable-science-envs/issues/110)).
- [ ] Donor–recipient discordance rates by population — waits for a public donor–recipient pair dataset; not computed from non-public data.

## Maintainer

Jason Brelsford — [hello@hlaverify.com](mailto:hello@hlaverify.com). Issues and pull requests welcome; corrections to upstream facts belong upstream (CPIC and PharmGKB both take feedback), and this repo will pick them up on the next build.
