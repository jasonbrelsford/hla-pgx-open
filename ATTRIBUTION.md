# Attribution

The data in this repository is derived from two upstream sources. Both publish under
CC BY-SA 4.0; this file is the attribution that license requires, and the citations their
authors ask for. If you redistribute or build on the data, carry this attribution with it.

## CPIC — Clinical Pharmacogenetics Implementation Consortium

- Content: gene–drug pairs, allele-status recommendations, guideline metadata and allele
  frequencies, read from the public CPIC API (`https://api.cpicpgx.org/v1/`).
- Website: <https://cpicpgx.org> · Licensing: <https://cpicpgx.org/?p=1245>
- Consortium citation: Relling MV, Klein TE. CPIC: Clinical Pharmacogenetics Implementation
  Consortium of the Pharmacogenomics Research Network. *Clin Pharmacol Ther.* 2011;89(3):464-467.
  PMID 21270786, doi:[10.1038/clpt.2010.279](https://doi.org/10.1038/clpt.2010.279).
- Guideline citations: each row of `cpic_hla_pairs.csv` carries the PMIDs of the guideline
  publication(s) it comes from (`citations_pmid`). At the first snapshot these were the CPIC
  guidelines for abacavir and HLA-B; carbamazepine and oxcarbazepine and HLA-A/HLA-B;
  allopurinol and HLA-B; and phenytoin/fosphenytoin and CYP2C9/HLA-B. Cite the guideline, not
  only this table, when a specific recommendation is used.
- CPIC® is a registered service mark of the U.S. Department of Health & Human Services. Use of
  the data does not imply endorsement by CPIC.

## PharmGKB / ClinPGx

- Content: clinical annotations (level of evidence, score, phenotype category, per-allele
  annotation text, evidence list) for HLA genes, from the public download
  `clinicalAnnotations.zip`.
- Website: <https://www.pharmgkb.org> (now also <https://www.clinpgx.org>) ·
  Data usage policy: <https://www.pharmgkb.org/page/dataUsagePolicy>
- Citation: Whirl-Carrillo M, Huddart R, Gong L, et al. An evidence-based framework for
  evaluating pharmacogenomics knowledge for personalized medicine. *Clin Pharmacol Ther.*
  2021;110(3):563-572. PMID 34216021, doi:[10.1002/cpt.2350](https://doi.org/10.1002/cpt.2350).
- The license file that ships inside the download is copied verbatim to
  `data/sources/pharmgkb_LICENSE.txt` at build time.

## IPD-IMGT/HLA (via HLA-Verify)

Allele names are checked against the current IPD-IMGT/HLA release through HLA-Verify's public
`/v1/normalize` endpoint; the release used is recorded in `manifest.json` and on every row
(`ipd_release`). IPD-IMGT/HLA citation: Barker DJ, Maccari G, Georgiou X, et al. The
IPD-IMGT/HLA Database. *Nucleic Acids Res.* 2023;51(D1):D1053-D1060. PMID 36350643, doi:[10.1093/nar/gkac1011](https://doi.org/10.1093/nar/gkac1011).

## This repository

Brelsford J. *hla-pgx-open: an open HLA–drug association table built from CPIC and PharmGKB.*
GitHub, 2026. <https://github.com/jasonbrelsford/hla-pgx-open>. (A Zenodo DOI per release is
on the roadmap; use it once it exists.)

## How to attribute in one line

> HLA–drug associations from CPIC (cpicpgx.org) and PharmGKB/ClinPGx (pharmgkb.org), packaged
> by hla-pgx-open (github.com/jasonbrelsford/hla-pgx-open), snapshot YYYY-MM-DD, CC BY-SA 4.0.
