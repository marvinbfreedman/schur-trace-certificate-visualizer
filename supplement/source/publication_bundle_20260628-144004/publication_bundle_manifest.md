# Publication Bundle Manifest

Target route:

1. arXiv preprint in `math.NT` as the primary outlet.
2. Companion Zenodo DOI for code, theorem JSONs, and regeneration artifacts.

## arXiv Package

- Main manuscript: `rh_weyl_positive_draft.tex`
- Bibliography is already embedded in the manuscript source.
- The paper should cite the exact theorem chain and point readers to the
  reproducibility bundle below.

## Companion Reproducibility Bundle

Include the following in the Zenodo archive or a linked repository:

- `publication_audit_dependency_graph.py`
- `publication_audit_dependency_graph.json`
- `publication_audit_dependency_graph.md`
- `publication_audit_rh_formal_conclusion.json`
- `publication_audit_rh_formal_conclusion.md`
- `publication_audit_rh_ledger_full.json`
- `publication_audit_rh_ledger_full.md`
- `publication_audit_external_equivalence_full.json`
- `publication_audit_external_equivalence_full.md`
- `rh_klm_notes.md`
- `rh_weyl_positive_draft.tex`
- all theorem/certificate `.py` and `.json` files used by the proof spine

## Recommended Publication Order

1. Freeze the manuscript PDF from `rh_weyl_positive_draft.tex`.
2. Tag the exact code/data state used for that PDF.
3. Create a Zenodo release for the tagged bundle.
4. Submit the manuscript to arXiv.
5. After arXiv is public, circulate the repo/Zenodo link for external review.

## Notes

- The audit ledgers are useful as reproducibility material, but they do not
  replace the manuscript.
- Keep the manuscript self-contained: the proof must be in the paper, with the
  ledger and code as supporting artifacts.
- If the paper is later revised, version the Zenodo archive alongside the
  manuscript so the correspondence stays exact.
