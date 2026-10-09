# Validation and review

Implementation baseline: `d166be532720abed66f1d9816b148207ccaa6064`.
Review scope: new CHD Quarto document, repo-local offline modules/tests and
its own rendered outputs. Existing website changes and other posts are excluded.

## Behavioral and integration checks

- All five test files pass via `Rscript offline/run-tests.R`, including fresh
  local input reconstruction, blocked scientific remote calls and full exports.
- All 14 R source/test/runner files parse successfully. R has no static
  typechecking configuration in this repository.
- The compact real GBC regression fixture matches pinned-reference estimates.
  An additional development check compared the full original local GBC bundle
  against the exact reference commit: all 33 pair scientific fields match.
- CHD retains 41 leads, 18 target SNPs, 23 pair rows, five pooled estimates and
  22 usable comparisons/decompositions. All 23 total associations remain usable.
  SBP/rs10139550 remains unavailable for indirect/residual inference.
- Outcome extraction coverage is 40/41 leads, 984/992 distinct instrument SNPs
  and 18/18 target SNPs. Expanded VCF duplicate records are retained for checks.
- CSV round trips preserve row counts. Geographic/covariance sensitivity CSVs
  have 46/69 rows, respectively, and 12 separate inference families are exported.
- A fresh report-only Quarto render succeeds with scientific remote access
  blocked. All 10 PNGs are embedded; native MathML needs no external scripts.
  All 35 local download links resolve.
- All exported pages were visually inspected: comparison PDF 5 pages,
  decomposition PDF 4 pages and pooled PDF 1 page. Dense LDL/SBP labels and
  unavailable components remain readable.

## Standards

Independent standards review found no actionable documented-standard violations
or material maintainability findings. Input reconstruction, statistical core,
inference, plotting and orchestration follow the documented offline workflow.

## Spec

Independent specification review identified one provenance defect: queried
counts were derived from successful outcome matches. It was corrected to use
independent inventory denominators, with regression assertions for 40/41 and
984/992 coverage. The report now distinguishes returned association records
from requested instrument identities. No other missing requirements, scientific
implementation errors or scope creep were identified. The user's new-QMD
instruction supersedes the original issue's notebook-edit delivery requirement.

Final findings: Standards 0; Spec 1 corrected, 0 unresolved.
