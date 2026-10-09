# Offline CHD follow-up

The current analysis is `../chd_followup.qmd`, replacing the historical
`../chd_example.ipynb` result sections. The notebook is preserved for acquisition
context and is not executed. Run from its directory:

```sh
Rscript offline/run-tests.R
Rscript offline/run-analysis.R
sh offline/render-report.sh
```

`run-analysis.R` optionally accepts input and output directories as positional
arguments. `reconstruct_chd_inputs(input_dir)` reads exactly
`chd_phewas_result.Rdata`, the five named mediator VCFs and `ukb-d-I9_IHD.vcf.gz`
with existing local tabix indices. Both command and report reconstruct fresh
associations. Results and compact derived input CSVs go in
`../chd-followup-results`; no result cache, network calls or installation is
required. The isolated Quarto wrapper avoids website-wide generated changes.

Dependencies: R, TwoSampleMR, ieugwasr, gwasvcf, GenomicRanges, Rsamtools, VariantAnnotation, knitr,
rmarkdown, testthat, and Quarto. Actual installed versions are exported as
`package_versions.csv`. VCF/Rdata source checksums and extraction counts are
exported in `provenance.csv` and `extraction_counts.csv`; source filenames are
relative to the selected input directory. Large VCFs and their indices are
ignored by the existing repository policy.

## Statistical source and adaptations

`R/offline-{core,pooled,decomposition,inference,plots}.R` are derived from
MRCIEU/gbc-mr, commit `3bbac89e254acfbfa3385aae9b4520e4e8b5932c`:
https://github.com/MRCIEU/gbc-mr/tree/3bbac89e254acfbfa3385aae9b4520e4e8b5932c.
They preserve its conservative allele policy, TwoSampleMR Wald/IVW calls,
heterogeneity formula, delta uncertainty and hypothetical covariance machinery.
Adaptations accept independent complete lead coordinates and outcome metadata,
reuse one computed comparator in pooled output, expose retained identities,
remove GBC-specific scale assumptions and paginate dense mediator figures.

`run_offline_analysis(inputs, output_dir, figures_dir)` is outcome-neutral.
`inputs` supplies `int_chd` (exact mediator associations in outcome-column
schema), `int_inst` (mediator instruments in exposure-column schema),
`outgwasf` (outcome associations), `lead_inventory`
(`target_snp/target_chr/target_pos`), and `outcome_metadata` (`name/id/units/build`).
Historical `phewas_result` does not define CHD's exclusion inventory. Missing
or conflicting coordinates fail closed. Build metadata are saved-source
claims, not independent validation or liftover.

The primary policy excludes all leads and inclusive ±1 Mb regions. Separate
±2 Mb and ±2 Mb/MHC comparison families and hypothetical equicorrelation
decomposition families are exported. Missing tests remain missing, all pair
identities are retained, and outcome units remain explicitly unknown. See
the report for interpretations, limits and every downloadable figure page.
