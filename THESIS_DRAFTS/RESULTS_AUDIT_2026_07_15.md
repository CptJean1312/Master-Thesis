# Results Audit - Chapter 5

Date: 2026-07-15

Active draft:

- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/THESIS_DRAFTS/THESIS DRAFT.md`

Purpose:

This memo records the check of Chapter 5 after writing the full Results prose. The aim was to verify whether the numerical claims, denominators, figure logic, and provider/protection interpretation are supported by the current analysis outputs.

## Bottom Line

Chapter 5 is numerically consistent with the current output tables and maps. The main results are reproducible from the exported analysis tables.

Three corrections were made during the audit:

1. The JRC flood-hazard dataset citation in Chapter 4 was corrected from `Baugh et al. (2026)` to `Baugh et al. (2024)`, following the official JRC catalogue entry.
2. Two Results formulations were softened from a broad "more vulnerable municipalities" phrasing to a more precise "upper vulnerability groups" / "suggests" wording, because the protection/loss pattern is clearer than the exposure pattern but not perfectly monotonic.
3. Chapter 4.10 was aligned with the actual Results workflow by removing the implied regression/model-check step. The final workflow is descriptive and distributive: maps, correlations, vulnerability quintiles, bivariate classes, exposure curves, and protection/loss summaries.

## Verified Core Numbers

### Sample and coverage

Checked against:

- `ANALYSE_FINAL_CORRIDOR/outputs/tables/corridor_analysis_rp100.csv`
- `ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/tables/corridor_protection_level_all_corridor_municipalities.csv`
- `ANALYSE_FINAL_CORRIDOR/PROTECTION/PROTECTION_DATA_PROVIDER_CLARIFICATION.md`

Verified:

- RP500 corridor municipalities: `835`
- Municipalities with final vulnerability index: `834`
- Missing vulnerability record: `16076094` Berga-Wuenschendorf
- Original provider protection/loss table: `301` municipalities
- Provider municipalities matching the RP500 corridor: `280`
- Corridor municipalities with no simulated loss event / zero modeled annual loss probability: `555`

### Exposure

Checked against:

- `corridor_analysis_rp100.csv`

Verified RP100 values:

- mean flood share: `25.4%`
- median flood share: `12.9%`
- Q25: `3.9%`
- Q75: `33.4%`
- zero RP100 exposure municipalities: `14`
- municipalities with RP100 share >= `99.9%`: `8`

Verified multi-return-period summary:

- RP10 mean / median: `20.3%` / `8.9%`
- RP500 mean / median: `27.0%` / `13.9%`
- zero exposure declines from `67` at RP10 to `0` at RP500

### Vulnerability

Checked against:

- `corridor_analysis_rp100.csv`
- `corridor_vulnerability_dimension_correlations.csv`

Verified:

- vulnerability index range: `-2.75` to `4.02`
- median: `-0.04`
- access vs demographic dimension correlation: `-0.215`
- access vs deprivation dimension correlation: `-0.053`
- demographic vs deprivation dimension correlation: `0.051`

### Exposure-vulnerability relationship

Checked against:

- `table_results_key_summary.csv`
- recalculation from `corridor_analysis_rp100.csv`

Verified:

- Pearson correlation vulnerability vs RP100 exposure: `-0.050`
- Spearman correlation vulnerability vs RP100 exposure: `-0.050`
- RP100 mean by vulnerability quintile:
  - Q1 lowest: `35.6%`
  - Q2: `21.9%`
  - Q3: `23.0%`
  - Q4: `19.2%`
  - Q5 highest: `27.2%`
- RP100 median by vulnerability quintile:
  - Q1 lowest: `17.3%`
  - Q5 highest: `16.3%`
- high RP100 exposure tercile counts by vulnerability tercile:
  - low vulnerability: `104`
  - medium vulnerability: `81`
  - high vulnerability: `93`

Interpretation check:

The text is correct to say that the exposure-vulnerability relationship is weak, slightly negative, and non-monotonic. It is also correct to say that high-high municipalities exist but do not dominate the high-exposure group.

### Exposure curves

Checked against:

- `ANALYSE_FINAL_CORRIDOR/EXPOSURE_CURVES/outputs/tables/corridor_exposure_curve_type_summary.csv`
- `ANALYSE_FINAL_CORRIDOR/EXPOSURE_CURVES/outputs/tables/corridor_exposure_curve_metric_summary_by_type.csv`
- `ANALYSE_FINAL_CORRIDOR/outputs/results_figures/table_mean_exposure_curves_by_vulnerability_quintile.csv`

Verified:

- early exposure: `772` municipalities / `92.46%`
- delayed jump: `46` municipalities / `5.51%`
- gradual increase: `17` municipalities / `2.04%`
- early-exposure municipalities realize on average `93.4%` of final RP500 exposure by RP100
- delayed-jump municipalities realize on average `15.6%` by RP100 and add `84.4%` after RP100
- Q1 exposure curve: `31.4%` at RP10 to `36.8%` at RP500
- Q5 exposure curve: `20.8%` at RP10 to `29.4%` at RP500
- Q4 exposure curve: `14.5%` at RP10, `19.2%` at RP100, `20.9%` at RP500

Interpretation check:

The text is correct to say that RP100 is not hiding a fundamentally different vulnerability pattern at rarer return periods. The vulnerability gradient remains non-linear.

### Protection and loss portfolio

Checked against:

- `table_portfolio_status_by_vulnerability_quintile.csv`
- `corridor_protection_level_all_corridor_municipalities.csv`
- `PROTECTION_STATUS_MEMO_CURRENT.md`
- `NO_LOSS_STREAM_CONTEXT_CHECK.md`

Verified positive modeled loss share by vulnerability quintile:

- Q1 lowest: `16.2%`
- Q2: `30.5%`
- Q3: `33.5%`
- Q4: `46.7%`
- Q5 highest: `41.0%`

Verified annual loss probability by vulnerability quintile:

- Q1 mean: `0.9%`
- Q2 mean: `1.3%`
- Q3 mean: `2.0%`
- Q4 mean: `4.8%`
- Q5 mean: `3.3%`
- median is `0` in all quintiles
- Q5 upper quartile: `1.88%`

Verified finite protection return period among positive-loss municipalities:

- Q1 median: `156` years
- Q2 median: `61` years
- Q3 median: `68` years
- Q4 median: `40` years
- Q5 median: `28` years

Verified correlations:

- vulnerability vs modeled annual loss probability, full corridor with no-event as zero:
  - Pearson: `0.118`
  - Spearman: `0.230`
- vulnerability vs finite return period, positive-loss subset:
  - Pearson: `-0.067`
  - Spearman: `-0.198`

Interpretation check:

The text is correct to say that the protection/loss layer shows a clearer social pattern than RP100 area exposure. The text is also correct to keep the claim cautious: the relationship is weak to moderate, not a clean monotonic staircase.

## Belegbarkeit and Citation Check

### Solid / supported

- All Chapter 5 numerical claims are supported by current CSV/GPKG-derived outputs.
- Figure paths in Chapter 5 point to existing files.
- The no-event interpretation is supported by local provider-clarification memos.
- The limitation that smaller streams are not represented is documented in the provider-clarification memos and should be connected to Sairam et al. (2021) for the broader simulation context.

### Corrected during audit

- JRC dataset citation:
  - corrected to `Baugh et al. (2024)`
  - DOI: `10.2905/1D128B6C-A4EE-4858-9E34-6210707F3C81`
- `DATA_SOURCES_CURRENT.md` was updated accordingly.

### Still needs final bibliography polish

Before submission, the final bibliography should explicitly include:

- Baugh et al. (2024) JRC dataset entry
- Dottori et al. (2022)
- Alfieri et al. (2014)
- Sairam et al. (2021)
- BBSR/INKAR data source and access/version information
- VG250/BKG municipal boundary source

The provider clarification should be cited in text as personal communication or documented in an appendix. Chapter 4 now mentions it as:

`N. Sairam and T. N. M. Phan, personal communication, July 2026`

## Remaining Non-Critical Layout Issue

Figure 5.1 duplicates the study-area map already used as Figure 2.1. This is not a data error. It is a final layout decision:

- keep Figure 2.1 only and refer back to it in Chapter 5, or
- keep Figure 5.1 as a Results-chapter reminder and remove/rename the Chapter 2 version.

## Audit Verdict

No substantive numerical error was found in Chapter 5. The empirical story is defensible:

1. RP100 area exposure is highly uneven but only weakly and slightly negatively associated with socio-economic vulnerability.
2. Exposure curves confirm that most corridor municipalities realize much of their RP500 exposure by RP100 or earlier.
3. The provider loss/protection layer shows a clearer social pattern, with positive modeled losses and lower finite return periods more common in upper vulnerability groups.
4. The no-event class is meaningful but must be interpreted cautiously because it combines high modeled protection cases and model-scope limitations, especially smaller streams.
