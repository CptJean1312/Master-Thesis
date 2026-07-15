# Protection Data Provider Clarification

Date noted: 2026-07-09

Updated: 2026-07-15 after final clarification from Nivedita.

This note records the provider clarification received after the first join between `elbe_protection_level_mun.csv` and the RP500 Elbe corridor municipalities.

## Current join result

- RP500 corridor municipalities: `835`
- Municipalities in `elbe_protection_level_mun.csv`: `301`
- Protection-table municipalities matching the RP500 corridor: `280`
- RP500 corridor municipalities absent from the protection table: `555`
- Protection-table municipalities outside the RP500 corridor: `21`

## Final provider clarification

The protection table contains municipalities that have damage/loss in the provider flood portfolio. According to Tan's clarification, municipalities that are not present in the table do not experience modeled losses in that model portfolio.

Nivedita and Tan later clarified that no municipalities were removed because of too few events. There is no 30-event threshold in the current protection table. All municipalities in the protection data are derived directly from the loss data without any additional event-count filter.

Nivedita later clarified the substantive interpretation of municipalities without loss events. For larger cities such as Hamburg and Magdeburg, she would interpret the absence of events as a protection level high enough that no simulated loss events are observed. She recommended treating the probability as zero. She also clarified that smaller streams are not considered in the simulation and that this should be written as a limitation. The simulation characteristics and limitations should be referenced through Sairam et al. (2021), `Process-Based Flood Risk Assessment for Germany`, DOI: `10.1029/2021EF002259`.

Tan also noted that he will send the list of municipalities and AGS in the provider exposure data so that it can be compared with the current RP500 corridor extent, if needed.

## Thesis interpretation

The absence of a municipality from `elbe_protection_level_mun.csv` should not be interpreted as observed lack of flood protection. It should be interpreted as no simulated loss event in the provider model output. Following the final provider clarification, modeled annual loss probability for these municipalities can be treated as `0`.

For finite protection return-period analysis, however, these no-event municipalities should not be assigned an arbitrary finite return period. They are better treated as no-event / zero-probability cases, while finite protection return periods are interpreted only for municipalities with positive modeled losses.

For mapping and descriptive summaries, the current defensible wording is:

- `Protection value available`: municipality has positive modeled loss in the provider output and an estimated finite protection return period.
- `No simulated loss event`: municipality has modeled annual loss probability `0` in the provider output; for larger cities this can reflect high protection levels, while smaller-stream effects are limited by model coverage.

## Follow-up needed

If Tan sends the provider exposure-data municipality list, compare it with the current RP500 corridor. This can clarify whether all no-event municipalities are inside the provider exposure domain but protected/no-loss, or whether some are absent because the represented river network excludes smaller streams.
