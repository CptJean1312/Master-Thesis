# Protection Data Provider Clarification

Date noted: 2026-07-09

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

Tan also noted that he will send the list of municipalities and AGS in the provider exposure data so that it can be compared with the current RP500 corridor extent.

## Thesis interpretation

The absence of a municipality from `elbe_protection_level_mun.csv` should not be interpreted as observed lack of flood protection. It should be interpreted as absence from the provider's modeled damage/loss portfolio.

For mapping and descriptive summaries, the current defensible wording is:

- `Protection value available`: municipality is included in the damage/loss portfolio and has an estimated protection return period.
- `No modeled loss in portfolio`: municipality is absent from the protection table because it does not experience modeled loss in the provider flood portfolio.

## Follow-up needed

Compare the provider exposure-data municipality list with the current RP500 corridor once Tan sends it. This will clarify whether the difference between the provider portfolio and the EFAS/JRC RP500 corridor is mainly due to different hazard/exposure definitions, modeled asset/loss representation, or catchment/corridor extent.
