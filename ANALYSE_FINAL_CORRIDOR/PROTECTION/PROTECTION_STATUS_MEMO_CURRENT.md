# Protection Status Memo - Current Working Interpretation

Date: 2026-07-13

Purpose: freeze the current interpretation of the newly received municipality-level protection/loss data before writing Chapter 4 and the Results chapter.

## Data files and scripts

Provider data:

- `/Users/maxi_161/Desktop/UNI/Master/THESIS/DATEN + GIS/PROTECTION/elbe_protection_level_mun.csv`

Main thesis-side scripts:

- `ANALYSE_FINAL_CORRIDOR/PROTECTION/BUILD_PROTECTION_MAPS.R`
- `ANALYSE_FINAL_CORRIDOR/PROTECTION/CHECK_NO_LOSS_STREAM_CONTEXT.R`

Main outputs:

- `ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/maps/map_corridor_protection_coverage.png`
- `ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/maps/map_corridor_protection_return_period.png`
- `ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/maps/map_no_loss_stream_context.png`
- `ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/tables/corridor_protection_level_all_corridor_municipalities.csv`
- `ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/tables/large_municipalities_portfolio_status.csv`
- `ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/tables/corridor_stream_context_by_municipality.csv`

## Provider clarification

The provider table contains all municipalities with modeled damage/loss in the provider flood portfolio. Municipalities absent from the table should currently be interpreted as municipalities without modeled losses in that portfolio.

Important clarification:

- No municipality was removed because of a minimum event-count threshold.
- There is currently no 30-event threshold in the provided table.
- The protection table is derived directly from the loss data without an additional event-count filter.

Therefore, absence from `elbe_protection_level_mun.csv` must not be interpreted as:

- no flood risk,
- no flood exposure,
- no flood protection,
- or a random missing value.

The safest current wording is:

`No modeled loss in the provider flood portfolio.`

## Join result with the EFAS/JRC RP500 corridor

Current thesis corridor:

- RP500 corridor municipalities: `835`

Provider protection/loss table:

- municipalities in `elbe_protection_level_mun.csv`: `301`
- protection-table municipalities matching the RP500 corridor: `280`
- RP500 corridor municipalities absent from the protection table: `555`
- protection-table municipalities outside the current RP500 corridor: `21`

The match rate within the RP500 corridor is therefore `280 / 835 = 33.5%`.

## What the maps add

The coverage map shows that the unmatched municipalities are not distributed randomly. The provider loss/protection portfolio is concentrated in a central and upper/middle Elbe-related corridor, while many municipalities in the lower Elbe / Hamburg / Schleswig-Holstein part of the EFAS corridor are coded as no modeled loss.

The city-labelled maps were added to avoid misreading the unmatched area as only low-population periphery. The large-municipality label rule is:

- all corridor municipalities with population `>= 50,000`, plus
- Elbe-near municipalities within `5 km` of the Elbe line with population `>= 25,000`.

Among the `33` labelled large or Elbe-near municipalities, `16` have no modeled loss in the provider portfolio. Examples include:

- Hamburg
- Magdeburg
- Kiel
- Potsdam
- Cottbus
- Norderstedt
- Lueneburg
- Elmshorn
- Cuxhaven
- Stade
- Seevetal
- Winsen (Luhe)
- Wedel
- Geesthacht

Several of these are directly on or near the Elbe and have non-trivial EFAS flood shares. Examples checked in the current table include:

- Hamburg: `flood_share_rp100 = 0.399`
- Magdeburg: `flood_share_rp100 = 0.367`
- Stade: `flood_share_rp100 = 0.357`
- Winsen (Luhe): `flood_share_rp100 = 0.425`
- Brunsbuettel: `flood_share_rp100 = 0.867`
- Jork: `flood_share_rp100 = 0.959`

This means the difference cannot be explained only by tiny or irrelevant municipalities.

## Stream-context diagnostic

The stream-context check was created in response to Nivedita's suggestion that absent municipalities might mainly be related to small streams not represented in the provider simulation.

Among the `555` corridor municipalities without modeled losses:

- `339` (`61.1%`) intersect only minor/other WFD water bodies.
- `16` (`2.9%`) do not intersect a WFD river body in the checked layer.
- Together, `355` (`64.0%`) are not directly associated with the Elbe main stem or a larger WFD river/tributary.
- `138` (`24.9%`) intersect a larger WFD river/tributary.
- `62` (`11.2%`) are directly on or within 1 km of the Elbe main stem.

Interpretation:

The small-stream explanation is partly supported, but it is not complete. Many absent municipalities plausibly fall outside the provider river/loss simulation logic because they relate to smaller streams or water bodies. However, a non-negligible group of absent municipalities lies directly on or close to the Elbe main stem, especially in the lower Elbe / tidal Elbe area, and some of these have high EFAS flood shares.

## Current thesis interpretation

The current working interpretation is a model-domain / portfolio-coverage mismatch between:

1. the thesis-defined EFAS/JRC RP500 exposure corridor, and
2. the provider flood loss/protection portfolio.

This should be treated as a methodological feature and limitation, not as a simple data error.

For the thesis, the defensible structure is:

- Main exposure-vulnerability analysis: all `835` EFAS/JRC RP500 corridor municipalities.
- Protection return-period analysis: only the `280` corridor municipalities with modeled loss/protection values.
- Remaining `555` municipalities: coded and discussed as `no modeled loss in provider flood portfolio`.

The protection table should therefore not replace the corridor definition. It is an additional loss/protection module with its own coverage boundary.

## Open questions after follow-up email

The following points remain open until Nivedita/Tan/Phillip reply:

1. Does the provider portfolio cover the lower Elbe / Hamburg / tidal or storm-surge-influenced section?
2. Is the absence of Hamburg/lower-Elbe municipalities mainly caused by model domain, river-network coverage, loss/exposure representation, or hazard-process definition?
3. Should the thesis describe the unmatched municipalities as outside the provider loss portfolio, no modeled loss, or both?
4. Is there provider-recommended wording for this limitation?
5. Will Tan's provider exposure-municipality list clarify which municipalities are in their exposure/loss domain?

## Chapter 4 consequence

Chapter 4 should include the protection data as a separate damage/loss-based module. It should explicitly state that protection return periods are available only for the matched loss-portfolio subset.

Suggested wording:

`Because the protection return-period table is derived from modeled losses, it is not available for every municipality in the EFAS/JRC RP500 corridor. Municipalities absent from the table are not interpreted as unprotected or risk-free. Following provider clarification, they are coded as municipalities with no modeled loss in the provider flood portfolio. Protection-return-period analyses are therefore restricted to the matched loss-portfolio subset, while the full corridor remains the basis for the exposure-vulnerability analysis.`

## Results and discussion consequence

The coverage pattern itself is a result worth reporting, but it should be framed carefully:

- It shows that the protection/loss module has a narrower effective coverage than the EFAS/JRC exposure corridor.
- It raises an important limitation for lower Elbe / Hamburg interpretation.
- It does not prove that those municipalities have no flood risk or no flood protection.
- It may reflect differences between riverine EFAS exposure, provider flood/loss portfolio, tidal/coastal processes, and represented assets/losses.

