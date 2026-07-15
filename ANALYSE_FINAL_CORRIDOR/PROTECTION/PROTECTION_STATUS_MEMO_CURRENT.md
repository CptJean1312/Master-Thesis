# Protection Status Memo - Current Working Interpretation

Date: 2026-07-13

Updated: 2026-07-15 after final provider clarification from Nivedita.

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

The provider table contains all municipalities with positive modeled damage/loss in the provider flood portfolio. Municipalities absent from the table should be interpreted as municipalities without simulated loss events in that portfolio.

Important clarification:

- No municipality was removed because of a minimum event-count threshold.
- There is currently no 30-event threshold in the provided table.
- The protection table is derived directly from the loss data without an additional event-count filter.
- In the final clarification, Nivedita recommended treating municipalities without simulated loss events as having modeled loss probability `0`.
- For larger cities such as Hamburg and Magdeburg, the suggested interpretation is that the protection level is high enough that no loss events occur in the simulation.
- For smaller streams, the limitation is that they are not considered in the provider simulation and this should be stated explicitly in the thesis.

Therefore, absence from `elbe_protection_level_mun.csv` must not be interpreted as:

- no flood risk,
- no flood exposure,
- no flood protection,
- or a random missing value.

The safest current wording is:

`No simulated loss event in the provider flood portfolio; modeled annual loss probability treated as zero.`

For finite protection return-period analysis, however, only municipalities with positive modeled loss probability and an estimated protection return period should be used. Municipalities with zero modeled loss probability should be treated as right-censored / no-event cases rather than assigned an arbitrary finite return period.

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

Among the `33` labelled large or Elbe-near municipalities, `16` have no simulated loss event in the provider portfolio. Examples include:

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

This means the difference cannot be explained only by tiny or irrelevant municipalities. After Nivedita's final clarification, these cases are best interpreted as protected/no-event cases within the provider model rather than as simple missing data.

## Stream-context diagnostic

The stream-context check was created in response to Nivedita's suggestion that absent municipalities might mainly be related to small streams not represented in the provider simulation.

Among the `555` corridor municipalities without simulated loss events:

- `339` (`61.1%`) intersect only minor/other WFD water bodies.
- `16` (`2.9%`) do not intersect a WFD river body in the checked layer.
- Together, `355` (`64.0%`) are not directly associated with the Elbe main stem or a larger WFD river/tributary.
- `138` (`24.9%`) intersect a larger WFD river/tributary.
- `62` (`11.2%`) are directly on or within 1 km of the Elbe main stem.

Interpretation:

The stream-context result should now be read together with the final provider clarification. Many absent municipalities plausibly fall outside the represented river/loss simulation logic because they relate to smaller streams or water bodies. This must be reported as a limitation. However, a non-negligible group of absent municipalities lies directly on or close to the Elbe main stem, including large cities. For these cases, the provider interpretation is that the protection level is high enough that no simulated loss events occur, implying modeled annual loss probability `0`.

## Current thesis interpretation

The current working interpretation is no longer a simple data-gap interpretation. Instead, the provider table should be read as a positive-loss table derived from a process-based flood-loss simulation:

1. municipalities in the table have simulated losses and finite protection return-period estimates;
2. municipalities absent from the table have no simulated losses and should be assigned modeled annual loss probability `0`;
3. the absence of smaller-stream municipalities partly reflects the fact that smaller streams are not represented in the provider simulation.

This should be treated as a methodological feature and limitation, not as a simple data error.

For the thesis, the defensible structure is:

- Main exposure-vulnerability analysis: all `835` EFAS/JRC RP500 corridor municipalities.
- Finite protection return-period analysis: only the `280` corridor municipalities with modeled loss/protection values.
- Modeled annual loss probability interpretation: the remaining `555` municipalities can be coded as `0` modeled annual loss probability / no simulated loss event.
- Limitations: smaller streams are not considered in the provider simulation; EFAS exposure and provider loss/protection outputs are not identical hazard/loss concepts.

The protection table should therefore not replace the corridor definition. It is an additional loss/protection module with its own coverage boundary.

## Resolved and remaining questions after final provider email

Resolved:

- For large cities such as Hamburg and Magdeburg, no loss event in the provider table can be interpreted as protection level high enough that no simulated loss event is observed.
- The modeled annual loss probability for municipalities without simulated loss events should be treated as `0`.
- Smaller streams are not considered and this should be written as a limitation.
- Sairam et al. (2021) should be cited for the simulation characteristics and limitations.

Still useful to check if Tan sends the exposure-municipality list:

1. Whether all no-event municipalities are inside the provider exposure domain but have zero loss, or whether some are absent because they are outside the represented river network.
2. Whether the lower Elbe / Hamburg section has any special tidal/coastal-process limitation beyond the general no-event interpretation.

## Chapter 4 consequence

Chapter 4 should include the protection data as a separate damage/loss-based module. It should explicitly state that finite protection return periods are available only for the positive-loss subset, while municipalities without loss events are assigned modeled annual loss probability `0`.

Suggested wording:

`Because the protection return-period table is derived from modeled losses, finite protection return periods are available only where at least one simulated loss event occurs. Municipalities absent from the table are not interpreted as unprotected or risk-free. Following provider clarification, they are treated as municipalities with no simulated loss event and modeled annual loss probability zero in the provider flood portfolio. Finite protection-return-period analyses are therefore restricted to the positive-loss subset, while the full corridor remains the basis for the exposure-vulnerability analysis. Smaller streams are not represented in the provider simulation and are reported as a limitation.`

## Results and discussion consequence

The coverage / no-event pattern itself is a result worth reporting, but it should be framed carefully:

- It shows where the provider simulation produces positive modeled losses and where it produces no simulated loss events.
- It supports the interpretation that large cities may have very high protection levels in the simulation.
- It does not prove that those municipalities have no flood risk or no flood protection.
- It also reflects a limitation: smaller streams are not considered in the provider simulation.
- It may reflect differences between riverine EFAS exposure, provider flood/loss portfolio, represented flood processes, and represented assets/losses.
