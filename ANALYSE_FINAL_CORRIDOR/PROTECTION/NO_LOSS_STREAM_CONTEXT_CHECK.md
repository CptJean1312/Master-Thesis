# No-Loss Stream Context Check

Date: 2026-07-13

Purpose: respond to Nivedita's question whether the RP500 corridor municipalities absent from `elbe_protection_level_mun.csv` mainly have only small streams flowing through them.

## Data used

- Corridor municipalities: RP500 EFAS/JRC corridor, `n = 835`
- Provider protection table: `elbe_protection_level_mun.csv`
- Elbe main-stem layer: `Elbe.gpkg`
- WFD river water bodies: `am_riverwaterbody_de-basin.gpkg`

## Classification logic

For each corridor municipality, the script checks:

1. Whether it is directly on or within 1 km of the Elbe main-stem line.
2. Whether it intersects a larger WFD river/tributary.
3. Whether it intersects only minor/other WFD river water bodies.
4. Whether no WFD river water body intersects the municipality.

This is a GIS diagnostic and not a hydraulic attribution of the modeled flood source.

## Main result

The small-stream explanation is partly supported, but it does not explain all absent municipalities.

Among the `555` corridor municipalities without modeled losses in the provider portfolio:

- `339` (`61.1%`) intersect only minor/other WFD water bodies.
- `16` (`2.9%`) do not intersect a WFD river water body in the checked layer.
- Together, `355` (`64.0%`) are not directly associated with the Elbe main stem or a larger WFD river/tributary.
- `138` (`24.9%`) intersect a larger WFD river/tributary.
- `62` (`11.2%`) are directly on or within 1 km of the Elbe main stem.

For comparison, among the `280` corridor municipalities with protection/loss values:

- `65` (`23.2%`) are directly on or within 1 km of the Elbe main stem.
- `136` (`48.6%`) intersect a larger WFD river/tributary.
- `78` (`27.9%`) intersect only minor/other WFD water bodies.
- `1` (`0.4%`) does not intersect a WFD river water body in the checked layer.

## Interpretation

The provider table is much more concentrated along the Elbe main stem and larger WFD river/tributary network than the EFAS RP500 corridor as a whole. This supports the idea that many absent municipalities may be related to smaller streams or water bodies not represented in the provider simulation.

However, the explanation is not complete. A non-negligible group of absent municipalities lies directly on or close to the Elbe main stem, especially in the lower Elbe / tidal Elbe area, and some of these have high EFAS flood shares. Therefore, the difference between the EFAS RP500 corridor and the provider loss portfolio likely also reflects differences in hazard model extent, river-network coverage, exposure/loss representation, or portfolio definition.

## Useful output files

- Map: `ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/maps/map_no_loss_stream_context.png`
- Municipality table: `ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/tables/corridor_stream_context_by_municipality.csv`
- Summary by portfolio status: `ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/tables/corridor_stream_context_summary_by_portfolio_status.csv`
- No-loss summary: `ANALYSE_FINAL_CORRIDOR/PROTECTION/outputs/tables/no_loss_stream_context_summary.csv`

## Suggested wording to Nivedita / Tan

Hi Nivedita, hi Tan,

I checked this spatially against the Elbe main-stem layer and the WFD river water-body layer.

The small-stream explanation is partly supported, but it does not explain all excluded municipalities. Of the 555 RP500 corridor municipalities without modeled losses in the provider portfolio, 339 (61.1%) intersect only minor/other WFD water bodies and another 16 (2.9%) do not intersect a WFD river body in the checked layer. So about 64% are not directly associated with the Elbe main stem or a larger WFD river/tributary.

However, 138 municipalities (24.9%) intersect a larger WFD river/tributary, and 62 municipalities (11.2%) are directly on or within 1 km of the Elbe main stem. Several of the latter are in the lower Elbe / tidal Elbe area and have high EFAS RP100/RP500 flood shares.

So my interpretation would be: many of the absent municipalities plausibly relate to smaller streams or water bodies not represented in the provider simulation, but the difference is not only a small-stream issue. It likely also reflects differences between my EFAS RP500 corridor definition and the provider flood/loss portfolio, especially in the lower Elbe area. Tan's exposure-municipality list will be very useful to compare these extents directly.

Best,
Max
