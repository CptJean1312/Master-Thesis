# Land-Cover and Figure Integration QA

Updated: 2026-07-21

## Status

The RP100 DLR land-cover extension is integrated into the active thesis draft and the active thesis figures have been rebuilt with one design system.

## Methodological Check

- The DLR source raster is a 10 m categorical land-cover product in `EPSG:3035`.
- The seven classes are converted to separate binary layers and projected to the 100 m EFAS/JRC template in `EPSG:25832` using area-weighted averaging.
- Class fractions are multiplied by cell area and intersected with the RP100 flood mask.
- All `835` corridor municipalities are retained.
- Aggregated land-cover areas reproduce municipal area within raster-processing precision.
- No calculated Artificial Land flood share lies outside `[0, 1]`.
- One municipality has no DLR Artificial Land and therefore no defined Artificial Land denominator.

This workflow is methodologically coherent for a corridor-wide robustness analysis. It preserves fractional class coverage when moving from the 10 m source to the 100 m flood grid and performs area calculations in the common projected CRS. The main interpretive boundary remains essential: Artificial Land is a built/artificial land-cover proxy, not population, buildings, residential use, or asset value.

## Empirical Check

- Classified RP100 flooded land: `8,355.2 km²`.
- Artificial Land: `499.4 km²` (`6.0%`).
- Vulnerability versus total-area RP100 exposure: Pearson `-0.050`, Spearman `-0.050`.
- Vulnerability versus RP100 Artificial Land exposure: Pearson `-0.088`, Spearman `-0.063`.
- Highest Artificial Land exposure tercile: `87` positive-loss and `191` no-event municipalities.

The land-cover result refines but does not reverse the thesis argument. It confirms that the total-area measure contains extensive non-built land, while also showing that no hidden positive vulnerability gradient appears when the denominator is restricted to Artificial Land. The provider loss/protection outcome remains analytically distinct from both exposure measures.

## Draft Integration

- Chapter 3 introduces the Artificial Land refinement and names Figure 3.1 in the prose.
- Chapter 4 documents the DLR source, processing, formula, group definitions, QA, statistical use, and reproducibility scripts.
- Chapter 5 includes a complete land-cover results subsection and Figures 5.7 to 5.9.
- Chapter 5 synthesis distinguishes total-area exposure, Artificial Land exposure, and provider loss/protection outcomes.
- Chapter 6 interprets the robustness result and documents its limits.
- Chapter 7 includes the result in the findings, contribution, outlook, and final remarks.

## Figure QA

- Embedded figure sequence: Figure 2.1, Figure 3.1, and Figures 5.1 to 5.17.
- Every embedded figure has a matching caption.
- Every embedded local image path exists.
- All final result graphics are available as 360 dpi PNG and vector PDF.
- Maps use a common font, color system, corridor geometry, Elbe line, federal-state context, major-city labels, north arrow, scale bar, legend, and source caption where applicable.
- Charts use the same typography, colors, legend position, and source-caption treatment.
- The conceptual model was rebuilt to show modeled protection/loss outcomes rather than a protection-infrastructure inventory.

Final figure outputs:

- PNG: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/THESIS_DRAFTS/FIGURES/FINAL/png`
- PDF: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/THESIS_DRAFTS/FIGURES/FINAL/pdf`
- Manifest: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/THESIS_DRAFTS/FIGURES/FINAL/FIGURE_MANIFEST.csv`
- Builder: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/BUILD_THESIS_FINAL_FIGURES.R`

## DLR Source

- Dataset page: https://geoservice.dlr.de/web/datasets/landcover
- Download directory: https://download.geoservice.dlr.de/LCC_DE/files/
- DOI: https://doi.org/10.15489/1ccmlap3mn39
- Recommended citation: German Aerospace Center (DLR). (2020). *Land Cover DE - Sentinel-2 - Germany, 2015*. https://doi.org/10.15489/1ccmlap3mn39

## Remaining Limits

- RP100 is the only return period processed for land cover.
- The DLR source is based on 2015 to 2017 imagery and does not perfectly match all other reference years.
- Municipality-level analysis cannot identify exposed households or within-municipality inequalities.
- Artificial Land combines residential, industrial, transport, and other artificial surfaces.
- Provider no-event status remains model-bound and must not be interpreted as real-world absence of flood risk.
