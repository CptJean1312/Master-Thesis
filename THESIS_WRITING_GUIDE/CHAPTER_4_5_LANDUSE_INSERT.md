# Chapter 4/5 Land-Use Exposure Insert

This file contains thesis-ready text blocks for adding the DLR land-use exposure extension to the draft. It should be inserted after the main EFAS/JRC municipal exposure method in Chapter 4 and after the main RP100 exposure-vulnerability results in Chapter 5.

## Chapter 4 Method Insert

To address the limitation that municipal flood shares measure flooded administrative area rather than settlement or asset exposure, an additional land-cover-sensitive exposure layer was constructed. This extension uses the DLR Land Cover DE 2015 raster as a spatial proxy for the type of land cover intersecting the modeled RP100 flood footprint. The dataset provides a 10 m categorical land-cover classification in `EPSG:3035` with seven classes: Artificial Land, Open Soil, High Seasonal Vegetation, High Perennial Vegetation, Low Seasonal Vegetation, Low Perennial Vegetation, and Water. For the present analysis, these classes were grouped into four thesis categories: `artificial`, `open_or_seasonal`, `perennial_vegetation`, and `water`. The `artificial` group is based on DLR class 1, Artificial Land, and is used as a built/artificial land-cover proxy. It is not interpreted as a cadastral building footprint or as a direct population-exposure measure.

The land-cover raster was harmonized with the existing flood-exposure workflow. The RP100 EFAS/JRC flood raster in `EPSG:25832` was used as the analysis template and cropped to the final RP500 corridor extent. Each DLR class was converted to a binary raster and reprojected to this corridor-cropped 100 m flood grid using area-weighted averaging. The resulting class-fraction rasters indicate the proportion of each 100 m analysis cell covered by a given DLR class. These fractions were multiplied by the cell area in `EPSG:25832` to obtain class-specific land-cover area per grid cell. Flooded class areas were then calculated by masking these class-area rasters with valid RP100 flood-depth cells and summing the result within municipality polygons. This produced, for each municipality, total land-cover area by class, flooded RP100 land-cover area by class, grouped land-cover exposure metrics, and the flooded share of artificial land.

This extension is used as a robustness and substantive refinement of the area-based RP100 exposure measure. It asks whether the weak association between socio-economic vulnerability and RP100 municipal flood share remains weak when exposure is measured only for artificial land-cover areas. The method reduces one important ambiguity of the total-area flood share, namely that high municipal exposure may reflect floodplain, agricultural, water, or low-density territory. At the same time, it does not fully solve the exposure problem. Artificial land cover includes multiple built or sealed surfaces and does not identify buildings, residents, assets, or social groups within the flood footprint.

Relevant outputs:

- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/tables/corridor_landuse_exposure_wide.csv`
- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/tables/corridor_landuse_protection_analysis_rp100.csv`
- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/gpkg/corridor_landuse_protection_analysis_rp100.gpkg`

## Chapter 5 Results Insert

The land-cover extension shows that the RP100 flooded area is not dominated by artificial land. Figure 5.x summarizes the land-cover composition of the modeled RP100 flooded area after aggregation to the four thesis groups. Artificial land accounts for `499.4 km²`, or `6.0%`, of all flooded land-cover area. By contrast, open or seasonal land cover accounts for `3,883.8 km²` (`46.5%`) and perennial vegetation for `3,289.5 km²` (`39.4%`). Water accounts for another `682.4 km²` (`8.2%`). This confirms that the original area-based flood share captures a broad floodplain geography rather than only built or settlement-related exposure.

![Figure 5.x Land-cover composition of RP100 flooded area](/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/maps/plot_rp100_flooded_landcover_composition.png)

*Figure 5.x. Land-cover composition of modeled RP100 flooded area in the Elbe corridor. Output path: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/maps/plot_rp100_flooded_landcover_composition.png`.*

The spatial pattern of artificial-land exposure is shown in Figure 5.x. The map displays, for each corridor municipality, the share of DLR Artificial Land intersecting the modeled RP100 flood footprint. High values occur in several municipalities along the Elbe corridor, but the pattern remains spatially heterogeneous. The highest artificial-land exposure values are not limited to the most socially vulnerable municipalities, and many municipalities with low or medium vulnerability also show substantial artificial-land exposure.

![Figure 5.x RP100 artificial-land exposure map](/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/maps/map_rp100_artificial_land_flood_share_context.png)

*Figure 5.x. RP100 artificial-land exposure in the RP500 Elbe corridor, measured as the flooded share of DLR Artificial Land in each municipality. Output path: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/maps/map_rp100_artificial_land_flood_share_context.png`.*

The correlation results support this interpretation. The relationship between the final vulnerability index and the original RP100 total-area flood share is weak and slightly negative (`Pearson = -0.050`, `Spearman = -0.050`). Replacing total municipal area exposure with artificial-land exposure does not produce a stronger positive justice gradient. The correlation between vulnerability and RP100 artificial-land flood share is also weak and slightly negative (`Pearson = -0.088`, `Spearman = -0.063`). Among the vulnerability dimensions, the deprivation and labour-market dimension has the most negative Pearson correlation with artificial-land exposure (`-0.165`), while the other dimensions remain close to zero.

Figure 5.x compares RP100 artificial-land exposure across vulnerability quintiles. The lowest-vulnerability quintile has the highest mean flooded share of artificial land (`26.4%`). The highest-vulnerability quintile has a mean of `16.7%`, while Q4 has the lowest mean (`11.9%`). Median values are much lower across all quintiles, ranging from `5.5%` to `7.5%`, which shows that high artificial-land exposure is driven by a smaller number of municipalities with very high values rather than by a broad shift across the whole vulnerability distribution.

![Figure 5.x Artificial-land exposure by vulnerability quintile](/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/plots/plot_rp100_artificial_land_exposure_by_vulnerability_quintile.png)

*Figure 5.x. RP100 artificial-land exposure by vulnerability quintile. Output path: `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/ANALYSE_FINAL_CORRIDOR/LANDUSE/outputs/plots/plot_rp100_artificial_land_exposure_by_vulnerability_quintile.png`.*

The land-use extension therefore refines rather than reverses the main exposure finding. The original RP100 municipal flood share clearly includes large non-artificial land-cover areas, so it should not be read as direct settlement or population exposure. However, when the analysis is restricted to artificial land cover, the relationship with socio-economic vulnerability remains weak and non-monotonic. The result strengthens the conclusion that the Elbe corridor does not show a simple pattern in which more vulnerable municipalities are generally more exposed at RP100, even when exposure is measured through a built/artificial land-cover proxy.

## Chapter 6 Discussion Insert

The DLR land-use extension partly addresses the strongest limitation of the area-based exposure metric. It shows that the main RP100 flood-share indicator does indeed combine very different types of land cover. Most modeled RP100 flooded land-cover area is not artificial land, but open/seasonal land cover, perennial vegetation, or water. This means that total municipal flood share should be interpreted as a spatial floodplain indicator rather than as a direct measure of people, buildings, or assets at risk.

At the same time, the extension also shows that correcting for land-cover type does not create the expected positive vulnerability gradient. Artificial-land exposure remains weakly and slightly negatively associated with the final vulnerability index. This is important because it prevents two opposite overinterpretations. The thesis should not claim that the original RP100 exposure result captures social exposure in a direct sense. But it should also not assume that a more settlement-related proxy would automatically reveal a hidden regressive pattern. Within the available municipality-level data, the weak exposure-vulnerability relationship remains robust.

The justice-relevant implication is that the empirical story still has to move beyond exposure alone. The land-use extension makes the exposure result more defensible, but it does not replace the protection/loss analysis. The clearer social pattern remains in the modeled loss/protection layer, where positive modeled losses and annual loss probability are more clearly associated with vulnerability than either total-area or artificial-land RP100 exposure.
