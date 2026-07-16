# Was Muss Gemacht Werden

Updated: 2026-07-16

Purpose: turn the current Elbe flood-justice thesis into a very strong, defensible submission draft. The realistic working target is to have a near-final thesis by the end of next week, but the quality target is higher than speed alone. If one methodological upgrade genuinely improves the thesis, it is worth doing.

Current strategic shift:

- The thesis draft already contains Chapters 1 to 7 as full prose.
- The main empirical story exists: RP100 area exposure is spatially uneven but not strongly socially regressive; modeled protection/loss outcomes show a clearer, but still cautious, social differentiation.
- The strongest remaining methodological critique is not the vulnerability index anymore. It is the difference between area-based flood share and socially/settlement-relevant exposure.
- The DLR Land Cover DE dataset has now been moved out of the Git repository and should be used as an external input from `DATEN + GIS`.

---

## Decision Rule

### Must go in

Something must go in if it changes one of these:

- the empirical answer to a research question;
- the defensibility of the protection/exposure interpretation;
- the ability to answer a hard reviewer;
- the formal submission quality.

### Should go in

Something should go in if it strengthens the argument without reopening the whole thesis.

### Cherry on top

Something is a cherry on top if it would be nice, but the thesis remains strong without it.

### Do not do now

Something should not be done now if it creates a new thesis instead of finishing this one.

---

# Priority 0: Hard Blockers Before Submission

These are not optional. If they are not done, the thesis is not submission-ready.

## P0.1 Bibliography and citation integrity

Status: open

Checklist:

- [ ] Build final bibliography/reference list.
- [ ] Ensure every in-text citation appears in the bibliography.
- [ ] Ensure every dataset has a citation or clear source note.
- [ ] Add DLR Land Cover DE citation if land-use analysis is included.
- [ ] Add BBSR/INKAR citation.
- [ ] Add BKG/VG250 citation.
- [ ] Add EFAS/JRC/Copernicus/JRC flood hazard map citations.
- [ ] Add provider protection data as supplied dataset/personal communication.
- [ ] Check spelling variants in citations, especially author names and years.

Why it matters:

This is a formal grading risk. A strong thesis can lose credibility very quickly if the references are incomplete.

## P0.2 Formal thesis apparatus

Status: open

Checklist:

- [ ] Replace manual placeholder table of contents with Word-generated TOC.
- [ ] Add figure list and table list in Word.
- [ ] Check figure numbering after final edits.
- [ ] Check captions and output paths.
- [ ] Remove output paths from final Word version if the supervisor prefers clean captions.
- [ ] Check title page names, date, reviewers, institute.
- [ ] Check page numbers, margins, line spacing, heading styles.

Why it matters:

This is not intellectual, but it affects whether the work feels finished.

## P0.3 No large raw data in Git

Status: partly done

Checklist:

- [x]irm DLR Land Cover DE is outside the Git repository.
- [x] Add or check `.gitignore` for large raster/zip patterns if needed.
- [x]Keep only scripts, metadata notes, and derived small outputs in the repo.
- [x] Document external raw-data path in the land-use script.

Why it matters:

The DLR raster is too large for normal Git workflow. The repo should remain reproducible but not contain giant raw data.

---

# Priority 1: Must-Have Intellectual and Methodological Work

These are the tasks that most directly move the thesis toward a very good grade.

## P1.1 Add DLR land-use exposure analysis

Status: urgent, must do

Core idea:

The current exposure metric is `flooded municipal area / total municipal area`. This is defensible but vulnerable to the criticism that high flood share may mainly represent floodplains, agricultural land, water, or low-density territory. The DLR land-cover layer lets us test whether flooded area also intersects artificial/built land-cover classes.

Input:

- DLR Land Cover DE 2015.
- External raw data path under `DATEN + GIS`.
- Local metadata:
  - 10 m x 10 m raster.
  - EPSG:3035.
  - class 1 = Artificial Land.
  - class 2 = Open Soil.
  - class 3 = High Seasonal Vegetation.
  - class 4 = High Perennial Vegetation.
  - class 5 = Low Seasonal Vegetation.
  - class 6 = Low Perennial Vegetation.
  - class 7 = Water.

Outputs to build:

- [ ] Land-cover area per municipality.
- [ ] Flooded land-cover area per municipality for RP100.
- [ ] If feasible, flooded land-cover area per municipality for RP10, RP20, RP50, RP100, RP200, RP500.
- [ ] `artificial_total_m2`.
- [ ] `artificial_flooded_rp100_m2`.
- [ ] `artificial_flood_share_rp100 = artificial_flooded_rp100_m2 / artificial_total_m2`.
- [ ] `flooded_area_artificial_share_rp100 = artificial_flooded_rp100_m2 / total_flooded_area_rp100_m2`.
- [ ] grouped land-cover categories:
  - `artificial`;
  - `open_or_seasonal`;
  - `perennial_vegetation`;
  - `water`.
- [ ] QA table with missing/no-artificial-land municipalities.

Main question:

Does the weak exposure-vulnerability result remain weak when exposure is measured through artificial/built land cover rather than total municipal area?

Expected thesis use:

- Chapter 4: add land-cover data and overlay method.
- Chapter 5: add results figure/table.
- Chapter 6: use it to answer the strongest exposure-metric critique.
- Chapter 7: update contribution and future work accordingly.

Decision after outputs:

- [ ] If the land-use result materially changes the story, include as a main Results subsection.
- [ ] If it confirms the existing story, include as a robustness/substantive refinement.
- [ ] If data quality is problematic, include only as limitation/future work and do not force it.

## P1.2 Rebuild or extend the master analysis table

Status: must do after land-use outputs

Checklist:

- [ ] Add land-use exposure indicators to the final corridor analysis table.
- [ ] Keep original RP100 flood share for comparability.
- [ ] Add artificial-land exposure metrics.
- [ ] Add grouped land-cover exposure metrics.
- [ ] Add vulnerability index and vulnerability dimensions.
- [ ] Add protection/loss status.
- [ ] Add no-event status.
- [ ] Add large-city labels if used in maps.

Required output:

- [ ] Updated CSV.
- [ ] Updated GPKG.
- [ ] Short data dictionary.

Why it matters:

The thesis needs one clean empirical backbone. No scattered tables that tell slightly different stories.

## P1.3 Rerun key Results with land-use exposure

Status: must do after P1.1 and P1.2

Checklist:

- [ ] Correlation: vulnerability index vs RP100 total-area flood share.
- [ ] Correlation: vulnerability index vs RP100 artificial-land flood share.
- [ ] Correlation: vulnerability dimensions vs artificial-land flood share.
- [ ] Quintile comparison using artificial-land flood share.
- [ ] High-high classification using vulnerability and artificial-land flood share.
- [ ] Compare high total-area exposure vs high artificial-land exposure.
- [ ] Check whether protection/loss positive cases align more strongly with artificial-land exposure than total-area exposure.

Minimum outputs:

- [ ] Scatterplot vulnerability vs artificial-land flood share.
- [ ] Boxplot artificial-land flood share by vulnerability quintile.
- [ ] Bivariate map: artificial-land exposure and vulnerability.
- [ ] Table comparing area-based and artificial-land-based results.

Why it matters:

This is the cleanest answer to the reviewer critique that the exposure metric may be socially weak because it includes fields, floodplain, water, and unused land.

## P1.4 Map redesign and visual consistency

Status: must do

Problem:

The existing maps are useful, but they do not yet feel like one thesis figure system. Orientation is sometimes hard without cities and federal states.

Checklist:

- [ ] Create one shared map theme for all maps.
- [ ] Use consistent fonts, legend positions, line widths, and color logic.
- [ ] Add Elbe river line where relevant.
- [ ] Add Bundesländer boundaries to the overview map.
- [ ] Add labels for major cities where useful:
  - Hamburg;
  - Magdeburg;
  - Dresden;
  - Leipzig;
  - Berlin;
  - possibly Wittenberge, Dessau-Roßlau, Lutherstadt Wittenberg, Riesa, Torgau, Cuxhaven depending on map extent.
- [ ] Avoid clutter: large cities on all maps only if they improve orientation.
- [ ] Use the same corridor outline and background style.
- [ ] Re-export all result maps at consistent resolution.

Minimum maps to standardize:

- [ ] Study area/corridor overview.
- [ ] RP100 exposure map.
- [ ] Vulnerability index map.
- [ ] Bivariate exposure-vulnerability map.
- [ ] Protection coverage map.
- [ ] Protection return-period map.
- [ ] Modeled annual loss probability map.
- [ ] New land-use/artificial exposure map if included.

Why it matters:

Maps are not decoration in this thesis. They carry the geographical argument. If they look inconsistent, the whole thesis feels less controlled.

## P1.5 Protection/no-event robustness and wording

Status: must do

Checklist:

- [ ] Keep no-event municipalities as `no simulated loss event`, not as real-world risk-free.
- [ ] Keep zero annual loss probability wording strictly model-bound.
- [ ] Add sensitivity table:
  - full corridor with no-event as zero modeled annual loss probability;
  - positive-loss subset only;
  - no-event municipalities discussed separately.
- [ ] Check whether no-event municipalities with high artificial-land exposure exist.
- [ ] If they exist, use them as important limitation/model-visibility cases.
- [ ] Keep provider clarification in Methods and Discussion.

Why it matters:

This is the second strongest reviewer attack after the area-exposure critique. The thesis can survive it if the wording is disciplined and the sensitivity view is explicit.

## P1.6 Title and research question alignment

Status: must discuss with supervisors

Current issue:

The current title uses `flood protection infrastructure`, but the empirical analysis measures modeled protection/loss outcomes, not an inventory of dikes, walls, or flood protection assets.

Checklist:

- [ ] Ask supervisors whether title can be changed.
- [ ] Prepare two or three alternative titles.
- [ ] If title changes, update Introduction, Research Questions, Chapter 4, Chapter 6, and Conclusion.
- [ ] If title cannot change, keep strong operationalization paragraph and do not overclaim.

Possible direction for a revised title:

- `Spatial distribution of modeled flood exposure, protection outcomes and socio-economic vulnerability in the Elbe flood corridor`
- `Flood exposure, modeled protection outcomes and socio-economic vulnerability in the Elbe flood corridor`
- `Distributive flood justice in the Elbe flood corridor: modeled exposure, protection outcomes and socio-economic vulnerability`

Why it matters:

This is not just wording. It affects what the reader expects the data to prove.

## P1.7 Strengthen the justice and capitalism-critical framing, but stay data-bound

Status: must do carefully

Aim:

The thesis should not become a neutral GIS exercise. It should show that flood exposure, protection, model visibility, land values, infrastructure priorities, and socio-economic vulnerability are politically and spatially produced. But it must not claim to have measured capital flows, investment decisions, or class relations directly.

Checklist:

- [ ] Keep distributive justice as the main empirical lens.
- [ ] Add or sharpen a short capitalism-critical/geographical reading in Chapter 2.
- [ ] Connect protection to:
  - land values;
  - asset concentration;
  - political visibility;
  - infrastructure investment priorities;
  - model visibility;
  - uneven capacity to attract attention and resources.
- [ ] In Chapter 6, interpret results as a question of uneven territorial visibility and protection-related outcomes.
- [ ] Avoid unsupported claims such as:
  - "capital intentionally protects rich municipalities";
  - "policy deliberately disadvantages vulnerable municipalities";
  - "no-event cities are safe";
  - "finite return period equals actual dike height".
- [ ] If adding new theory sources, keep them few and directly useful.

Why it matters:

This is part of the intellectual identity of the thesis. It can make the work stronger, but only if it remains tied to what the data can actually show.

---

# Priority 2: Should-Have Improvements

These are important, but they should not delay the thesis indefinitely.

## P2.1 Vulnerability dimensions as separate checks

Status: should do if quick

Checklist:

- [ ] Run correlations between each vulnerability dimension and:
  - RP100 total-area flood share;
  - RP100 artificial-land flood share;
  - annual loss probability;
  - finite protection return period among positive-loss municipalities.
- [ ] Check whether one dimension drives the justice signal more than the overall index.
- [ ] Add one compact table if useful.
- [ ] Do not over-expand Results.

Why it matters:

The dimensions are weakly correlated. A reviewer may ask whether the composite index hides specific patterns.

## P2.2 High-high municipality table

Status: should do

Checklist:

- [ ] Identify municipalities with:
  - high vulnerability;
  - high total-area flood share;
  - high artificial-land flood share;
  - positive modeled loss/protection concern.
- [ ] Export a table of top cases.
- [ ] Use as a practical interpretation, not as the whole thesis.

Why it matters:

This gives the thesis a concrete geographical payoff. It also helps supervisors and readers see why the analysis matters.

## P2.3 Figure 3.1 final polish

Status: should do

Checklist:

- [ ] Check whether Figure 3.1 still matches final design after land-use addition.
- [ ] Do not turn it into a workflow diagram.
- [ ] Keep exposure, vulnerability, protection/loss outcomes, and distributive justice logic visible.
- [ ] Export clean PNG for Word.

Why it matters:

This is the conceptual anchor. It should not look weaker than the empirical maps.

## P2.4 Chapter 4 update after land-use module

Status: should do immediately after land-use outputs

Checklist:

- [ ] Add DLR Land Cover DE to data sources.
- [ ] Explain class grouping.
- [ ] Explain CRS harmonization from EPSG:3035 to analysis CRS or vice versa.
- [ ] Explain overlay method.
- [ ] Explain why `artificial land` is used as built/settlement-related exposure proxy.
- [ ] State limitations:
  - 2015 land-cover timing vs later exposure/vulnerability data;
  - land-cover class is not cadastral building footprint;
  - artificial land includes roads/industrial/urban surfaces, not only residential buildings.

Why it matters:

If the land-use result goes into Results, Methods must make it reproducible.

## P2.5 Chapter 5 to 7 update after new outputs

Status: should do

Checklist:

- [ ] Insert land-use result where it logically belongs in Chapter 5.
- [ ] Update Discussion to reflect whether artificial-land exposure changes the story.
- [ ] Update Limitations so area-based exposure critique is handled as partly addressed, not only acknowledged.
- [ ] Update Conclusion with the final refined finding.

Why it matters:

New analysis cannot just be appended. It has to change the story where appropriate.

## P2.6 Supervisor-ready memo

Status: should do

Checklist:

- [ ] Prepare a 1-2 page memo with:
  - final research questions;
  - current key findings;
  - protection-data interpretation;
  - land-use exposure extension;
  - title question;
  - remaining limitations.
- [ ] Attach 3-5 key maps/figures.

Why it matters:

This helps get supervisor feedback fast without forcing them to read the full draft immediately.

---

# Priority 3: Cherry on Top

These would be nice, but only after P0-P2 are under control.

## P3.1 Alternative settlement dataset comparison

Examples:

- DLR World Settlement Footprint 2019.
- Copernicus Imperviousness.
- CORINE Land Cover.

Decision:

- [ ] Do only if DLR Land Cover DE gives unclear or suspicious results.
- [ ] Do not spend a week comparing land-cover products unless the thesis depends on it.

## P3.2 More advanced statistical modeling

Examples:

- Regression models.
- Spatial autocorrelation models.
- Controls for municipality size, population density, artificial-land share.

Decision:

- [ ] Only add if simple descriptive results become too ambiguous.
- [ ] Do not turn the thesis into a modeling paper at the last minute.

## P3.3 More Marxist theory expansion

Possible sources/angles:

- uneven development;
- urban political ecology;
- infrastructure and state spatiality;
- risk, land value, and protection as territorial investment.

Decision:

- [ ] Add only a compact, useful theoretical bridge.
- [ ] Do not rewrite Chapter 2 into a theory thesis.
- [ ] Keep claims attached to the data.

## P3.4 Full official hazard-map comparison

Decision:

- [ ] Too large for now unless data are already ready.
- [ ] Keep as future research/validation.

## P3.5 Interactive maps or appendix atlas

Decision:

- [ ] Nice for presentation.
- [ ] Not necessary for thesis grade unless all main work is finished.

---

# What Should Not Be Done Now

- [ ] Do not restart the vulnerability index from scratch.
- [ ] Do not chase a perfect physical flood-protection infrastructure inventory unless supervisors provide one.
- [ ] Do not convert the study back to full Elbe basin.
- [ ] Do not make causal claims about policy intent.
- [ ] Do not add a large theory section that the empirical analysis cannot carry.
- [ ] Do not include every possible map just because it exists.
- [ ] Do not let formatting wait until the last night.

---

# Practical Work Order From 2026-07-16

## Step 1: Land-use module

Target: 2026-07-16 to 2026-07-18

Checklist:

- [ ] Locate final external DLR raw-data path in `DATEN + GIS`.
- [ ] Write reproducible R script in `ANALYSE_FINAL_CORRIDOR/LANDUSE`.
- [ ] Create municipality-level land-cover totals.
- [ ] Create flooded land-cover totals for RP100.
- [ ] Extend to all return periods if runtime is acceptable.
- [ ] Export CSV/GPKG.
- [ ] Make first diagnostic maps/plots.
- [ ] Decide how strongly land-use results enter the thesis.

## Step 2: Rebuild figures and maps

Target: 2026-07-18 to 2026-07-20

Checklist:

- [ ] Create shared map theme.
- [ ] Add Bundesländer to overview map.
- [ ] Add major city labels.
- [ ] Re-export study area map.
- [ ] Re-export exposure/vulnerability/protection maps with consistent style.
- [ ] Add land-use/artificial exposure map if results support it.

## Step 3: Update empirical story

Target: 2026-07-20 to 2026-07-22

Checklist:

- [ ] Update Chapter 4 methods.
- [ ] Update Chapter 5 results.
- [ ] Update Chapter 6 discussion.
- [ ] Update Chapter 7 conclusion.
- [ ] Make sure all numbers match output tables.
- [ ] Keep the final argument cautious but confident.

## Step 4: Title and supervisor package

Target: 2026-07-22 to 2026-07-23

Checklist:

- [ ] Draft title options.
- [ ] Draft supervisor memo.
- [ ] Include key maps.
- [ ] Ask directly whether title change is acceptable.
- [ ] Ask whether land-use exposure extension is sufficient to address area-exposure critique.

## Step 5: Formal finalization

Target: 2026-07-23 to 2026-07-26

Checklist:

- [ ] Build bibliography.
- [ ] Convert/import into Word.
- [ ] Generate TOC, figure list, table list.
- [ ] Check all captions.
- [ ] Final language pass in Word.
- [ ] Final PDF export.
- [ ] One cold read from start to finish.

---

# Current Best Thesis Story

The thesis should not claim that the Elbe corridor shows one simple environmental-justice gradient. That would be too easy and not supported by the current results.

The stronger story is:

1. Flood exposure is spatially uneven.
2. Area-based exposure alone is not strongly socially regressive at municipality level.
3. This may partly reflect the fact that area exposure includes floodplains, agricultural land, water, and low-density land.
4. Land-cover-weighted exposure tests whether the socially relevant part of exposure is hidden inside the area metric.
5. Modeled protection/loss outcomes show a clearer social differentiation than area exposure alone.
6. No-event municipalities reveal both high modeled protection and model-visibility limits.
7. Flood justice is therefore not one map of where water goes, but a question of how exposure, land use, vulnerability, protection, loss, and model visibility are brought into relation.

This is a stronger, more mature argument than forcing a simple "poor municipalities are more flooded" claim.

---

# Immediate Next Move

Do this next:

- [ ] Build the DLR land-use exposure script.
- [ ] Produce the first land-use exposure table and a quick diagnostic summary.
- [ ] Decide from the actual output whether land-use exposure becomes a main Results block or a robustness block.

