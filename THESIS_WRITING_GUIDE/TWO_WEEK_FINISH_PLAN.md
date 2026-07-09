# Two-Week Thesis Finish Plan

Created: 2026-07-09

Purpose: turn the current Elbe corridor thesis into a complete, defensible submission draft within two weeks.

Current methodological turning point: the protection dataset is now interpretable. According to provider clarification, `elbe_protection_level_mun.csv` contains all municipalities with modeled damage/loss in the provider flood portfolio. Municipalities absent from the table are not missing because of an event-count threshold; they do not experience modeled losses in that model. There is no 30-event filter in the delivered table.

---

# Current Stand

## Already solid

- RP500 Elbe flood corridor defined: `835` municipalities.
- EFAS/JRC multi-return-period exposure pipeline complete for `RP10`, `RP20`, `RP50`, `RP100`, `RP200`, `RP500`.
- RP100 exposure benchmark exists.
- Exposure quality checks exist.
- Exposure curves and typology exist.
- INKAR corridor inventory exists with `176` available indicators.
- Variable review and PCA comparison exist.
- Chapter 1-3 exist as integrated draft text.
- Chapter 4 exists as a strong methods draft with placeholders.
- Protection data are now available and interpretable.

## Still to lock

- Final vulnerability index specification.
- Final integrated analysis table combining exposure, vulnerability, exposure curves, and protection.
- Final results figures/tables.
- Chapter 4 update for direct protection data.
- Chapter 5 Results.
- Chapter 6 Discussion.
- Chapter 7 Conclusion.
- Figure 3.1 final polish.
- Final consistency pass across terminology, citations, captions, and limitations.

---

# Week 1: Freeze Methods and Produce Final Results

## Day 1: Protection data integration

Goal: make the protection module final and reproducible.

Tasks:
- keep raw `elbe_protection_level_mun.csv` unchanged;
- use `ANALYSE_FINAL_CORRIDOR/PROTECTION/BUILD_PROTECTION_MAPS.R` as the reproducible protection script;
- document provider clarification in `PROTECTION_DATA_PROVIDER_CLARIFICATION.md`;
- compare Tan's forthcoming exposure municipality list against the RP500 corridor once received;
- decide final wording for absent municipalities: `No modeled loss in provider portfolio`;
- keep maps:
  - `map_corridor_protection_coverage.png`
  - `map_corridor_protection_return_period.png`

Deliverable:
- final protection data subsection notes for Chapter 4.

## Day 2: Final vulnerability index decision

Goal: stop reopening the index question.

Recommended decision:
- main index: `thesis_candidate_23`;
- sensitivity checks: `original_51`, `student_sensitivity_26`, and optionally `curated_17`;
- do not use all 176 as final index.

Tasks:
- rerun PCA/main analysis using `thesis_candidate_23`;
- export final PCA diagnostics:
  - variable list;
  - scree table;
  - top loadings;
  - component interpretation;
  - index direction check;
- create final vulnerability map.

Deliverable:
- final vulnerability index table and diagnostics.

## Day 3: Integrated final analysis table

Goal: create one master table for all Results.

Table should include:
- `AGS`, municipality name;
- municipality area;
- flood shares for all return periods;
- final vulnerability index and sensitivity indices;
- exposure curve metrics;
- protection availability / portfolio status;
- protection return period;
- annual loss probability;
- key classes:
  - vulnerability tertile/quintile;
  - RP100 exposure class;
  - protection RP class;
  - loss-portfolio status.

Deliverable:
- one final CSV/GPKG for Chapter 5.

## Day 4: Results figures and tables

Goal: make the visual and tabular backbone of Chapter 5.

Minimum figure set:
- study area / corridor map;
- RP100 exposure map;
- vulnerability index map;
- protection coverage map;
- protection return-period map;
- exposure vs vulnerability plot;
- protection return period vs exposure/vulnerability plot for matched municipalities;
- grouped bar/table for vulnerability-exposure-protection overlap.

Minimum table set:
- study sample summary;
- exposure summary by return period;
- PCA/index diagnostics;
- protection coverage and return-period summary;
- key distributive comparison table.

Deliverable:
- all Chapter 5 figures/tables exported.

## Day 5: Update Chapter 4

Goal: methods chapter matches what was actually done.

Tasks:
- add protection dataset as core data source;
- replace old placeholder language in section 4.8;
- explain provider interpretation:
  - included = modeled loss portfolio with estimated protection RP;
  - absent = no modeled loss in provider portfolio;
  - no event-count filter applied;
- update statistical/spatial analysis section with final model/table structure;
- update PCA section after final rerun.

Deliverable:
- Chapter 4 no longer has unresolved methodological placeholders except source/citation details.

## Day 6: Write Chapter 5 Results

Goal: complete first full Results draft.

Recommended structure:
- 5.1 Study corridor and exposure profile;
- 5.2 Socio-economic vulnerability patterns;
- 5.3 Protection portfolio coverage and return-period patterns;
- 5.4 Distributional overlap: exposure, vulnerability, and protection/loss portfolio;
- 5.5 Sensitivity checks.

Rule:
- describe results first, interpret cautiously, save broader meaning for Discussion.

Deliverable:
- Chapter 5 complete as rough but coherent text.

## Day 7: Buffer and internal consistency check

Goal: catch broken joins, inconsistent counts, and unstable wording.

Tasks:
- verify all n-values:
  - corridor `835`;
  - INKAR matched `834`;
  - protection table `301`;
  - protection matched to corridor `280`;
  - no-loss/absent in corridor `555`;
- check all figures have captions;
- check Chapter 4 and 5 use the same terminology.

Deliverable:
- Week-1 complete empirical backbone.

---

# Week 2: Write, Polish, Defend

## Day 8: Discussion core argument

Goal: write the analytical heart of the thesis.

Recommended structure:
- what the Elbe corridor results show;
- what the protection data changes;
- how exposure, vulnerability, and modeled loss/protection relate;
- what is justice-relevant and what cannot be claimed;
- relation to German flood-risk management and flood justice literature.

Deliverable:
- Chapter 6 first full draft.

## Day 9: Limitations and robustness

Goal: make the thesis defensible.

Must include:
- municipality-level MAUP/ecological inference limits;
- EFAS/JRC modeled hazard not legal flood maps;
- area-based exposure, not population/damage exposure;
- protection dataset represents modeled loss portfolio, not observed protection infrastructure;
- absence from protection table means no modeled loss in provider portfolio;
- PCA index sensitivity and variable-selection limits;
- no causal policy-intent claims.

Deliverable:
- strong limitations subsection.

## Day 10: Conclusion

Goal: close the thesis cleanly.

Include:
- answer RQ1 and RQ2;
- summarize methodological contribution;
- state policy/planning relevance;
- propose future work:
  - settlement/population-weighted exposure;
  - deeper protection infrastructure validation;
  - comparison with provider exposure list;
  - temporal/climate scenarios.

Deliverable:
- Chapter 7 full draft.

## Day 11: Front-half polish

Goal: make Chapters 1-3 align with the final empirical story.

Tasks:
- update Introduction if protection is now more concrete;
- update Chapter 3 scope/assumptions;
- ensure no outdated language says protection is only future/inferred;
- finalize Figure 3.1 caption and placement.

Deliverable:
- Chapters 1-3 consistent with final method.

## Day 12: Citation and bibliography pass

Goal: reduce citation risk.

Tasks:
- verify EFAS/JRC source citation;
- verify INKAR/BBSR citation;
- verify VG250/BKG citation;
- cite provider protection data as personal communication / supplied dataset if no formal citation exists;
- check all in-text citations have reference-list entries.

Deliverable:
- citation-clean draft.

## Day 13: Full read-through and formatting

Goal: make the draft readable as one thesis.

Tasks:
- remove duplicated paragraphs;
- standardize terms:
  - `RP500 flood corridor`;
  - `RP100 benchmark exposure`;
  - `socio-economic vulnerability`;
  - `modeled damage/loss portfolio`;
  - `protection return period`;
- check figure/table numbering;
- check chapter transitions.

Deliverable:
- near-final manuscript.

## Day 14: Final supervisor/sendable package

Goal: prepare a version that can go to supervisors or into final revision.

Package:
- thesis draft PDF or DOCX;
- key maps/figures;
- final analysis summary table;
- short memo explaining protection-data interpretation and remaining limits.

Deliverable:
- complete supervisor-ready thesis package.

---

# Immediate Next Three Moves

1. Integrate Tan's forthcoming exposure municipality list and compare it to the RP500 corridor.
2. Finalize and rerun the vulnerability index using `thesis_candidate_23`.
3. Build the final integrated Results table and start Chapter 5.
