# Full Draft Review Audit

Date: 2026-07-15

Active draft:

- `/Users/maxi_161/Desktop/UNI/Master/THESIS/Master-Thesis/THESIS_DRAFTS/THESIS DRAFT.md`

Purpose:

This memo records a hard reviewer pass over the current thesis draft. The review focused on thesis logic, evidence, citation support, Results consistency, figure integration, language, overclaim control, and the justice framing.

## Reviewer Verdict

The thesis logic is defensible. The current empirical story is coherent and supported by the available outputs:

- RP100 municipal area exposure is uneven, but not simply regressive across the vulnerability distribution.
- The protection/loss portfolio shows a clearer social differentiation than the area-exposure layer alone.
- The no-event class is meaningful, but must remain model-dependent and limitation-sensitive.
- The vulnerability index is methodologically defensible because it is dimension-balanced and uses PCA as a diagnostic rather than as a mechanical weighting engine.

The draft is not yet submission-ready. The most important remaining problems are not the Results logic. They are final-form problems: bibliography, title alignment, table of contents, figure duplication, and final citation formatting.

## Fixes Applied In The Draft

The following issues were corrected directly in the active draft:

- Fixed incorrect citation punctuation: `(Sairam et al., 2021) & (Thieken et al., 2016)` is now `(Sairam et al., 2021; Thieken et al., 2016)`.
- Corrected `social hirarchies` to `social hierarchies`.
- Replaced the questionable `Alfieri et al. (2016)` scope sentence with a more accurate EFAS/JRC hazard-map citation chain: `Alfieri et al. (2014); Baugh et al. (2024); Dottori et al. (2022)`.
- Corrected Fekete citations from split semicolon form to `Fekete (2009, 2010)`.
- Replaced a too-strong Results synthesis phrase, `social gradient`, with `clearer social differentiation`.
- Standardized visible spelling inconsistencies around `modeling`, `standardized`, `defense`, and `analyzing`.
- Added a critical flood-justice paragraph in Chapter 2 that connects protection to land values, infrastructure, political visibility, investment priorities, and model visibility, while explicitly avoiding unsupported causal claims.

## Results Check

Chapter 5 remains consistent with the verified Results audit:

- Corridor municipalities: `835`
- Vulnerability index records: `834`
- Missing INKAR municipality: `16076094` Berga-Wünschendorf
- Provider table total: `301`
- Provider positive-loss municipalities in corridor: `280`
- No simulated loss event / zero modeled annual loss probability: `555`
- RP100 exposure-vulnerability correlation: Pearson `-0.050`, Spearman `-0.050`
- RP100 mean exposure by vulnerability quintile: Q1 `35.6%`, Q2 `21.9%`, Q3 `23.0%`, Q4 `19.2%`, Q5 `27.2%`
- Positive modeled loss share: Q1 `16.2%`, Q2 `30.5%`, Q3 `33.5%`, Q4 `46.7%`, Q5 `41.0%`
- Mean annual loss probability: Q1 `0.9%`, Q2 `1.3%`, Q3 `2.0%`, Q4 `4.8%`, Q5 `3.3%`
- Median finite return period among positive-loss municipalities: Q1 `156`, Q2 `61`, Q3 `68`, Q4 `40`, Q5 `28`

Interpretation is now correctly cautious:

- The thesis does not claim that vulnerable municipalities are generally more exposed.
- It does not claim that rich municipalities are simply better protected.
- It does not claim policy intent.
- It does claim that the modeled loss/protection layer is more socially differentiated than municipal area exposure alone.

## Figure Check

The figure paths inserted in the draft point to existing PNG files. The Results figures are introduced in the surrounding text and are not floating randomly.

Remaining layout decision:

- Figure 2.1 and Figure 5.1 currently use the same study-area map.
- This is not a data error, but it is a final layout choice. Either keep both for readability, or keep Figure 2.1 and have Chapter 5 refer back to it.

## Major Remaining Submission Risks

### 1. Title alignment under fixed-title constraint

The current title is:

`Spatial distribution of flood protection infrastructure across socially relevant factors in the Elbe river basin.`

Fixed constraint:

- The title is the title of the thesis and should not be changed casually.
- The task is therefore not to replace the title, but to make the whole thesis clearly derive from it.

Hard reviewer concern under this constraint:

- The empirical thesis is no longer mainly about a physical inventory of flood protection infrastructure.
- It is about modeled flood exposure, socio-economic vulnerability, modeled loss occurrence, annual loss probability, and finite protection return periods.
- It is empirically corridor-based, not full-basin-wide in the analytical sample.

Resolution:

- Keep the title.
- Make the operationalization explicit in the Introduction and Methods.
- State clearly that the thesis does not inventory individual dikes, walls, or local protection assets.
- Interpret `flood protection infrastructure` through modeled protection-related outcomes: positive modeled loss occurrence, modeled annual loss probability, finite protection return periods, and no-event municipalities.
- Interpret `socially relevant factors` through the INKAR-based vulnerability dimensions.
- Interpret `Elbe river basin` as the broader hydrological and governance frame, while the empirical analysis focuses on the RP500 flood-relevant corridor.

### 2. Missing bibliography

The draft contains a `LITERATURE` entry in the table of contents, but no actual bibliography section at the end.

This is a real blocker. Before submission, the bibliography must include all in-text citations and all data sources. At minimum, it must include:

- White (1945)
- Jongman et al. (2012)
- Rojas et al. (2013)
- Thieken et al. (2016)
- Grabs (2016)
- Jüpner (2018)
- Barendrecht et al. (2018)
- Lindenschmidt et al. (2006)
- Tovar Reaños (2021)
- Odersky and Löffler (2024)
- Qiang (2019)
- Poussard et al. (2021)
- Kephart et al. (2025)
- Fekete (2009, 2010)
- Rufat et al. (2015, 2019)
- Schmidtlein et al. (2008)
- Tate (2013)
- de Moel et al. (2015)
- de Goër de Herve (2022)
- Thaler and Hartmann (2016)
- Thaler (2021)
- Johnson et al. (2007)
- Liao et al. (2019)
- O'Hare and White (2018)
- Moulds et al. (2021)
- Fielding (2012)
- Shao et al. (2022)
- Sivapalan et al. (2012)
- Tobin (1995)
- Ferdous et al. (2020)
- Serra-Llobet et al. (2022)
- Fu et al. (2023)
- Kreibich et al. (2017)
- Osberghaus (2015)
- Bubeck et al. (2013)
- Dillenardt and Thieken (2025)
- Kuhlicke et al. (2020)
- Surminski and Thieken (2017)
- Sairam et al. (2021)
- Alfieri et al. (2014)
- Dottori et al. (2022)
- Baugh et al. (2024)
- BKG/VG250
- BBSR/INKAR, including edition and territorial status
- Provider clarification as personal communication or appendix documentation

### 3. Table of contents

The current table of contents is still a working table. It contains fixed page numbers and page-count notes. For final submission, use an automatically generated Word/LibreOffice table of contents, or remove manual page numbers from the Markdown draft.

### 4. Personal communication and provider data

The provider clarification is correctly integrated conceptually, but the final thesis needs a clean citation convention:

- Either cite as `N. Sairam and T. N. M. Phan, personal communication, July 2026`
- Or document the exchange in an appendix and cite the appendix in Chapter 4 and Chapter 6.

### 5. Data-source citation details

Chapter 4 now names the core data sources. The final bibliography still needs precise entries for:

- BKG VG250
- BBSR INKAR 07/2025, territorial status 31.12.2023
- JRC flood hazard dataset: Baugh et al. (2024), DOI `10.2905/1D128B6C-A4EE-4858-9E34-6210707F3C81`

## Logic And Red Thread

The red thread is now coherent:

1. Flood risk is socially and politically produced, not just hydrological.
2. The Elbe is a flood-experienced, protected, socially uneven river corridor.
3. The thesis asks a distributive question at municipality level.
4. The method builds three layers: exposure, vulnerability, protection/loss.
5. The Results show that these layers do not tell the same story.
6. The Discussion argues that this difference is the main flood-justice insight.
7. The Conclusion closes with a layered, cautious justice claim.

The most important strategic strength is that the thesis does not force the expected environmental-justice result. It accepts the weak exposure-vulnerability association and then shows why protection/loss outcomes tell a more socially differentiated story.

## Justice And Critical Perspective

The draft now contains a more critical justice reading, but it remains data-bound. It connects protection to land values, administrative capacity, visibility, and investment priorities without claiming to have measured capital flows, policy intent, or historical protection investment.

This is the right level for the available data. A stronger capitalism critique would require additional evidence on land markets, fiscal capacity, protection investment histories, insurance, or political decision-making. Without that, the thesis should stay with the more defensible claim:

Flood protection and modeled loss visibility are not neutral technical categories. They are embedded in socio-economic, territorial, and institutional priorities.

## Next Best Fixes

1. Build the final bibliography from the cited literature and data sources.
2. Keep the fixed title and make sure the Introduction, Methods, and Discussion consistently explain how it is operationalized.
3. Remove or regenerate the manual table of contents.
4. Decide whether to keep both Figure 2.1 and Figure 5.1.
5. Do one final line edit after the bibliography is inserted, because reference formatting often exposes remaining style inconsistencies.
