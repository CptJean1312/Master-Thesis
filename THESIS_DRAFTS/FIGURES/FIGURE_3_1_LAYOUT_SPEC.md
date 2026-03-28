# Figure 3.1 Layout Specification

## Purpose of the figure

This figure should be a **conceptual model**, not a methods workflow.

It should show:

- the core analytical components of the thesis
- how they relate conceptually
- how the distributive justice perspective emerges from these relationships

It should **not** show:

- data source names such as `VG250`, `INKAR`, `EFAS`
- GIS preprocessing steps
- reprojection, clipping, raster extraction, PCA steps, or quality checks
- any wording such as `later`, `staged`, `to be done`

The figure must stay one level above the methods chapter.

---

## Overall composition

Recommended layout: **five-box conceptual model**

- one broad contextual banner at the top
- three core boxes in the middle/left-to-right analytical chain
- one supporting box below
- one final assessment box on the right

This gives the figure a stable reading order:

1. study setting
2. hazard
3. exposure
4. vulnerability / protection
5. distributive assessment

---

## Box placement

### 1. Top banner

**Position**
- top center, spanning most of the figure width

**Function**
- anchors the figure spatially without becoming methodological

**Suggested heading**
- `Study setting`

**Suggested text**
- `Municipality-level analysis in the Elbe flood corridor`

This banner should be visually quieter than the main analytical boxes.

---

### 2. Left box

**Position**
- left side, upper-middle

**Heading**
- `Flood Hazard`

**Function**
- starting point of the conceptual logic

**Suggested body text**
- `Physical flood potential across return periods`
- `Spatial flood intensity and extent`

**Optional small line**
- `Hazard exists independently of social conditions`

This box should communicate that hazard is the physical basis of the analysis.

---

### 3. Center box

**Position**
- center of the figure
- largest and visually most important box

**Heading**
- `Flood Exposure`

**Function**
- this is the central analytical hinge of the figure
- it connects the physical side to the social side

**Suggested body text**
- `Realized municipal flood exposure`
- `Spatial flood outcomes at municipality level`
- `Can vary across return periods and locations`

**Optional small line**
- `Central empirical expression of flood-related outcomes`

This should be the box the eye goes to after hazard.

---

### 4. Lower-left box

**Position**
- bottom left, below hazard

**Heading**
- `Socio-Economic Vulnerability`

**Function**
- introduces the social conditioning dimension

**Suggested body text**
- `Multidimensional social and demographic conditions`
- `Shapes sensitivity and adaptive capacity`

**Optional small line**
- `Social differentiation of flood-related disadvantage`

This box should not be as visually dominant as exposure, but it should be clearly legible.

---

### 5. Upper-middle or upper-right box

**Position**
- above exposure, slightly right of center

**Heading**
- `Flood Protection`

**Function**
- show that protection modifies the translation from hazard into realized exposure
- and also matters for residual risk and uneven benefit distribution

**Suggested body text**
- `Modifies how hazard becomes exposure`

**Optional small line**
- `Reduces frequent impacts, but does not remove all risk`

This box should not look like a data input box.
It should look conceptually connected to hazard and exposure.

---

### 6. Right box

**Position**
- right side, middle to lower-middle

**Heading**
- `Distributive Assessment`

**Function**
- endpoint of the conceptual logic
- shows what the thesis actually evaluates

**Suggested body text**
- `Compare the geography of exposure`
- `with the geography of vulnerability`
- `Assess whether flood-related burdens`
- `and protection-related benefits are unevenly distributed`

This box should read as the final interpretive step.

---

## Arrows

Use only a few arrows.
Too many arrows make the figure look like a workflow chart.

### Arrow 1

**From**
- `Flood Hazard`

**To**
- `Flood Exposure`

**Line type**
- solid

**Label**
- `hazard translated into exposure`

**Meaning**
- physical flood potential becomes realized municipal flood exposure

This is the main arrow of the figure.

---

### Arrow 2

**From**
- `Flood Protection`

**To**
- `Flood Exposure`

**Line type**
- dashed

**Label**
- `modifies realized exposure`

**Meaning**
- protection changes how much of the physical hazard becomes actual exposure

Dashed is useful here because protection is conceptually important, but not identical to hazard itself.

---

### Arrow 3

**From**
- `Flood Exposure`

**To**
- `Distributive Assessment`

**Line type**
- solid

**Label**
- `observed flood-related outcome`

**Meaning**
- exposure is one of the main outcomes the thesis evaluates

---

### Arrow 4

**From**
- `Socio-Economic Vulnerability`

**To**
- `Distributive Assessment`

**Line type**
- solid

**Label**
- `social interpretation of flood outcomes`

**Meaning**
- exposure becomes justice-relevant when read together with vulnerability

This arrow is crucial because it visually expresses the core argument of the thesis.

---

### Arrow 5

**From**
- `Flood Protection`

**To**
- `Distributive Assessment`

**Line type**
- dashed

**Label**
- `shapes residual risk and benefit distribution`

**Meaning**
- protection does not only affect exposure directly
- it also matters for how benefits and remaining burdens are distributed

This arrow should be thinner or visually lighter than the core exposure arrows.

---

## Caption logic

The caption should stay conceptual and not become methodological.

### Recommended caption

`Figure 3.1. Conceptual model of the study. Flood hazard creates the physical potential for inundation, while flood exposure captures how this hazard is realized at municipality level. Socio-economic vulnerability represents the social conditions that shape sensitivity and adaptive capacity, and flood protection modifies how physical hazard is translated into realized exposure and residual risk. The study evaluates these dimensions through a distributive perspective that compares the geography of exposure and protection-related benefits with the geography of socio-economic vulnerability.`

---

## What should NOT be written inside the figure

Avoid:

- `RP100`
- `RP500`
- `INKAR`
- `VG250`
- `PCA`
- `flood share`
- `municipality boundaries`
- `return-period curve metrics`
- `effective protection proxy`

Those belong in Chapter 4 or in the surrounding Chapter 3 text, not in the conceptual model.

If these terms appear in the figure, it will start looking like a methods diagram instead of a conceptual one.

---

## Visual guidance

### General style

- very light fills
- thin to medium outlines
- serif typography is fine
- minimal color contrast
- generous spacing

### Hierarchy

- `Flood Exposure` should be most visually central
- `Distributive Assessment` should be clearly the endpoint
- `Flood Protection` should be visibly connected, but not overpower the core chain

### Preferred emphasis

Visual hierarchy should suggest:

`Hazard -> Exposure -> Distributive Assessment`

with:

- `Vulnerability` giving social meaning
- `Protection` modifying the translation from hazard to exposure

---

## Short version for manual rebuilding in Inkscape

If rebuilding manually:

1. Put a horizontal banner at the top:
   - `Study setting: Municipality-level analysis in the Elbe flood corridor`

2. Put `Flood Hazard` on the left.

3. Put `Flood Exposure` in the middle, slightly larger than all other boxes.

4. Put `Socio-Economic Vulnerability` below left.

5. Put `Flood Protection` above exposure.

6. Put `Distributive Assessment` on the right.

7. Draw these arrows:
   - `Flood Hazard -> Flood Exposure`
   - `Flood Protection -> Flood Exposure`
   - `Flood Exposure -> Distributive Assessment`
   - `Socio-Economic Vulnerability -> Distributive Assessment`
   - `Flood Protection -> Distributive Assessment`

8. Use these labels:
   - `hazard translated into exposure`
   - `modifies realized exposure`
   - `observed flood-related outcome`
   - `social interpretation of flood outcomes`
   - `shapes residual risk and benefit distribution`

That is the cleanest conceptual version.
