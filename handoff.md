# Dexterity × DTI Study — Handoff

_Last updated: 2026-06-09_

## Goal

Test the a priori hypothesis that **higher manual dexterity is associated with higher
white-matter quality in the occipital gyri** (inferior / middle / superior occipital,
L+R), using regional FA and MD from the ABC cohort. Grip strength is included as a
**contrast / specificity control**, not a co-headline — the occipital story is a
**dexterity** story.

Framing for any eventual manuscript: a *confirmatory a priori test* of an occipital
visuomotor-integration hypothesis, situated within (not claimed to be unique to) a
broader sensorimotor–parietal–occipital network.

## Current state

Exploratory pass is **done and stands up well**. Within the pre-specified 6-ROI
occipital family (× FA/MD × dominant/non-dominant dexterity = 24 tests):

- **21 / 24 significant at uncorrected p<0.05, every one in the coherent direction**
  (dexterity↑ → FA↑, MD↓).
- Under **Bonferroni within the 6-ROI family** (α≈.008, the lab's standard primary
  correction), **3–4 ROIs survive per metric** for dominant-hand dexterity.
- Strongest: inferior occipital R (FA p=.0001), middle occipital R (FA p=.0005),
  middle occipital L (FA p=.001).
- Consistent weak spot: **superior occipital RIGHT** (the 3 non-sig cells are here) —
  report honestly, don't hide.
- **Grip strength does NOT track occipital** (its few MD hits are limbic/subcortical:
  hippocampus, parahippocampal, temporal pole, thalamus, caudate). Good — it works as
  a specificity contrast.

Covariate model: OLS, controlling **age, sex, BMI, race (one-hot race___1–6)**.
n ≈ 285–287 per model after merge. **No FDR/Bonferroni applied in the code yet** — all
correction is currently done by eye downstream. Exploratory by design.

### Mechanistic disambiguation (multi-metric) — DONE, and it's the novelty hook

Per the lit search, a plain FA correlation is NOT novel (occipital/visuomotor WM ↔
dexterity is published in aging/MS cohorts). The differentiator is **what drives the
FA effect**, using b-robust metrics most prior work lacks. Result on the occipital
family × dominant dexterity (`occipital_microstructure.py`):

- **FA↑ (6/6), MD↓ (5/6), RD↓ (5/6)**, with **AD mostly flat (1/6)** → the FA effect is
  **radial-diffusivity-driven → a myelin signature**, not axonal.
- **FSUM↑ (6/6)** → more total fiber volume = genuine fiber content, not just geometry.
- **MO flat (≤1/6)** and **F2↑ (6/6, i.e. *more* secondary fiber, opposite to a naive
  "less crossing" guess)** → the FA increase is **NOT a crossing-fiber/tensor-shape
  artifact**; better dexterity goes with richer, better-myelinated architecture despite
  more crossing. Report F2 transparently — it's unexpected but *supports* the anti-
  artifact argument.
- **NODDI (experimental, n~214):** NDI↑, ODI↓, FWF↓ — all 6/6 in the "better quality"
  direction (2–3/6 individually significant given smaller n). **FWF↓ directly rebuts the
  gray-matter partial-volume/CSF critique.** Keep these clearly labeled experimental
  (b2000 not b3000); use as converging support, not primary tests.

Bottom line for the manuscript: frame the paper around **microstructural quality
(myelination + fiber/neurite content), not a partial-volume or crossing-fiber artifact**
— that's the defensible new contribution.

### Supplementary sex analysis — DONE (`occipital_sex_analysis.py`)

Cohort is ~75% female (males n~67, females n~219 after merge — male stratum is the
power limit). On the occipital family × dominant dexterity:

- Association is present and **same-direction in BOTH sexes**.
- **Numerically larger in males** across every metric (FA r≈+.20–.33 vs +.11–.25; FSUM
  r≈+.26–.36 vs +.11–.23) despite the much smaller male n.
- **But the formal sex×dexterity interaction is essentially null** (sig in only 2/24
  ROIs, uncorrected; would not survive correction). So: **treat sex as a covariate, not
  a moderator** — the honest claim is generalizability across sex, with a male-stronger
  *pattern* flagged as a hypothesis for a sex-balanced replication. Do NOT overclaim a
  sex difference.

## Files

### In this folder (study-specific, moved here)
- `grip_dexterity_fa_md.py` — main analysis. Reads raw CSVs from project ROOT, writes
  `grip_dexterity_fa_md_results.csv` **into this folder**. Predictors = NIH Toolbox
  *uncorrected standard* grip + 9-hole-peg scores (dominant + non-dominant).
- `grip_dexterity_stratify.py` — reads the results CSV from this folder; prints
  (1) FA∩MD overlap, (2) all sig ROIs, (3) ROIs common to all 4 predictors.
- `grip_dexterity_fa_md_results.csv` — full per-test output (1512 rows: metric ×
  predictor × 189 ROIs, with beta, partial r, p, n).
- `occipital_microstructure.py` / `occipital_microstructure_results.csv` — the
  mechanistic disambiguation (FA/MD/AD/RD/MO/F2/FSUM + experimental NDI/ODI/FWF) on the
  6-ROI occipital family × dexterity. **This is the novelty engine.**
- `occipital_sex_analysis.py` / `occipital_sex_results.csv` — sex×dexterity interaction +
  within-sex stratified estimates on the occipital family.
- `Dexterity_DTI_abstract.docx` — short structured abstract (predates the mechanistic +
  sex results; **needs a refresh** to fold them in — see next steps).
- `handoff.md` — this file.

### Reference CSVs — LEFT IN PROJECT ROOT (shared across studies, do NOT move/edit)
- `gripstrength-dexterity.csv` — behavioral predictors (subjects with data use ABC####
  IDs; integer-ID rows are empty placeholders, filtered in code).
- `DTIstuff/baseline/FA_jhu.csv`, `DTIstuff/baseline/MD_jhu.csv` — 189 JHU ROIs, 293 subj.
- `ABCStudyBehavioralDa-05282025ABCDemograph_DATA_2026-05-24_1638.csv` — age
  (`age_dayofday1testing`), `sex`, `race___1..6`.
- `BMI.csv` — `bmi_calculated`.

## Files in flight / what changed this session

- Created folder `DexterityDTI study/`; moved the analysis files in; repointed script
  path constants (`base` = ROOT for raw data, `here` = this folder for derived files).
  All scripts verified to run from the new location.
- Ran a literature check (3 web searches) → established the plain FA correlation is NOT
  novel; pivoted the study's contribution to the multi-metric mechanistic disambiguation.
- Added & ran `occipital_microstructure.py` (AD/RD/MO/F2/FSUM + NODDI) and
  `occipital_sex_analysis.py`. Both wrote results CSVs. See Current state for findings.
- `Dexterity_DTI_abstract.docx` is now **stale** — written before these two analyses.

## Failed attempts / gotchas (so the next person doesn't repeat them)

- **ID mismatch trap:** `gripstrength-dexterity.csv` contains placeholder rows with
  plain integer subject_ids (1,2,3…, all NIH-toolbox-incomplete). Real data rows use
  `ABC####`. Code filters to `startswith("ABC")`. Don't merge on the raw column blindly.
- **Don't double-correct age:** we deliberately use the **uncorrected_standard** NIH
  Toolbox columns, NOT the `_agecorr` ones, because age is already a regression
  covariate. Switching to agecorr columns would correct for age twice.
- **race dummies with no variance:** within a merged subset some `race___X` columns are
  all-zero; the code drops zero-variance covariates per model to avoid singular design.
- Trivial `%`-format bug in a print line (fixed) — `0.05*len` was hitting `%o`; cosmetic.
- **`grip_strength_dominant_score` vs `_uncorrected_standard_dominant`:** we chose the
  standard score. The raw `_score` (kg force / seconds) is an untested alternative — see
  next steps.

## Next things to try (in priority order)

1. **Visual-acuity covariate (highest value).** The grip CSV has `visual_acuity_*`
   columns. Re-run the occipital family adding visual acuity as a covariate. If occipital
   survives, it's visuomotor white-matter integrity, not "can they see the pegs" — this
   preempts the #1 reviewer objection and is arguably the single most important analysis.
2. **Specificity / global-brain control (still the make-or-break analysis).** Dexterity
   hit 40+ ROIs, so we CANNOT yet claim occipital is special. Add whole-brain mean FA/MD
   as a covariate (and/or a random-ROI permutation null) and show the occipital effect
   survives / exceeds baseline. Until this is run, the occipital-*specific* claim is
   unproven. **Do this before any more writing.**
3. ~~NDI / microstructure layer~~ — **DONE** (`occipital_microstructure.py`). See Current
   state. This is now the paper's core, not a future step.
4. **Refresh the abstract docx** to lead with the mechanistic (RD/FSUM/NDI) disambiguation
   and add the sex supplementary; the current `Dexterity_DTI_abstract.docx` predates both.
4. **Formalize the a priori family + correction.** Lock the 6-ROI (or per-lab 12-ROI)
   family, apply Bonferroni within-family as primary in code (currently by-eye). Decide
   whether superior-occipital-R stays in the family or is reported as the predicted miss.
5. **Sensitivity check:** re-run with raw `_score` columns instead of uncorrected
   standard scores; confirm the occipital pattern is robust to that choice.
6. **Handedness:** consider modeling dominant-hand dexterity as primary (it carries the
   strongest signal) and treating non-dominant as confirmation.

## How to reproduce

```powershell
cd "c:\Users\bandd\Downloads\52426_DepressionWork\DexterityDTI study"
python grip_dexterity_fa_md.py      # regenerates results CSV
python grip_dexterity_stratify.py   # prints stratified ROI breakdown
```
