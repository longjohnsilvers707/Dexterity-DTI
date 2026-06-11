"""
Supplementary sex analysis on the a priori occipital family.

Two questions:
  (1) MODERATION  - does sex moderate the dexterity-occipital association?
                    (sex x dexterity interaction term in the full-sample model)
  (2) STRATIFIED  - is the association present within each sex separately?
                    (male-only and female-only partial r, each adjusting age, BMI, race)

Outcomes: FA, MD (primary) + RD, FSUM (mechanistic leads). Predictor: dominant-hand
dexterity (the stronger signal). Cohort is ~75% female, so the male stratum is the
power-limiting one - interpret male estimates with that caveat.
"""
import pandas as pd, numpy as np, statsmodels.api as sm, re, warnings
warnings.filterwarnings("ignore")

base = r"c:\Users\bandd\Downloads\52426_DepressionWork"
here = r"c:\Users\bandd\Downloads\52426_DepressionWork\DexterityDTI study"
OCC = ["IOG_L", "IOG_R", "MOG_L", "MOG_R", "SOG_L", "SOG_R"]
OCC_LONG = {"IOG_L": "inferior occipital gyrus left", "IOG_R": "inferior occipital gyrus right",
            "MOG_L": "middle occipital gyrus left", "MOG_R": "MIDDLE OCCIPITAL GYRUS right",
            "SOG_L": "superior occipital gyrus left", "SOG_R": "superior occipital gyrus right"}
PRED = "dexterity_dominant"

grip = pd.read_csv(f"{base}/gripstrength-dexterity.csv", dtype={"subject_id": str})
grip = grip[["subject_id", "nine_hole_peg_dexterity_uncorrected_standard_dominant"]].rename(
    columns={"subject_id": "PatientID", "nine_hole_peg_dexterity_uncorrected_standard_dominant": PRED})
grip = grip[grip["PatientID"].str.startswith("ABC", na=False)]

demo = pd.read_csv(f"{base}/ABCStudyBehavioralDa-05282025ABCDemograph_DATA_2026-05-24_1638.csv",
                   dtype={"subject_id": str})
race_cols = [c for c in demo.columns if re.fullmatch(r"race___\d", c)]
cov = demo[["subject_id", "age_dayofday1testing", "sex"] + race_cols].rename(
    columns={"subject_id": "PatientID", "age_dayofday1testing": "age"})
bmi = pd.read_csv(f"{base}/BMI.csv", dtype={"subject_id": str})[["subject_id", "bmi_calculated"]].rename(
    columns={"subject_id": "PatientID", "bmi_calculated": "BMI"})

def load(metric):
    if metric in ("FA", "MD"):
        df = pd.read_csv(f"{base}/DTIstuff/baseline/{metric}_jhu.csv")
        keep = {OCC_LONG[r]: r for r in OCC if OCC_LONG[r] in df.columns}
        return df[["PatientID"] + list(keep)].rename(columns=keep)
    df = pd.read_csv(f"{base}/{metric}_JHU_ROI.csv")
    return df[["PatientID"] + [r for r in OCC if r in df.columns]]

def partial_r(m, term):
    t = m.tvalues[term]
    return np.sign(m.params[term]) * np.sqrt(t**2 / (t**2 + m.df_resid))

rows = []
for metric in ["FA", "MD", "RD", "FSUM"]:
    mdf = load(metric)
    df = grip.merge(cov, on="PatientID").merge(bmi, on="PatientID").merge(mdf, on="PatientID")
    for c in ["age", "sex", "BMI"] + race_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df[df["sex"].isin([1, 2])].copy()
    df["female"] = (df["sex"] == 2).astype(int)
    for roi in [r for r in OCC if r in mdf.columns]:
        base_cols = [PRED, roi, "age", "BMI", "female"] + race_cols
        d = df[base_cols].apply(pd.to_numeric, errors="coerce").dropna()
        d["dex_c"] = d[PRED] - d[PRED].mean()
        usecov = [c for c in ["age", "BMI"] + race_cols if d[c].nunique() > 1]
        # (1) interaction model
        X = sm.add_constant(d[["dex_c", "female"] + usecov].assign(dexXfem=d["dex_c"] * d["female"]))
        mi = sm.OLS(d[roi], X).fit()
        # (2) stratified
        strat = {}
        for lab, sub in (("M", d[d.female == 0]), ("F", d[d.female == 1])):
            uc = [c for c in ["age", "BMI"] + race_cols if sub[c].nunique() > 1]
            if len(sub) < len(uc) + 5:
                strat[lab] = (np.nan, np.nan, len(sub)); continue
            ms = sm.OLS(sub[roi], sm.add_constant(sub[[PRED] + uc])).fit()
            strat[lab] = (partial_r(ms, PRED), ms.pvalues[PRED], int(len(sub)))
        rows.append(dict(metric=metric, roi=roi,
                         interaction_p=mi.pvalues["dexXfem"],
                         rM=strat["M"][0], pM=strat["M"][1], nM=strat["M"][2],
                         rF=strat["F"][0], pF=strat["F"][1], nF=strat["F"][2]))

R = pd.DataFrame(rows)
R.to_csv(f"{here}/occipital_sex_results.csv", index=False)

print("=" * 84)
print("SEX ANALYSIS - occipital family x dominant dexterity (covariates: age, BMI, race)")
print("=" * 84)
print(f"Cohort after dropping non-1/2 sex & NaN: males n~{int(R.nM.median())}, females n~{int(R.nF.median())}")
print("int_p = sex x dexterity interaction p (full sample) | M/F blocks = within-sex partial r\n")
hdr = f"{'metric':>5} {'roi':>6} | {'int_p':>7} | {'r_Male':>7} {'p_M':>6} | {'r_Female':>8} {'p_F':>6}"
print(hdr); print("-" * len(hdr))
for _, x in R.iterrows():
    sM = "*" if x.pM < 0.05 else " "
    sF = "*" if x.pF < 0.05 else " "
    sI = "*" if x.interaction_p < 0.05 else " "
    print(f"{x.metric:>5} {x.roi:>6} | {x.interaction_p:>6.3f}{sI} | {x.rM:>+7.2f} {x.pM:>6.3f}{sM} | {x.rF:>+8.2f} {x.pF:>6.3f}{sF}")

print("\nSummary:")
for metric in ["FA", "MD", "RD", "FSUM"]:
    s = R[R.metric == metric]
    print(f"  {metric:>4}: interaction sig in {int((s.interaction_p<0.05).sum())}/6 ROIs | "
          f"within-sex sig: M {int((s.pM<0.05).sum())}/6, F {int((s.pF<0.05).sum())}/6")
