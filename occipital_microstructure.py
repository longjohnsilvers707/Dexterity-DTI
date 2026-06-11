"""
Mechanistic disambiguation of the occipital dexterity effect.

Question: when higher dexterity tracks higher FA / lower MD in the occipital gyri,
is that (a) genuine white-matter quality (-> RD down, FSUM up, NDI up) or
(b) merely reduced crossing-fiber / dispersion geometry, i.e. a partial-volume /
tensor-shape artifact (-> F2 down, ODI down, MO up, with FSUM/NDI flat)?

Primary (b-robust, 295 subj): FA, MD, AD, RD, MO, F2, FSUM
Supplementary (experimental NODDI, b2000 not b3000, 217 subj): NDI, ODI, FWF

Same covariate model as the main study: age, sex, BMI, race (one-hot).
A priori occipital family = inferior/middle/superior occipital gyri, L+R.
"""
import pandas as pd, numpy as np, statsmodels.api as sm, re, warnings
warnings.filterwarnings("ignore")

base = r"c:\Users\bandd\Downloads\52426_DepressionWork"
here = r"c:\Users\bandd\Downloads\52426_DepressionWork\DexterityDTI study"

OCC = ["IOG_L", "IOG_R", "MOG_L", "MOG_R", "SOG_L", "SOG_R"]   # short-name family
# long-name equivalents in the FA/MD JHU files (note caps quirk for MOG_R)
OCC_LONG = {
    "IOG_L": "inferior occipital gyrus left", "IOG_R": "inferior occipital gyrus right",
    "MOG_L": "middle occipital gyrus left",   "MOG_R": "MIDDLE OCCIPITAL GYRUS right",
    "SOG_L": "superior occipital gyrus left", "SOG_R": "superior occipital gyrus right",
}

# directional hypothesis: sign of beta expected if "higher dexterity = better quality"
EXPECT = {"FA": +1, "MD": -1, "AD": None, "RD": -1, "MO": +1,
          "F2": -1, "FSUM": +1, "NDI": +1, "ODI": -1, "FWF": -1}
TIER = {m: "primary" for m in ["FA", "MD", "AD", "RD", "MO", "F2", "FSUM"]}
TIER.update({m: "experimental" for m in ["NDI", "ODI", "FWF"]})

# ---- predictors ----
grip = pd.read_csv(f"{base}/gripstrength-dexterity.csv", dtype={"subject_id": str})
preds = {"dexterity_dominant":    "nine_hole_peg_dexterity_uncorrected_standard_dominant",
         "dexterity_nondominant": "nine_hole_peg_dexterity_uncorrected_standard_nondominant"}
grip = grip[["subject_id"] + list(preds.values())].rename(columns={v: k for k, v in preds.items()})
grip = grip[grip["subject_id"].str.startswith("ABC", na=False)].rename(columns={"subject_id": "PatientID"})

# ---- covariates ----
demo = pd.read_csv(f"{base}/ABCStudyBehavioralDa-05282025ABCDemograph_DATA_2026-05-24_1638.csv",
                   dtype={"subject_id": str})
race_cols = [c for c in demo.columns if re.fullmatch(r"race___\d", c)]
cov = demo[["subject_id", "age_dayofday1testing", "sex"] + race_cols].rename(
    columns={"subject_id": "PatientID", "age_dayofday1testing": "age"})
bmi = pd.read_csv(f"{base}/BMI.csv", dtype={"subject_id": str})[["subject_id", "bmi_calculated"]].rename(
    columns={"subject_id": "PatientID", "bmi_calculated": "BMI"})
covars = ["age", "sex", "BMI"] + race_cols

# ---- metric file sources ----
short_files = {"AD": f"{base}/AD_JHU_ROI.csv", "RD": f"{base}/RD_JHU_ROI.csv",
               "MO": f"{base}/MO_JHU_ROI.csv", "F2": f"{base}/F2_JHU_ROI.csv",
               "FSUM": f"{base}/FSUM_JHU_ROI.csv",
               "NDI": f"{base}/NDI_JHU_ROI_experimental.csv",
               "ODI": f"{base}/ODI_JHU_ROI_experimental.csv",
               "FWF": f"{base}/FWF_JHU_ROI_experimental.csv"}
long_files = {"FA": f"{base}/DTIstuff/baseline/FA_jhu.csv",
              "MD": f"{base}/DTIstuff/baseline/MD_jhu.csv"}

def load_metric(metric):
    if metric in long_files:
        df = pd.read_csv(long_files[metric])
        keep = {OCC_LONG[r]: r for r in OCC if OCC_LONG[r] in df.columns}
        df = df[["PatientID"] + list(keep)].rename(columns=keep)
    else:
        df = pd.read_csv(short_files[metric])
        df = df[["PatientID"] + [r for r in OCC if r in df.columns]]
    return df

def fit(df, pred, roi):
    sub = df[[pred, roi] + covars].apply(pd.to_numeric, errors="coerce").dropna()
    usecov = [c for c in covars if sub[c].nunique() > 1]
    if len(sub) < len(usecov) + 5:
        return None
    m = sm.OLS(sub[roi], sm.add_constant(sub[[pred] + usecov])).fit()
    t = m.tvalues[pred]
    return dict(n=int(len(sub)), beta=m.params[pred], p=m.pvalues[pred],
                r_partial=np.sign(m.params[pred]) * np.sqrt(t**2 / (t**2 + m.df_resid)))

rows = []
for metric in ["FA", "MD", "AD", "RD", "MO", "F2", "FSUM", "NDI", "ODI", "FWF"]:
    mdf = load_metric(metric)
    df = grip.merge(cov, on="PatientID").merge(bmi, on="PatientID").merge(mdf, on="PatientID")
    for c in covars:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for pred in preds:
        for roi in [r for r in OCC if r in mdf.columns]:
            res = fit(df, pred, roi)
            if res:
                rows.append(dict(metric=metric, tier=TIER[metric], predictor=pred, roi=roi, **res))

R = pd.DataFrame(rows)
R.to_csv(f"{here}/occipital_microstructure_results.csv", index=False)

# ===================== report =====================
def block(title): print("\n" + "=" * 76 + f"\n{title}\n" + "=" * 76)

block("OCCIPITAL FAMILY (6 ROIs) x DEXTERITY  |  covariates: age, sex, BMI, race")
print("dir-OK = beta sign matches the 'higher dexterity = better quality' hypothesis")
for metric in ["FA", "MD", "AD", "RD", "MO", "F2", "FSUM", "NDI", "ODI", "FWF"]:
    for pred in preds:
        sub = R[(R.metric == metric) & (R.predictor == pred)]
        if sub.empty:
            continue
        nsig = int((sub.p < 0.05).sum())
        exp = EXPECT[metric]
        ndir = "n/a" if exp is None else int((np.sign(sub.beta) == exp).sum())
        nsig_dir = "n/a" if exp is None else int(((sub.p < 0.05) & (np.sign(sub.beta) == exp)).sum())
        n_typ = int(sub.n.median())
        print(f"  {metric:>4} [{TIER[metric][:4]}] {pred:<22} "
              f"n={n_typ:<4} sig={nsig}/6  same-direction={ndir}/6  sig&expected={nsig_dir}/6")

block("PER-ROI DETAIL (dominant-hand dexterity; the stronger signal)")
for metric in ["FA", "MD", "AD", "RD", "MO", "F2", "FSUM", "NDI", "ODI", "FWF"]:
    sub = R[(R.metric == metric) & (R.predictor == "dexterity_dominant")].set_index("roi")
    if sub.empty:
        continue
    cells = []
    for roi in OCC:
        if roi in sub.index:
            x = sub.loc[roi]
            star = "*" if x.p < 0.05 else " "
            cells.append(f"{roi}:{x.r_partial:+.2f}{star}")
        else:
            cells.append(f"{roi}:--")
    print(f"  {metric:>4} [{TIER[metric][:4]}]  " + "  ".join(cells))
print("\n  (* p<0.05; value = partial r)")
