import pandas as pd, numpy as np, statsmodels.api as sm, re, warnings
warnings.filterwarnings("ignore")

# reference/raw CSVs stay in the project root; this study's outputs live alongside this script
base = r"c:\Users\bandd\Downloads\52426_DepressionWork"
here = r"c:\Users\bandd\Downloads\52426_DepressionWork\DexterityDTI study"

# --- behavioral predictors (grip strength + 9-hole peg dexterity) ---
grip = pd.read_csv(f"{base}/gripstrength-dexterity.csv", dtype={"subject_id": str})
predictors = {
    "grip_dominant":        "grip_strength_uncorrected_standard_dominant",
    "grip_nondominant":     "grip_strength_uncorrected_standard_nondominant",
    "dexterity_dominant":   "nine_hole_peg_dexterity_uncorrected_standard_dominant",
    "dexterity_nondominant":"nine_hole_peg_dexterity_uncorrected_standard_nondominant",
}
grip = grip[["subject_id"] + list(predictors.values())].rename(columns={v: k for k, v in predictors.items()})
grip = grip[grip["subject_id"].str.startswith("ABC", na=False)]
grip = grip.rename(columns={"subject_id": "PatientID"})

# --- covariates: age, race, sex, BMI ---
demo = pd.read_csv(f"{base}/ABCStudyBehavioralDa-05282025ABCDemograph_DATA_2026-05-24_1638.csv",
                   dtype={"subject_id": str})
race_cols = [c for c in demo.columns if re.fullmatch(r"race___\d", c)]
cov = demo[["subject_id", "age_dayofday1testing", "sex"] + race_cols].rename(
    columns={"subject_id": "PatientID", "age_dayofday1testing": "age"})

bmi = pd.read_csv(f"{base}/BMI.csv", dtype={"subject_id": str})[["subject_id", "bmi_calculated"]]
bmi = bmi.rename(columns={"subject_id": "PatientID", "bmi_calculated": "BMI"})

# --- imaging ---
fa = pd.read_csv(f"{base}/DTIstuff/baseline/FA_jhu.csv")
md = pd.read_csv(f"{base}/DTIstuff/baseline/MD_jhu.csv")
roi_cols = [c for c in fa.columns if c != "PatientID"]

covars = ["age", "sex", "BMI"] + race_cols

def run_metric(img, metric_name):
    df = (grip.merge(cov, on="PatientID")
              .merge(bmi, on="PatientID")
              .merge(img, on="PatientID"))
    for c in ["age", "sex", "BMI"] + race_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    results = []
    for pred in predictors:
        for roi in roi_cols:
            sub = df[[pred, roi] + covars].apply(pd.to_numeric, errors="coerce").dropna()
            # drop covariate columns with no variance (e.g. race dummy all 0)
            usecov = [c for c in covars if sub[c].nunique() > 1]
            if len(sub) < len(usecov) + 5:
                continue
            X = sm.add_constant(sub[[pred] + usecov])
            y = sub[roi]
            m = sm.OLS(y, X).fit()
            results.append({
                "metric": metric_name, "predictor": pred, "roi": roi,
                "n": int(len(sub)), "beta": m.params[pred],
                "p": m.pvalues[pred], "r_partial": np.sign(m.params[pred]) * np.sqrt(
                    m.tvalues[pred] ** 2 / (m.tvalues[pred] ** 2 + m.df_resid))
            })
    return pd.DataFrame(results)

all_res = pd.concat([run_metric(fa, "FA"), run_metric(md, "MD")], ignore_index=True)
all_res.to_csv(f"{here}/grip_dexterity_fa_md_results.csv", index=False)

sig = all_res[all_res["p"] < 0.05]
print(f"Total tests run: {len(all_res)}")
print(f"Merged sample n (typical): {all_res['n'].median():.0f}")
print(f"Covariates: age, sex, BMI, race (one-hot) | NO FDR/Bonferroni\n")
print(f"=== SIGNIFICANT (uncorrected p<0.05): {len(sig)} ===\n")

print("Breakdown by metric x predictor:")
tab = all_res.assign(sig=all_res.p < 0.05).groupby(["metric", "predictor"]).agg(
    tests=("p", "size"), n_sig=("sig", "sum")).reset_index()
print(tab.to_string(index=False))

print("\nExpected by chance (5%% of tests): %.1f" % (0.05 * len(all_res)))
print("\nStrongest 15 associations:")
print(sig.reindex(sig.p.abs().sort_values().index)
      .head(15)[["metric", "predictor", "roi", "n", "beta", "r_partial", "p"]]
      .to_string(index=False))
