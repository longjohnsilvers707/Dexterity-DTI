import pandas as pd

here = r"c:\Users\bandd\Downloads\52426_DepressionWork\DexterityDTI study"
r = pd.read_csv(f"{here}/grip_dexterity_fa_md_results.csv")
r["sig"] = r["p"] < 0.05
sig = r[r["sig"]].copy()
preds = ["grip_dominant", "grip_nondominant", "dexterity_dominant", "dexterity_nondominant"]

def line(s=""): print(s)

# ============================================================
# 1) ROIs common to BOTH FA and MD (same ROI, same predictor)
# ============================================================
line("="*70)
line("1) ROIs significant in BOTH FA *and* MD (same predictor)")
line("="*70)
both_all = set()
for p in preds:
    fa = set(sig[(sig.predictor==p)&(sig.metric=="FA")].roi)
    md = set(sig[(sig.predictor==p)&(sig.metric=="MD")].roi)
    common = sorted(fa & md)
    both_all |= set((p, roi) for roi in common)
    line(f"\n{p}: {len(common)} ROIs in both FA & MD")
    for roi in common:
        sub = r[(r.predictor==p)&(r.roi==roi)]
        fb = sub[sub.metric=="FA"].iloc[0]; mb = sub[sub.metric=="MD"].iloc[0]
        line(f"   {roi:<48} FA r={fb.r_partial:+.3f} p={fb.p:.4f} | MD r={mb.r_partial:+.3f} p={mb.p:.4f}")

# ============================================================
# 2) ALL significant ROIs, listed per metric x predictor
# ============================================================
line("\n"+"="*70)
line("2) ALL significant ROIs (uncorrected p<0.05), per metric x predictor")
line("="*70)
for metric in ["FA", "MD"]:
    for p in preds:
        sub = sig[(sig.metric==metric)&(sig.predictor==p)].sort_values("p")
        line(f"\n--- {metric} / {p}  ({len(sub)} ROIs) ---")
        for _, row in sub.iterrows():
            line(f"   {row.roi:<48} r={row.r_partial:+.3f}  p={row.p:.4f}")

# ============================================================
# 3) ROIs consistent across ALL 4 predictors
# ============================================================
line("\n"+"="*70)
line("3) ROIs significant across ALL 4 grip/dexterity variables")
line("="*70)
for metric in ["FA", "MD", "either(FA or MD)"]:
    if metric == "either(FA or MD)":
        sets = [set(sig[sig.predictor==p].roi) for p in preds]
    else:
        sets = [set(sig[(sig.metric==metric)&(sig.predictor==p)].roi) for p in preds]
    common4 = sorted(set.intersection(*sets))
    line(f"\n[{metric}] ROIs hitting all 4 predictors: {len(common4)}")
    for roi in common4:
        line(f"   {roi}")

# also: ROIs significant in all 4 AND in both metrics
line("\n-- Bonus: ROI x ... that are sig in all 4 predictors in BOTH FA and MD --")
fa4 = set.intersection(*[set(sig[(sig.metric=="FA")&(sig.predictor==p)].roi) for p in preds])
md4 = set.intersection(*[set(sig[(sig.metric=="MD")&(sig.predictor==p)].roi) for p in preds])
line(f"   {sorted(fa4 & md4) if (fa4&md4) else 'none'}")
