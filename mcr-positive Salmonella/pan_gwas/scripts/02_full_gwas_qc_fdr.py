#!/usr/bin/env python3
import numpy as np
import pandas as pd

f = "../pyseer_results/panGWAS_full_pyseer.tsv"
df = pd.read_csv(f, sep="\t", low_memory=False)

df["notes_clean"] = df["notes"].fillna("none").replace("", "none")
df["lrt-pvalue"] = pd.to_numeric(df["lrt-pvalue"], errors="coerce")
df["beta"] = pd.to_numeric(df["beta"], errors="coerce")

df["p_valid"] = (
    df["lrt-pvalue"].notna()
    & np.isfinite(df["lrt-pvalue"])
)

print("=== notes × p_valid ===")
print(pd.crosstab(df["notes_clean"], df["p_valid"]))

# Primary QC: exclude diagnostic-note rows and nonfinite LRT P values
main = df[
    (df["notes_clean"] == "none") &
    (df["p_valid"])
].copy()

# Benjamini-Hochberg FDR correction
p = main["lrt-pvalue"].to_numpy()
m = len(p)

order = np.argsort(p)
ranked = p[order]

q_ranked = ranked * m / np.arange(1, m + 1)
q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
q_ranked = np.minimum(q_ranked, 1.0)

q = np.empty(m)
q[order] = q_ranked

main["q_value_BH_mainQC"] = q
main["odds_ratio"] = np.exp(
    np.clip(main["beta"], -700, 700)
)

main = main.sort_values(
    ["q_value_BH_mainQC", "lrt-pvalue"]
)

main.to_csv(
    "panGWAS_full_pyseer_mainQC_notesNone_with_FDR.tsv",
    sep="\t", index=False
)

main[
    main["q_value_BH_mainQC"] < 0.05
].to_csv(
    "panGWAS_full_pyseer_mainQC_FDR0.05.tsv",
    sep="\t", index=False
)

main[
    main["q_value_BH_mainQC"] < 0.10
].to_csv(
    "panGWAS_full_pyseer_mainQC_FDR0.10.tsv",
    sep="\t", index=False
)

flagged = df[df["notes_clean"] != "none"].copy()
flagged.to_csv(
    "panGWAS_full_pyseer_flagged_notes.tsv",
    sep="\t", index=False
)

invalid = df[~df["p_valid"]].copy()
invalid.to_csv(
    "panGWAS_full_pyseer_invalid_pvalue.tsv",
    sep="\t", index=False
)

print("Main QC passed:", len(main))
print("BH-FDR < 0.05:",
      (main["q_value_BH_mainQC"] < 0.05).sum())
print("BH-FDR < 0.10:",
      (main["q_value_BH_mainQC"] < 0.10).sum())
