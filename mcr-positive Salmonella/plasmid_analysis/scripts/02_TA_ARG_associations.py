#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def load(path):
    df = pd.read_csv(path, sep="\t", dtype={"accession": str})
    if df["accession"].isna().any() or df["accession"].duplicated().any():
        raise ValueError(f"Nonunique/missing accessions: {path}")
    return df

def spearman_pair(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return np.nan, np.nan
    result = spearmanr(x, y)
    return float(result.statistic), float(result.pvalue)

def bh(p):
    p = np.asarray(p, dtype=float)
    good = np.isfinite(p)
    q = np.full(len(p), np.nan)
    if not good.any():
        return q
    idx = np.flatnonzero(good)
    ordered = idx[np.argsort(p[idx])]
    vals = p[ordered]
    adjusted = np.minimum.accumulate(
        (vals * len(vals) / np.arange(1, len(vals)+1))[::-1]
    )[::-1]
    q[ordered] = np.clip(adjusted, 0, 1)
    return q

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    data, out = Path(args.data), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta = load(data/"plasmid_metadata_937.tsv")
    arg = load(data/"plasmid_non_mcr_ARG_937.tsv")
    ta = load(data/"plasmid_TA_937.tsv")
    ids = set(meta.accession)
    if ids != set(arg.accession) or ids != set(ta.accession):
        raise ValueError("937 plasmid accession sets do not match.")

    d = (meta[["accession","sequence_length_bp","bacterial_species"]]
         .merge(arg[["accession","non_mcr_ARG_count"]], on="accession",
                validate="one_to_one")
         .merge(ta[["accession","figure_backbone_group",
                    "primary_TA_modules","primary_TA_family_list"]],
                on="accession", validate="one_to_one"))
    for col in ("sequence_length_bp", "non_mcr_ARG_count", "primary_TA_modules"):
        d[col] = pd.to_numeric(d[col], errors="raise")

    metrics = [
        ("TA_modules_vs_non_mcr_ARGs", "primary_TA_modules","non_mcr_ARG_count"),
        ("TA_modules_vs_plasmid_length", "primary_TA_modules","sequence_length_bp"),
        ("non_mcr_ARGs_vs_plasmid_length","non_mcr_ARG_count","sequence_length_bp"),
    ]
    summary = []
    for label, x, y in metrics:
        rho, p = spearman_pair(d[x], d[y])
        summary.append({"analysis":label, "n":len(d),
                        "rho":rho, "p_two_sided":p})
    pd.DataFrame(summary).to_csv(out/"TA_ARG_overall_spearman.tsv",
                                 sep="\t",index=False)

    strata = []
    for name, subset in d.groupby("figure_backbone_group", dropna=False):
        rho, p = spearman_pair(subset["primary_TA_modules"],
                               subset["non_mcr_ARG_count"])
        strata.append({"backbone":name, "n":len(subset),
                       "rho":rho, "p_two_sided":p})
    strata = pd.DataFrame(strata).sort_values("backbone")
    # BH applied only to evaluable (finite-p) tests, as in manuscript Methods.
    strata["q_BH_evaluable_strata"] = bh(strata["p_two_sided"])
    strata.to_csv(out/"TA_ARG_by_backbone.tsv",sep="\t",index=False)


    print(f"Plasmids: {len(d)}")
    print(f"Evaluable backbone strata: {strata.p_two_sided.notna().sum()}")
    print(f"Results: {out.resolve()}")

if __name__ == "__main__":
    main()
