#!/usr/bin/env python3

import argparse
from pathlib import Path
import re
import pandas as pd

def load(path):
    frame = pd.read_csv(path, sep="\t", dtype={"accession":str})
    if "accession" not in frame:
        raise ValueError(f"Missing accession column in {path}")
    if frame["accession"].isna().any() or frame["accession"].duplicated().any():
        raise ValueError(f"Missing or duplicated plasmid accession in {path}")
    return frame

def split_terms(raw):
    if pd.isna(raw) or str(raw).strip().lower() in ("", "nan", "none"):
        return []
    return [x.strip() for x in str(raw).split(";") if x.strip()]

def has_gene_family(genes, family):
    # e.g. blaNDM-1 -> NDM, blaKPC-2 -> KPC.
    pattern = re.compile(r"^bla" + family + r"(?:$|[^a-z])", flags=re.I)
    return any(pattern.search(x) for x in genes)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data", default="data")
    p.add_argument("--out", default="results")
    args = p.parse_args()
    data, out = Path(args.data), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    meta = load(data / "plasmid_metadata_937.tsv")
    arg = load(data / "plasmid_non_mcr_ARG_937.tsv")
    if set(meta["accession"]) != set(arg["accession"]):
        raise ValueError("Plasmid accessions are not identical in metadata and ARG data.")
    d = meta.merge(arg, on="accession", validate="one_to_one", suffixes=("", "_arg"))
    d["non_mcr_ARG_count"] = pd.to_numeric(d["non_mcr_ARG_count"], errors="raise")
    d["sequence_length_bp"] = pd.to_numeric(d["sequence_length_bp"], errors="raise")

    # Carrying at least one or ten non-mcr ARGs, grouped by recorded host species.
    d["has_additional_ARG"] = d["non_mcr_ARG_count"] >= 1
    d["has_ge10_additional_ARG"] = d["non_mcr_ARG_count"] >= 10
    grp = d.groupby("bacterial_species", dropna=False)
    summary = grp.agg(
        n_plasmids=("accession", "size"),
        n_with_additional_ARG=("has_additional_ARG", "sum"),
        n_with_ge10_ARG=("has_ge10_additional_ARG", "sum"),
        median_additional_ARG=("non_mcr_ARG_count", "median"),
    ).reset_index()
    summary["pct_additional_ARG"] = (100*summary["n_with_additional_ARG"] /
                                     summary["n_plasmids"]).round(2)
    summary["pct_ge10_ARG"] = (100*summary["n_with_ge10_ARG"] /
                               summary["n_plasmids"]).round(2)
    summary = summary.sort_values(["n_plasmids", "bacterial_species"],
                                  ascending=[False, True])
    summary.to_csv(out/"plasmid_host_species_summary.tsv", sep="\t", index=False)

    
    families = ["NDM", "VIM", "IMP", "KPC"]
    for family in families:
        d["bla" + family] = d["non_mcr_ARG_combination"].map(
            lambda x: has_gene_family(split_terms(x), family))
    d["any_four_carbapenemase_families"] = d[["bla"+x for x in families]].any(axis=1)
    d["any_tetX"] = d["non_mcr_ARG_combination"].map(
        lambda raw: any(re.match(r"^tet\(X[0-9]*\)", gene, flags=re.I)
                        for gene in split_terms(raw)))
    selected = d[d["any_four_carbapenemase_families"] | d["any_tetX"]].copy()
    fields = ["accession", "bacterial_species", "target_mcr_families",
              "known_replicon_combination", "non_mcr_ARG_count",
              "non_mcr_ARG_combination", "blaNDM", "blaVIM", "blaIMP",
              "blaKPC", "any_tetX"]
    selected[fields].to_csv(out/"selected_co_resistance_plasmids.tsv",
                            sep="\t", index=False)

    by_backbone = d.groupby("broad_backbone_group", dropna=False).agg(
        n_plasmids=("accession", "size"),
        median_length_bp=("sequence_length_bp", "median"),
        median_non_mcr_ARG=("non_mcr_ARG_count", "median")
    ).reset_index()
    by_backbone.to_csv(out/"plasmid_backbone_summary.tsv",sep="\t",index=False)
    lines = [
        f"Plasmids: {len(d)}",
        f"Unique accessions: {d.accession.nunique()}",
        f"Non-mcr ARG >=1: {d.has_additional_ARG.sum()}",
        f"Non-mcr ARG >=10: {d.has_ge10_additional_ARG.sum()}",
        f"Any of blaNDM/blaVIM/blaIMP/blaKPC: {d.any_four_carbapenemase_families.sum()}",
        f"Any tet(X) marker: {d.any_tetX.sum()}",
    ] + [f"bla{family}: {d['bla'+family].sum()}" for family in families]
    (out/"summary_counts.txt").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"Output: {out.resolve()}")

if __name__ == "__main__":
    main()
