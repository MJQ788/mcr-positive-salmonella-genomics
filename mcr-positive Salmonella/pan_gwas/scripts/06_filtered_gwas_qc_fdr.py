#!/usr/bin/env python3
"""BH-FDR processing of existing conservative/strict Pyseer output files.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ap=argparse.ArgumentParser()
ap.add_argument("--project-root",default=".")
args=ap.parse_args()
root=Path(args.project_root).resolve()
folder=root/"pyseer_results/filtered_nondefinition_gwas"

def bh_fdr(p):
    p=np.asarray(p,dtype=float)
    n=len(p)
    if n==0: return np.array([],dtype=float)
    order=np.argsort(p)
    ranks=np.arange(1,n+1,dtype=float)
    raw=p[order]*n/ranks
    adj=np.minimum.accumulate(raw[::-1])[::-1]
    q=np.empty(n,dtype=float)
    q[order]=np.minimum(adj,1.0)
    return q

for mode in ("conservative","strict"):
    path=folder/f"panGWAS_filtered_{mode}_pyseer.tsv"
    df=pd.read_csv(path,sep="\t",low_memory=False)
    if "lrt-pvalue" not in df or "notes" not in df:
        raise SystemExit(f"{path} lacks lrt-pvalue or notes")
    p=pd.to_numeric(df["lrt-pvalue"],errors="coerce")
    notes=df["notes"].fillna("").astype(str).str.strip()
    good=p.notna() & np.isfinite(p) & (notes.eq("")|notes.str.lower().eq("none"))
    df["mainQC_pass"]=good
    df["q_value_BH_mainQC"]=np.nan
    df.loc[good,"q_value_BH_mainQC"]=bh_fdr(p[good].to_numpy())
    df["FDR05"]=good & (df["q_value_BH_mainQC"]<0.05)
    df.to_csv(folder/f"panGWAS_filtered_{mode}_mainQC_FDR.tsv",sep="\t",index=False)
    print(mode, "rows",len(df),"QC tests",int(good.sum()),"q<0.05",int(df["FDR05"].sum()))
