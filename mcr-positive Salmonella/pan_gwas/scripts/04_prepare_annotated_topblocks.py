#!/usr/bin/env python3
"""Prepare Stage-1 top-block representatives with Panaroo annotations.
"""
from pathlib import Path
import argparse
import pandas as pd

ap=argparse.ArgumentParser()
ap.add_argument("--project-root",default=".")
args=ap.parse_args()
root=Path(args.project_root).resolve()
results=root/"pyseer_results"
stage1=results/"Candidate_Prioritization_Stage1_Hybrid_Enriched.tsv"
panaroo=root/"panaroo_result_moderate/gene_presence_absence.csv"
df=pd.read_csv(stage1,sep="\t",low_memory=False)
required=["variant","is_block_representative","freq_diff_pct","Hybrid_freq_pct"]
missing=[c for c in required if c not in df]
if missing: raise SystemExit("Missing Stage-1 columns: "+", ".join(missing))
df["freq_diff_pct"]=pd.to_numeric(df["freq_diff_pct"],errors="coerce")
df["Hybrid_freq_pct"]=pd.to_numeric(df["Hybrid_freq_pct"],errors="coerce")
mask=df["is_block_representative"]
if mask.dtype!="bool":
    mask=mask.astype(str).str.lower().eq("true")
top=df[mask & (df["freq_diff_pct"]>=20) & (df["Hybrid_freq_pct"]>=30)].copy()
topfile=results/"Candidate_Prioritization_Stage1_Hybrid_Enriched_TopBlocks.tsv"
top.to_csv(topfile,sep="\t",index=False)
pa=pd.read_csv(panaroo,low_memory=False)
for required_name in ["Gene","Non-unique Gene name","Annotation"]:
    if required_name not in pa.columns:
        raise SystemExit(f"Panaroo table lacks required column: {required_name}")
ann=pa[["Gene","Non-unique Gene name","Annotation"]].copy()
ann["Gene"]=ann["Gene"].astype(str)
if ann["Gene"].duplicated().any():
    raise SystemExit("Duplicate Gene keys in Panaroo table; verify original annotation mapping.")
ann=ann.rename(columns={
   "Gene":"variant",
   "Non-unique Gene name":"panaroo_nonunique_gene_name",
   "Annotation":"panaroo_annotation"})
out=top.merge(ann,on="variant",how="left",validate="many_to_one")
outfile=results/"Candidate_Prioritization_Stage1_Hybrid_Enriched_TopBlocks_Annotated.tsv"
out.to_csv(outfile,sep="\t",index=False)
print("Stage-1 top representatives:",len(top),
      "annotations matched:",int(out["panaroo_annotation"].notna().sum()))
print("Wrote:",topfile)
print("Wrote:",outfile)
