#!/usr/bin/env python3
"""Stage-1 candidate prioritization: frequency calculations and exact co-presence block representatives."""
import argparse
import pandas as pd
import numpy as np
import hashlib
import sys
from pathlib import Path

parser=argparse.ArgumentParser()
parser.add_argument("--project-root",default=".")
args=parser.parse_args()
base_dir=Path(args.project_root).resolve()
pyseer_dir=base_dir/"pyseer_results"
gwas_input_dir=base_dir/"panGWAS_input"
annotated_sig_file=gwas_input_dir/"panGWAS_full_pyseer_mainQC_FDR0.05.tsv"
pheno_file=gwas_input_dir/"pyseer_Hybrid_phenotypes.tsv"
rtab_file=gwas_input_dir/"pyseer_gene_presence_absence_min1pct.Rtab"
pyseer_dir.mkdir(parents=True,exist_ok=True)
out_all=pyseer_dir/"Candidate_Prioritization_Stage1_All.tsv"
out_enriched=pyseer_dir/"Candidate_Prioritization_Stage1_Hybrid_Enriched.tsv"
out_depleted=pyseer_dir/"Candidate_Prioritization_Stage1_Hybrid_Depleted.tsv"
out_block_summary=pyseer_dir/"Candidate_CoPresence_Block_Summary.tsv"

print("=== Stage-1 candidate prioritization started ===")
print(f"Python executable: {sys.executable}")
print(f"pandas version: {pd.__version__}")
print(f"GWAS input file: {annotated_sig_file}")
for f in [annotated_sig_file,pheno_file,rtab_file]:
    if not f.is_file(): raise SystemExit(f"FATAL ERROR: file not found: {f}")

sig_df=pd.read_csv(annotated_sig_file,sep="\t",low_memory=False)
pheno=pd.read_csv(pheno_file,sep="\t")
rtab=pd.read_csv(rtab_file,sep="\t",low_memory=False)
required_sig_cols=["variant","beta","lrt-pvalue"]
required_pheno_cols=["samples","Hybrid"]
for col in required_sig_cols:
    if col not in sig_df.columns: raise SystemExit(f"FATAL ERROR: missing column in GWAS table: {col}")
for col in required_pheno_cols:
    if col not in pheno.columns: raise SystemExit(f"FATAL ERROR: missing column in phenotype file: {col}")
sig_df["variant"]=sig_df["variant"].astype(str)
sig_df["beta"]=pd.to_numeric(sig_df["beta"],errors="coerce")
sig_df["lrt-pvalue"]=pd.to_numeric(sig_df["lrt-pvalue"],errors="coerce")
sig_df=sig_df.dropna(subset=["variant","beta","lrt-pvalue"]).copy()

pheno["samples"]=pheno["samples"].astype(str)
pheno["Hybrid"]=pd.to_numeric(pheno["Hybrid"],errors="coerce")
pheno=pheno[pheno["Hybrid"].isin([0,1])].copy()
hybrid_samples_raw=pheno.loc[pheno["Hybrid"]==1,"samples"].tolist()
non_hybrid_samples_raw=pheno.loc[pheno["Hybrid"]==0,"samples"].tolist()

gene_col=rtab.columns[0]
rtab[gene_col]=rtab[gene_col].astype(str)
sample_cols=[str(x) for x in rtab.columns[1:]]
rtab.columns=[gene_col]+sample_cols
rtab[sample_cols]=rtab[sample_cols].apply(pd.to_numeric,errors="coerce").fillna(0)
rtab[sample_cols]=(rtab[sample_cols]>0).astype(np.int8)
rtab_sig=rtab[rtab[gene_col].isin(sig_df["variant"])].copy()
missing_variants=sorted(set(sig_df["variant"])-set(rtab_sig[gene_col]))
if missing_variants:
    print(f"WARNING: {len(missing_variants)} significant variants not found in Rtab.")
    print("First ten missing:",missing_variants[:10])
rtab_sig=rtab_sig.set_index(gene_col)
hybrid_samples=sorted(s for s in hybrid_samples_raw if s in rtab_sig.columns)
non_hybrid_samples=sorted(s for s in non_hybrid_samples_raw if s in rtab_sig.columns)
n_hybrid=len(hybrid_samples)
n_non_hybrid=len(non_hybrid_samples)
if n_hybrid==0 or n_non_hybrid==0:
    raise SystemExit("FATAL ERROR: no matched Hybrid or non-Hybrid samples between phenotype and Rtab.")

freq_hybrid=rtab_sig[hybrid_samples].sum(axis=1)
freq_non_hybrid=rtab_sig[non_hybrid_samples].sum(axis=1)
freq_df=pd.DataFrame({
    "variant":freq_hybrid.index,
    "Hybrid_count":freq_hybrid.values,
    "NonHybrid_count":freq_non_hybrid.values,
    "Hybrid_freq_pct":np.round(freq_hybrid.values/n_hybrid*100,2),
    "NonHybrid_freq_pct":np.round(freq_non_hybrid.values/n_non_hybrid*100,2)
})
freq_df["freq_diff_pct"]=np.round(freq_df["Hybrid_freq_pct"]-freq_df["NonHybrid_freq_pct"],2)
sig_df=sig_df.merge(freq_df,on="variant",how="left")
sig_df["association_direction"]=np.where(
    sig_df["beta"]>0,"Hybrid-enriched",
    np.where(sig_df["beta"]<0,"Hybrid-depleted","No-direction"))
if "odds_ratio" not in sig_df.columns:
    sig_df["odds_ratio"]=np.exp(sig_df["beta"])
all_matched_samples=hybrid_samples+non_hybrid_samples
def hash_presence_vector(row):
    vec="".join(row[all_matched_samples].astype(int).astype(str).tolist())
    return hashlib.md5(vec.encode("utf-8")).hexdigest()
presence_for_blocks=rtab_sig.loc[
    rtab_sig.index.intersection(sig_df["variant"]),all_matched_samples].copy()
block_signature_df=pd.DataFrame({
    "variant":presence_for_blocks.index,
    "presence_vector_hash":presence_for_blocks.apply(hash_presence_vector,axis=1)})
sig_df=sig_df.merge(block_signature_df,on="variant",how="left")
block_rank=(sig_df.dropna(subset=["presence_vector_hash"])
    .groupby("presence_vector_hash")["lrt-pvalue"].min()
    .sort_values().reset_index())
block_rank["CoPresence_Block_ID"]=[
    f"CoPresence_Block_{i+1:04d}" for i in range(len(block_rank))]
sig_df=sig_df.merge(block_rank,on="presence_vector_hash",how="left",suffixes=("","_block_min"))
block_counts=(sig_df.groupby("CoPresence_Block_ID",dropna=True)["variant"]
    .count().reset_index().rename(columns={"variant":"Genes_in_CoPresence_Block"}))
sig_df=sig_df.merge(block_counts,on="CoPresence_Block_ID",how="left")
sig_df["is_block_representative"]=False
idx_rep=(sig_df.dropna(subset=["CoPresence_Block_ID"])
         .sort_values(["CoPresence_Block_ID","lrt-pvalue","variant"])
         .groupby("CoPresence_Block_ID").head(1).index)
sig_df.loc[idx_rep,"is_block_representative"]=True
sig_df["candidate_stage1_class"]=np.select([
    (sig_df["association_direction"]=="Hybrid-enriched") & sig_df["is_block_representative"],
    (sig_df["association_direction"]=="Hybrid-enriched") & (~sig_df["is_block_representative"]),
    (sig_df["association_direction"]=="Hybrid-depleted") & sig_df["is_block_representative"],
    (sig_df["association_direction"]=="Hybrid-depleted") & (~sig_df["is_block_representative"])
],[
    "Hybrid-enriched block representative",
    "Hybrid-enriched co-presence block member",
    "Hybrid-depleted block representative",
    "Hybrid-depleted co-presence block member"
],default="Unclassified")
sig_df=sig_df.sort_values(["association_direction","lrt-pvalue","CoPresence_Block_ID","variant"])
enriched_df=sig_df[sig_df["association_direction"]=="Hybrid-enriched"].copy()
depleted_df=sig_df[sig_df["association_direction"]=="Hybrid-depleted"].copy()
block_summary=(sig_df.dropna(subset=["CoPresence_Block_ID"])
   .sort_values(["CoPresence_Block_ID","lrt-pvalue","variant"])
   .groupby("CoPresence_Block_ID")
   .agg(Genes_in_Block=("variant","count"),
        Representative_gene=("variant",lambda x:x.iloc[0]),
        Min_lrt_pvalue=("lrt-pvalue","min"),
        Mean_beta=("beta","mean"),
        Direction=("association_direction",lambda x:x.mode().iloc[0] if len(x.mode())>0 else "mixed"),
        Hybrid_count=("Hybrid_count","max"),
        NonHybrid_count=("NonHybrid_count","max"),
        Hybrid_freq_pct=("Hybrid_freq_pct","max"),
        NonHybrid_freq_pct=("NonHybrid_freq_pct","max"),
        freq_diff_pct=("freq_diff_pct","max"))
   .reset_index().sort_values("Min_lrt_pvalue"))
sig_df.to_csv(out_all,sep="\t",index=False)
enriched_df.to_csv(out_enriched,sep="\t",index=False)
depleted_df.to_csv(out_depleted,sep="\t",index=False)
block_summary.to_csv(out_block_summary,sep="\t",index=False)
print("=== Stage-1 candidate prioritization finished ===")
print("Total significant genes:",len(sig_df))
print("Co-presence blocks:",sig_df["CoPresence_Block_ID"].nunique())
print("Hybrid-enriched genes:",len(enriched_df))
print("Hybrid-depleted genes:",len(depleted_df))
print("Block representatives:",sig_df["is_block_representative"].sum())
