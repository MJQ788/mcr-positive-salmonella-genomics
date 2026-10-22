#!/usr/bin/env python3

import argparse
import hashlib
import re
from pathlib import Path

import pandas as pd
import numpy as np

p=argparse.ArgumentParser()
p.add_argument("--project-root",default=".")
args=p.parse_args()
base=Path(args.project_root).resolve()
rtab_file=base/"panGWAS_input/pyseer_gene_presence_absence_min1pct.Rtab"
panaroo_gpa_file=base/"panaroo_result_moderate/gene_presence_absence.csv"
gene_data_file=base/"panaroo_result_moderate/gene_data.csv"
outdir=base/"pyseer_results/filtered_nondefinition_gwas"
outdir.mkdir(parents=True,exist_ok=True)
feature_annot_out=outdir/"AllRtab_FeatureAnnotation_DefinitionCargoFlags.tsv"
exclude_out=outdir/"DefinitionCargo_ExcludedFeatures.tsv"
summary_out=outdir/"DefinitionCargo_FilteringSummary.txt"
filtered_conservative_out=base/"panGWAS_input/pyseer_gene_presence_absence_min1pct.filtered_conservative.Rtab"
filtered_strict_out=base/"panGWAS_input/pyseer_gene_presence_absence_min1pct.filtered_strict.Rtab"
for f in [rtab_file,panaroo_gpa_file,gene_data_file]:
    if not f.exists():
        raise SystemExit(f"FATAL ERROR: file not found: {f}")
rtab=pd.read_csv(rtab_file,sep="\t",low_memory=False)
feature_col=rtab.columns[0]
rtab[feature_col]=rtab[feature_col].astype(str)
features=rtab[feature_col].tolist()
sample_cols=list(rtab.columns[1:])
gpa=pd.read_csv(panaroo_gpa_file,low_memory=False)
gd=pd.read_csv(gene_data_file,low_memory=False)

possible_gene_cols=["gene","group","cluster","clustering_id"]
gpa_gene_col=next((c for c in gpa.columns if c.lower() in possible_gene_cols),gpa.columns[0])
gpa_meta_cols=[]
for c in gpa.columns:
    lc=c.lower()
    if c==gpa_gene_col or lc in [
       "annotation","description","non-unique gene name","gene name",
       "no. isolates","no. sequences","avg sequences per isolate",
       "genome fragment","order within fragment","accessory fragment",
       "accessory order with fragment","qc","min group size nuc",
       "max group size nuc","avg group size nuc"]:
        gpa_meta_cols.append(c)
gpa_small=gpa[list(dict.fromkeys(gpa_meta_cols))].copy()
gpa_small[gpa_gene_col]=gpa_small[gpa_gene_col].astype(str)
gpa_annot_map={}
for _,row in gpa_small.iterrows():
    gene=str(row[gpa_gene_col])
    gpa_annot_map[gene]=" | ".join(
        f"{c}={row[c]}" for c in gpa_small.columns if c!=gpa_gene_col and pd.notna(row[c])
    )

if "clustering_id" not in gd.columns:
    raise SystemExit("FATAL ERROR: gene_data.csv does not contain clustering_id.")
for c in ["clustering_id","gene_name","description","annotation_id"]:
    if c not in gd.columns:
        gd[c]=""
gd["clustering_id"]=gd["clustering_id"].astype(str)
gd_summary=(
 gd.groupby("clustering_id")
 .agg({c:lambda x:";".join(sorted(set(
     str(v) for v in x if pd.notna(v) and str(v)!="nan")))[:20]
       for c in ["gene_name","description","annotation_id"]})
 .reset_index())
gd_map={}
for _,row in gd_summary.iterrows():
    gd_map[str(row["clustering_id"])]=(
        f"gene_name={row['gene_name']} | "
        f"description={row['description']} | "
        f"annotation_id={row['annotation_id']}"
    )

def pattern_hash(vals):
    s="".join(map(str,vals.astype(int).tolist()))
    return hashlib.md5(s.encode()).hexdigest()
rtab_matrix=rtab[sample_cols].fillna(0).astype(int)
pattern_ids=rtab_matrix.apply(pattern_hash,axis=1)
rtab["CoPresence_Block_AllRtab"]=pattern_ids
block_sizes=rtab.groupby("CoPresence_Block_AllRtab").size().to_dict()

user_arg_cargo_terms=["fosA","oqxB","oqx","bepE","hph","mphL","mph","neo"]
user_rep_methylation_terms=["dnaC","dcm"]
amr_terms=[
    r"\bfosA\b",r"\boqx",r"\bhph\b",r"\bmph",r"\bneo\b",
    r"\baph",r"\baad",r"\baac",r"\bant",r"\bstrA\b",r"\bstrB\b",
    r"\bbla",r"\btet",r"\bsul",r"\bdfr",r"\bqnr",r"\bcat\b",r"\bcml",r"\bfloR\b",
    r"\berm",r"\bere",r"\barr",r"\bmdf",r"\bmef",
    r"aminoglycoside",r"beta-lactam",r"betalactam",r"tetracycline",
    r"sulfonamide",r"trimethoprim",r"chloramphenicol",r"macrolide",
    r"quinolone",r"fosfomycin",r"hygromycin",r"neomycin",
    r"antibiotic resistance",r"drug resistance",r"multidrug",r"efflux"
]
vf_terms=[
    r"virulence",r"\bspi[-_ ]?1\b",r"\bspi[-_ ]?2\b",
    r"\binv[A-Z]?",r"\bsip[A-Z]?",r"\bsop[A-Z]?",r"\bsse[A-Z]?",
    r"\bssa[A-Z]?",r"\bssr[A-Z]?",r"\bprg[A-Z]?",r"\bhil[A-Z]?",
    r"\bstn\b",r"\bspv[A-Z]?",r"\biro[A-Z]?",r"\biut[A-Z]?"
]
mge_rep_methyl_terms=[
    r"\bdnaC\b",r"\bdcm\b",r"methyltransferase",r"methylase",
    r"restriction",r"modification",r"replication protein",
    r"replicase",r"plasmid replication",r"partition protein"
]
mge_context_terms=[
    r"transposase",r"integrase",r"recombinase",r"insertion sequence",
    r"\bIS[0-9A-Za-z_-]*\b",r"plasmid",r"conjug",r"relaxase",
    r"phage",r"prophage",r"bacteriophage",r"tail fiber",r"capsid",
    r"terminase",r"portal protein",r"toxin",r"antitoxin"
]
def match_any(text,patterns):
    return any(re.search(regex,str(text),flags=re.IGNORECASE) for regex in patterns)
def exact_name_match(feature,terms):
    parts=[p.strip() for p in re.split(r"~~~|;|,|\s+",str(feature)) if p.strip()]
    return any(term.lower()==v.lower() or term.lower() in v.lower()
               for term in terms for v in parts)
rows=[]
for i,feature in enumerate(features):
    pieces=[f"feature={feature}"]
    if feature in gpa_annot_map: pieces.append(gpa_annot_map[feature])
    if feature in gd_map: pieces.append(gd_map[feature])
    for comp in re.split(r"~~~",feature):
        comp=comp.strip()
        if comp and comp!=feature:
            if comp in gpa_annot_map: pieces.append(f"component_{comp}: {gpa_annot_map[comp]}")
            if comp in gd_map: pieces.append(f"component_{comp}: {gd_map[comp]}")
    annot_text=" | ".join(pieces)
    rows.append({
        "feature_id":feature,
        "CoPresence_Block_AllRtab":rtab.loc[i,"CoPresence_Block_AllRtab"],
        "CoPresence_Block_Size_AllRtab":block_sizes[rtab.loc[i,"CoPresence_Block_AllRtab"]],
        "annotation_text":annot_text,
        "flag_user_ARG_cargo_terms":exact_name_match(feature,user_arg_cargo_terms),
        "flag_user_replication_methylation_terms":exact_name_match(feature,user_rep_methylation_terms),
        "flag_AMR_resistance_efflux":match_any(annot_text,amr_terms),
        "flag_VF_virulence":match_any(annot_text,vf_terms),
        "flag_MGE_replication_methylation_review":match_any(annot_text,mge_rep_methyl_terms),
        "flag_MGE_context_review":match_any(annot_text,mge_context_terms)
    })
feat=pd.DataFrame(rows)
feat["conservative_seed_exclude"]=(
    feat["flag_user_ARG_cargo_terms"]|
    feat["flag_AMR_resistance_efflux"]|
    feat["flag_VF_virulence"])
feat["strict_seed_exclude"]=(
    feat["conservative_seed_exclude"]|
    feat["flag_user_replication_methylation_terms"]|
    feat["flag_MGE_replication_methylation_review"])
conservative_blocks=set(feat.loc[feat["conservative_seed_exclude"],"CoPresence_Block_AllRtab"])
strict_blocks=set(feat.loc[feat["strict_seed_exclude"],"CoPresence_Block_AllRtab"])
feat["conservative_block_exclude"]=feat["CoPresence_Block_AllRtab"].isin(conservative_blocks)
feat["strict_block_exclude"]=feat["CoPresence_Block_AllRtab"].isin(strict_blocks)
def reason_list(row,mode):
    reasons=[]
    for flag,reason in [
       ("flag_user_ARG_cargo_terms","user_listed_ARG_AMR_efflux_cargo"),
       ("flag_AMR_resistance_efflux","AMR_resistance_efflux_keyword"),
       ("flag_VF_virulence","VF_virulence_keyword"),
       ("flag_user_replication_methylation_terms","user_listed_replication_methylation_review"),
       ("flag_MGE_replication_methylation_review","MGE_replication_methylation_keyword")
    ]:
        if row[flag]: reasons.append(reason)
    if mode=="conservative" and row["conservative_block_exclude"] and not row["conservative_seed_exclude"]:
        reasons.append("co_presence_block_with_conservative_seed")
    if mode=="strict" and row["strict_block_exclude"] and not row["strict_seed_exclude"]:
        reasons.append("co_presence_block_with_strict_seed")
    return ";".join(reasons)
feat["conservative_exclusion_reason"]=feat.apply(lambda r:reason_list(r,"conservative"),axis=1)
feat["strict_exclusion_reason"]=feat.apply(lambda r:reason_list(r,"strict"),axis=1)
feat.to_csv(feature_annot_out,sep="\t",index=False)
excluded=feat[feat["conservative_block_exclude"]|feat["strict_block_exclude"]].copy()
excluded.to_csv(exclude_out,sep="\t",index=False)
conservative_keep=set(feat.loc[~feat["conservative_block_exclude"],"feature_id"])
strict_keep=set(feat.loc[~feat["strict_block_exclude"],"feature_id"])
rtab_cons=rtab[rtab[feature_col].isin(conservative_keep)].drop(columns=["CoPresence_Block_AllRtab"])
rtab_strict=rtab[rtab[feature_col].isin(strict_keep)].drop(columns=["CoPresence_Block_AllRtab"])
rtab_cons.to_csv(filtered_conservative_out,sep="\t",index=False)
rtab_strict.to_csv(filtered_strict_out,sep="\t",index=False)
with open(summary_out,"w") as f:
    f.write("=== Definition/cargo filtering summary ===\n\n")
    f.write(f"Input Rtab: {rtab_file}\nTotal features: {len(feat)}\nSamples: {len(sample_cols)}\n\n")
    f.write("=== Seed flags ===\n")
    for col in ["flag_user_ARG_cargo_terms",
                "flag_user_replication_methylation_terms",
                "flag_AMR_resistance_efflux","flag_VF_virulence",
                "flag_MGE_replication_methylation_review",
                "flag_MGE_context_review"]:
        f.write(f"{col}: {int(feat[col].sum())}\n")
    f.write("\n=== Conservative filtering ===\n")
    f.write(f"Conservative seed excluded features: {int(feat['conservative_seed_exclude'].sum())}\n")
    f.write(f"Conservative excluded co-presence blocks: {len(conservative_blocks)}\n")
    f.write(f"Conservative block-excluded features: {int(feat['conservative_block_exclude'].sum())}\n")
    f.write(f"Conservative retained features: {len(rtab_cons)}\n\n")
    f.write("=== Strict filtering ===\n")
    f.write(f"Strict seed excluded features: {int(feat['strict_seed_exclude'].sum())}\n")
    f.write(f"Strict excluded co-presence blocks: {len(strict_blocks)}\n")
    f.write(f"Strict block-excluded features: {int(feat['strict_block_exclude'].sum())}\n")
    f.write(f"Strict retained features: {len(rtab_strict)}\n")
print("Conservative",len(rtab_cons),"Strict",len(rtab_strict),"Samples",len(sample_cols))
