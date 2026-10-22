#!/usr/bin/env python3
import pandas as pd
import re

infile = (
    "pyseer_results/"
    "Candidate_Prioritization_Stage1_Hybrid_Enriched_"
    "TopBlocks_Annotated.tsv"
)

outfile = (
    "pyseer_results/"
    "Candidate_Prioritization_Stage1_Hybrid_Enriched_"
    "TopBlocks_Annotated_Reclassified_v2.tsv"
)

df = pd.read_csv(infile, sep="\t", low_memory=False)

for col in [
    "variant",
    "panaroo_nonunique_gene_name",
    "panaroo_annotation"
]:
    if col not in df.columns:
        df[col] = ""

text = (
    df["variant"].fillna("").astype(str) + " " +
    df["panaroo_nonunique_gene_name"].fillna("").astype(str) + " " +
    df["panaroo_annotation"].fillna("").astype(str)
)

# AMR/efflux: use gene-name boundaries
amr_gene_pattern = re.compile(
    r"(^|[^A-Za-z0-9])("
    r"fosA|oqx[A-Za-z0-9_]*|tet[A-Za-z0-9_]*|bla[A-Za-z0-9_]*|"
    r"aac[A-Za-z0-9_()'\-]*|aph[A-Za-z0-9_()'\-]*|"
    r"aad[A-Za-z0-9_]*|ant[A-Za-z0-9_()'\-]*|"
    r"sul[A-Za-z0-9_]*|dfr[A-Za-z0-9_]*|floR|"
    r"qnr[A-Za-z0-9_]*|mcr[A-Za-z0-9_]*|"
    r"erm[A-Za-z0-9_]*|mph[A-Za-z0-9_]*|"
    r"lnu[A-Za-z0-9_]*|hph|neo[A-Za-z0-9_]*|"
    r"str[A-Za-z0-9_]*|cml[A-Za-z0-9_]*|"
    r"arr[A-Za-z0-9_]*|rmt[A-Za-z0-9_]*|"
    r"arm[A-Za-z0-9_]*|ble[A-Za-z0-9_]*"
    r")([^A-Za-z0-9]|$)",
    re.IGNORECASE
)

amr_annotation_pattern = re.compile(
    r"(aminoglycoside|hygromycin|bleomycin resistance|"
    r"macrolide|multidrug efflux|RND transporter|"
    r"efflux pump|FosA|phosphotransferase|"
    r"acetyltransferase-like protein)",
    re.IGNORECASE
)

mge_pattern = re.compile(
    r"(transposase|integrase|recombinase|resolvase|"
    r"insertion sequence|IS6 family|IS4 family|"
    r"IS15|IS26|IS1006|Tn3|TnAs|TnShfr|"
    r"prophage|phage|plasmid|conjug|"
    r"mobilization|relaxase)",
    re.IGNORECASE
)

replication_pattern = re.compile(
    r"(DNA replication protein|DnaC|"
    r"DNA-cytosine methyltransferase|methyltransferase)",
    re.IGNORECASE
)

hypothetical_pattern = re.compile(
    r"hypothetical protein",
    re.IGNORECASE
)

df["annotation_flag_AMR_or_efflux_v2"] = text.apply(
    lambda x: bool(
        amr_gene_pattern.search(x)
        or amr_annotation_pattern.search(x)
    )
)

df["annotation_flag_MGE_transposon_phage_v2"] = (
    text.apply(lambda x: bool(mge_pattern.search(x)))
)

df["annotation_flag_replication_methylation_v2"] = (
    text.apply(lambda x: bool(replication_pattern.search(x)))
)

df["annotation_flag_hypothetical_v2"] = (
    text.apply(lambda x: bool(hypothetical_pattern.search(x)))
)

amr = df["annotation_flag_AMR_or_efflux_v2"]
mge = df["annotation_flag_MGE_transposon_phage_v2"]
rep = df["annotation_flag_replication_methylation_v2"]
hyp = df["annotation_flag_hypothetical_v2"]

df["annotation_aware_priority_v2"] = "review"

df.loc[
    amr,
    "annotation_aware_priority_v2"
] = "possible_AMR_efflux_or_definition_related"

df.loc[
    (~amr) & mge,
    "annotation_aware_priority_v2"
] = "possible_MGE_transposon_phage_or_cargo_context"

df.loc[
    (~amr) & (~mge) & rep,
    "annotation_aware_priority_v2"
] = "possible_replication_or_methylation_context"

# High-priority non-obvious hypothetical candidates
df.loc[
    (~amr) & (~mge) & (~rep) & hyp &
    (df["freq_diff_pct"] >= 25) &
    (df["Hybrid_freq_pct"] >= 40),
    "annotation_aware_priority_v2"
] = "high_priority_hypothetical_nonobvious"

# High-effect low-frequency hypothetical candidates
df.loc[
    (~amr) & (~mge) & (~rep) & hyp &
    (df["freq_diff_pct"] >= 25) &
    (df["Hybrid_freq_pct"] < 40) &
    (df["NonHybrid_freq_pct"] <= 10),
    "annotation_aware_priority_v2"
] = "high_effect_low_frequency_hypothetical"

# Functional non-AMR candidates
df.loc[
    (~amr) & (~mge) & (~rep) & (~hyp),
    "annotation_aware_priority_v2"
] = "functional_non_AMR_candidate_review"

df.to_csv(outfile, sep="\t", index=False)

print("=== All category counts ===")
print(df["annotation_aware_priority_v2"].value_counts())

selected_categories = [
    "high_priority_hypothetical_nonobvious",
    "high_effect_low_frequency_hypothetical",
    "functional_non_AMR_candidate_review"
]

selected = df[
    df["annotation_aware_priority_v2"].isin(
        selected_categories
    )
].copy()

selected = selected.sort_values(
    ["annotation_aware_priority_v2", "q_value_BH_mainQC"]
)

selected.to_csv(
    "pyseer_results/Candidate15_Initial_Selection_Audit.tsv",
    sep="\t",
    index=False
)

selected["variant"].to_csv(
    "pyseer_results/Candidate15_Initial_Gene_IDs.txt",
    sep="\t",
    index=False,
    header=False
)

print("Prioritized candidate count:", len(selected))
print(selected[[
    "variant",
    "annotation_aware_priority_v2",
    "q_value_BH_mainQC",
    "odds_ratio",
    "freq_diff_pct"
]].to_string(index=False))
