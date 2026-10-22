#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


BLAST_FILE = Path(
    "06_TA_systems/03_blast/results/"
    "TADB2_vs_mcr139_proteins.raw.tsv"
)

QUERY_METADATA = Path(
    "06_TA_systems/02_TADB2_reference/processed/"
    "TADB2_typeII_component_metadata.tsv"
)

TARGET_METADATA = Path(
    "06_TA_systems/01_plasmid_proteins/"
    "mcr139_complete_plasmid_protein_map.tsv"
)

OUT_DIR = Path(
    "06_TA_systems/03_blast/results"
)

FILTERED_FILE = (
    OUT_DIR / "TADB2_component_hits.filtered.tsv"
)

BEST_FILE = (
    OUT_DIR / "TADB2_component_best_hit_by_target.tsv"
)

AMBIGUOUS_FILE = (
    OUT_DIR / "TADB2_component_ambiguous_targets.tsv"
)

CROSS_COMPONENT_FILE = (
    OUT_DIR / "TADB2_toxin_antitoxin_cross_hits.tsv"
)

PLASMID_SUMMARY_FILE = (
    OUT_DIR / "TADB2_component_hits_by_plasmid.tsv"
)

SUMMARY_FILE = (
    OUT_DIR / "TADB2_component_screening_summary.txt"
)


columns = [
    "qseqid",
    "sseqid",
    "pident",
    "alignment_length",
    "qlen",
    "qcovhsp",
    "qstart",
    "qend",
    "sstart",
    "send",
    "slen",
    "evalue",
    "bitscore",
]

if not BLAST_FILE.exists():
    raise SystemExit(
        f"ERROR: BLAST file does not exist: {BLAST_FILE}"
    )

if BLAST_FILE.stat().st_size == 0:
    raise SystemExit(
        "ERROR: BLAST output is empty."
    )

hits = pd.read_csv(
    BLAST_FILE,
    sep="\t",
    names=columns,
)

numeric_columns = [
    "pident",
    "alignment_length",
    "qlen",
    "qcovhsp",
    "qstart",
    "qend",
    "sstart",
    "send",
    "slen",
    "evalue",
    "bitscore",
]

for column in numeric_columns:
    hits[column] = pd.to_numeric(
        hits[column],
        errors="raise",
    )

# Strict thresholds: identity >80% and query coverage >80%, following the reference.
filtered = hits.loc[
    (hits["pident"] > 80)
    & (hits["qcovhsp"] > 80)
    & (hits["evalue"] < 1e-5)
].copy()

query_metadata = pd.read_csv(
    QUERY_METADATA,
    sep="\t",
    dtype={
        "query_id": str,
        "pair_id": int,
        "component": str,
    },
)

target_metadata = pd.read_csv(
    TARGET_METADATA,
    sep="\t",
    dtype={
        "target_protein_id": str,
        "accession": str,
    },
)

filtered = filtered.merge(
    query_metadata,
    left_on="qseqid",
    right_on="query_id",
    how="left",
    validate="many_to_one",
)

if filtered["component"].isna().any():
    missing = filtered.loc[
        filtered["component"].isna(),
        "qseqid",
    ].drop_duplicates()

    raise SystemExit(
        "ERROR: query metadata missing for: "
        + ",".join(missing.head(20))
    )

filtered = filtered.merge(
    target_metadata,
    left_on="sseqid",
    right_on="target_protein_id",
    how="left",
    validate="many_to_one",
)

if filtered["accession"].isna().any():
    missing = filtered.loc[
        filtered["accession"].isna(),
        "sseqid",
    ].drop_duplicates()

    raise SystemExit(
        "ERROR: target metadata missing for: "
        + ",".join(missing.head(20))
    )

filtered = filtered.sort_values(
    [
        "accession",
        "sseqid",
        "component",
        "bitscore",
        "pident",
        "qcovhsp",
        "alignment_length",
        "evalue",
        "qseqid",
    ],
    ascending=[
        True,
        True,
        True,
        False,
        False,
        False,
        False,
        True,
        True,
    ],
)

filtered.to_csv(
    FILTERED_FILE,
    sep="\t",
    index=False,
)

# Keep the best reference hit per target protein and component.
best = filtered.drop_duplicates(
    subset=[
        "sseqid",
        "component",
    ],
    keep="first",
).copy()

best.to_csv(
    BEST_FILE,
    sep="\t",
    index=False,
)

# Target protein and component combinations matching multiple TA pair IDs.
ambiguity = (
    filtered
    .groupby(
        [
            "accession",
            "sseqid",
            "component",
        ],
        as_index=False,
    )
    .agg(
        distinct_pair_ids=(
            "pair_id",
            "nunique",
        ),
        pair_ids=(
            "pair_id",
            lambda values: ",".join(
                map(
                    str,
                    sorted(set(values)),
                )
            ),
        ),
        passing_reference_hits=(
            "qseqid",
            "nunique",
        ),
        maximum_bitscore=(
            "bitscore",
            "max",
        ),
    )
)

ambiguity = ambiguity.loc[
    ambiguity["distinct_pair_ids"] > 1
].copy()

ambiguity.to_csv(
    AMBIGUOUS_FILE,
    sep="\t",
    index=False,
)

# Target proteins receiving both toxin and antitoxin best hits.
cross_component = (
    best
    .groupby(
        [
            "accession",
            "sseqid",
        ],
        as_index=False,
    )
    .agg(
        component_count=(
            "component",
            "nunique",
        ),
        components=(
            "component",
            lambda values: ",".join(
                sorted(set(values))
            ),
        ),
        pair_ids=(
            "pair_id",
            lambda values: ",".join(
                map(
                    str,
                    sorted(set(values)),
                )
            ),
        ),
    )
)

cross_component = cross_component.loc[
    cross_component["component_count"] > 1
].copy()

cross_component.to_csv(
    CROSS_COMPONENT_FILE,
    sep="\t",
    index=False,
)

plasmid_summary = (
    best
    .groupby(
        [
            "accession",
            "component",
        ],
        as_index=False,
    )
    .agg(
        best_hit_proteins=(
            "sseqid",
            "nunique",
        ),
        distinct_pair_ids=(
            "pair_id",
            "nunique",
        ),
    )
)

plasmid_summary.to_csv(
    PLASMID_SUMMARY_FILE,
    sep="\t",
    index=False,
)

best_toxin = best.loc[
    best["component"] == "toxin"
]

best_antitoxin = best.loc[
    best["component"] == "antitoxin"
]

toxin_plasmids = best_toxin[
    "accession"
].nunique()

antitoxin_plasmids = best_antitoxin[
    "accession"
].nunique()

all_positive_plasmids = best[
    "accession"
].nunique()

summary_lines = [
    "===== TADB2 TYPE II COMPONENT SCREENING =====",
    f"Raw BLAST HSPs: {len(hits)}",
    (
        "Threshold-passing HSPs "
        "(identity >80%, query coverage >80%, "
        f"E-value <1e-5): {len(filtered)}"
    ),
    (
        "Unique TADB reference queries with passing hits: "
        f"{filtered['qseqid'].nunique()}"
    ),
    (
        "Unique target proteins with passing hits: "
        f"{filtered['sseqid'].nunique()}"
    ),
    (
        "Best toxin-hit target proteins: "
        f"{len(best_toxin)}"
    ),
    (
        "Best antitoxin-hit target proteins: "
        f"{len(best_antitoxin)}"
    ),
    (
        "Plasmids with at least one toxin homolog: "
        f"{toxin_plasmids}/937"
    ),
    (
        "Plasmids with at least one antitoxin homolog: "
        f"{antitoxin_plasmids}/937"
    ),
    (
        "Plasmids with any Type II TA component homolog: "
        f"{all_positive_plasmids}/937"
    ),
    (
        "Target/component combinations matching "
        "multiple TADB pair IDs: "
        f"{len(ambiguity)}"
    ),
    (
        "Target proteins receiving both toxin and "
        "antitoxin classifications: "
        f"{len(cross_component)}"
    ),
    "",
    (
        "These are component-level homologs, "
        "not validated complete TA systems."
    ),
    "",
    "STATUS: PASSED",
]

SUMMARY_FILE.write_text(
    "\n".join(summary_lines) + "\n",
    encoding="utf-8",
)

print("\n".join(summary_lines))
