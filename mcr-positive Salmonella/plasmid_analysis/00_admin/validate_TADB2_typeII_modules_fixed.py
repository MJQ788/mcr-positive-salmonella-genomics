#!/usr/bin/env python3

from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
import math

import pandas as pd


HIT_FILE = Path(
    "06_TA_systems/03_blast/results/"
    "TADB2_component_hits.filtered.tsv"
)

PROTEIN_MAP_FILE = Path(
    "06_TA_systems/01_plasmid_proteins/"
    "mcr139_complete_plasmid_protein_map.tsv"
)

PLASMID_SUMMARY_FILE = Path(
    "06_TA_systems/01_plasmid_proteins/"
    "protein_count_by_plasmid.tsv"
)

OUT_DIR = Path(
    "06_TA_systems/04_pair_validation"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)

ALL_CANDIDATES_FILE = (
    OUT_DIR / "TADB2_physical_pair_candidates.tsv"
)

STRICT_CANDIDATES_FILE = (
    OUT_DIR / "TADB2_strict_adjacent_pair_candidates.tsv"
)

SELECTED_MODULES_FILE = (
    OUT_DIR / "TADB2_selected_nonoverlapping_modules.tsv"
)

PAIR_SUMMARY_FILE = (
    OUT_DIR / "TADB2_selected_pair_id_summary.tsv"
)

PLASMID_MODULE_SUMMARY_FILE = (
    OUT_DIR / "TADB2_module_count_by_plasmid.tsv"
)

REUSED_PROTEIN_FILE = (
    OUT_DIR / "TADB2_proteins_in_multiple_candidates.tsv"
)

SUMMARY_FILE = (
    OUT_DIR / "TADB2_module_validation_summary.txt"
)


MAX_INTERGENIC_GAP_BP = 500


def first_existing(columns, candidates):
    for candidate in candidates:
        if candidate in columns:
            return candidate
    return None


hits = pd.read_csv(
    HIT_FILE,
    sep="\t",
    low_memory=False,
)

proteins = pd.read_csv(
    PROTEIN_MAP_FILE,
    sep="\t",
    dtype={
        "target_protein_id": str,
        "accession": str,
    },
)

plasmids = pd.read_csv(
    PLASMID_SUMMARY_FILE,
    sep="\t",
    dtype={"accession": str},
)

required_hit_columns = [
    "sseqid",
    "pair_id",
    "component",
    "pident",
    "qcovhsp",
    "evalue",
    "bitscore",
]

missing = [
    column
    for column in required_hit_columns
    if column not in hits.columns
]

if missing:
    raise SystemExit(
        "ERROR: missing hit columns: "
        + ", ".join(missing)
    )

required_protein_columns = [
    "target_protein_id",
    "accession",
    "protein_index",
    "start",
    "end",
    "strand",
    "gene",
    "product",
]

missing = [
    column
    for column in required_protein_columns
    if column not in proteins.columns
]

if missing:
    raise SystemExit(
        "ERROR: missing protein-map columns: "
        + ", ".join(missing)
    )

for column in [
    "pair_id",
    "pident",
    "qcovhsp",
    "evalue",
    "bitscore",
]:
    hits[column] = pd.to_numeric(
        hits[column],
        errors="raise",
    )

for column in [
    "protein_index",
    "start",
    "end",
]:
    proteins[column] = pd.to_numeric(
        proteins[column],
        errors="raise",
    )

plasmids["translated_proteins"] = pd.to_numeric(
    plasmids["translated_proteins"],
    errors="raise",
)

plasmids["plasmid_length_bp"] = pd.to_numeric(
    plasmids["plasmid_length_bp"],
    errors="raise",
)

description_column = first_existing(
    hits.columns,
    [
        "original_description_aliases",
        "original_description",
        "original_id_aliases",
        "original_id",
    ],
)

if description_column is None:
    hits["reference_description"] = ""
    description_column = "reference_description"



best_pair_hits = (
    hits
    .sort_values(
        [
            "sseqid",
            "component",
            "pair_id",
            "bitscore",
            "pident",
            "qcovhsp",
            "evalue",
        ],
        ascending=[
            True,
            True,
            True,
            False,
            False,
            False,
            True,
        ],
    )
    .drop_duplicates(
        subset=[
            "sseqid",
            "component",
            "pair_id",
        ],
        keep="first",
    )
    .copy()
)

hit_lookup = defaultdict(dict)

for row in best_pair_hits.itertuples(index=False):
    key = (
        str(row.sseqid),
        str(row.component),
    )
    hit_lookup[key][int(row.pair_id)] = row


positive_protein_ids = set(
    best_pair_hits["sseqid"].astype(str)
)

positive_proteins = proteins.loc[
    proteins["target_protein_id"].isin(
        positive_protein_ids
    )
].copy()

protein_metadata = (
    proteins
    .set_index("target_protein_id")
    .to_dict("index")
)

plasmid_metadata = (
    plasmids
    .set_index("accession")
    .to_dict("index")
)


def linear_gap(row_a, row_b):
    if row_a["end"] < row_b["start"]:
        return int(
            row_b["start"]
            - row_a["end"]
            - 1
        )

    if row_b["end"] < row_a["start"]:
        return int(
            row_a["start"]
            - row_b["end"]
            - 1
        )

    return 0


def circular_gap(row_a, row_b, plasmid_length):
    leftmost = min(
        row_a["start"],
        row_b["start"],
    )

    rightmost = max(
        row_a["end"],
        row_b["end"],
    )

    return int(
        plasmid_length
        - rightmost
        + leftmost
        - 1
    )


def proteins_between(
    index_a,
    index_b,
    total_proteins,
    use_circular_gap,
):
    difference = abs(
        int(index_a) - int(index_b)
    )

    if use_circular_gap:
        return max(
            0,
            int(total_proteins)
            - difference
            - 1,
        )

    return max(0, difference - 1)


def component_pairs(
    protein_a,
    protein_b,
):
    possibilities = []

    for role_a, role_b in [
        ("toxin", "antitoxin"),
        ("antitoxin", "toxin"),
    ]:
        hits_a = hit_lookup.get(
            (protein_a, role_a),
            {},
        )

        hits_b = hit_lookup.get(
            (protein_b, role_b),
            {},
        )

        shared_ids = sorted(
            set(hits_a) & set(hits_b)
        )

        if not shared_ids:
            continue

        scored_ids = []

        for pair_id in shared_ids:
            hit_a = hits_a[pair_id]
            hit_b = hits_b[pair_id]

            scored_ids.append({
                "pair_id": pair_id,
                "combined_bitscore": (
                    float(hit_a.bitscore)
                    + float(hit_b.bitscore)
                ),
                "minimum_identity": min(
                    float(hit_a.pident),
                    float(hit_b.pident),
                ),
                "minimum_query_coverage": min(
                    float(hit_a.qcovhsp),
                    float(hit_b.qcovhsp),
                ),
                "maximum_evalue": max(
                    float(hit_a.evalue),
                    float(hit_b.evalue),
                ),
                "hit_a": hit_a,
                "hit_b": hit_b,
            })

        scored_ids.sort(
            key=lambda item: (
                -item["combined_bitscore"],
                -item["minimum_identity"],
                -item["minimum_query_coverage"],
                item["maximum_evalue"],
                item["pair_id"],
            )
        )

        possibilities.append({
            "role_a": role_a,
            "role_b": role_b,
            "shared_pair_ids": shared_ids,
            "best": scored_ids[0],
        })

    possibilities.sort(
        key=lambda item: (
            -item["best"]["combined_bitscore"],
            -item["best"]["minimum_identity"],
            -item["best"]["minimum_query_coverage"],
            item["best"]["maximum_evalue"],
        )
    )

    return possibilities


candidate_rows = []


for accession, group in positive_proteins.groupby(
    "accession",
    sort=True,
):
    group = group.sort_values(
        ["start", "end", "protein_index"]
    )

    plasmid_length = int(
        plasmid_metadata[accession][
            "plasmid_length_bp"
        ]
    )

    total_proteins = int(
        plasmid_metadata[accession][
            "translated_proteins"
        ]
    )

    records = group.to_dict("records")

    for row_a, row_b in combinations(
        records,
        2,
    ):
        protein_a = str(
            row_a["target_protein_id"]
        )

        protein_b = str(
            row_b["target_protein_id"]
        )

        possibilities = component_pairs(
            protein_a,
            protein_b,
        )

        if not possibilities:
            continue

        linear_distance = linear_gap(
            row_a,
            row_b,
        )

        wrap_distance = circular_gap(
            row_a,
            row_b,
            plasmid_length,
        )

        use_circular = (
            wrap_distance < linear_distance
        )

        genomic_gap = min(
            linear_distance,
            wrap_distance,
        )

        between = proteins_between(
            row_a["protein_index"],
            row_b["protein_index"],
            total_proteins,
            use_circular,
        )

        best_orientation = possibilities[0]
        best = best_orientation["best"]

        hit_a = best["hit_a"]
        hit_b = best["hit_b"]

        same_strand = (
            str(row_a["strand"])
            == str(row_b["strand"])
        )

        strict_adjacent = (
            genomic_gap
            <= MAX_INTERGENIC_GAP_BP
            and between == 0
        )

        candidate_rows.append({
            "accession": accession,
            "protein_a": protein_a,
            "protein_b": protein_b,
            "role_a": (
                best_orientation["role_a"]
            ),
            "role_b": (
                best_orientation["role_b"]
            ),
            "protein_a_index": int(
                row_a["protein_index"]
            ),
            "protein_b_index": int(
                row_b["protein_index"]
            ),
            "protein_a_start": int(
                row_a["start"]
            ),
            "protein_a_end": int(
                row_a["end"]
            ),
            "protein_b_start": int(
                row_b["start"]
            ),
            "protein_b_end": int(
                row_b["end"]
            ),
            "protein_a_strand": row_a["strand"],
            "protein_b_strand": row_b["strand"],
            "protein_a_gene": row_a["gene"],
            "protein_b_gene": row_b["gene"],
            "protein_a_product": row_a["product"],
            "protein_b_product": row_b["product"],
            "best_pair_id": int(
                best["pair_id"]
            ),
            "compatible_pair_id_count": len(
                best_orientation[
                    "shared_pair_ids"
                ]
            ),
            "compatible_pair_ids": ",".join(
                map(
                    str,
                    best_orientation[
                        "shared_pair_ids"
                    ],
                )
            ),
            "orientation_possibility_count": len(
                possibilities
            ),
            "combined_bitscore": (
                best["combined_bitscore"]
            ),
            "minimum_identity_percent": (
                best["minimum_identity"]
            ),
            "minimum_query_coverage_percent": (
                best[
                    "minimum_query_coverage"
                ]
            ),
            "maximum_evalue": (
                best["maximum_evalue"]
            ),
            "protein_a_reference_query": (
                hit_a.qseqid
            ),
            "protein_b_reference_query": (
                hit_b.qseqid
            ),
            "protein_a_reference_description": (
                getattr(
                    hit_a,
                    description_column,
                    "",
                )
            ),
            "protein_b_reference_description": (
                getattr(
                    hit_b,
                    description_column,
                    "",
                )
            ),
            "linear_gap_bp": linear_distance,
            "circular_gap_bp": wrap_distance,
            "genomic_gap_bp": genomic_gap,
            "circular_origin_pair": use_circular,
            "translated_proteins_between": (
                between
            ),
            "same_strand": same_strand,
            "strict_adjacent_module": (
                strict_adjacent
            ),
        })


candidates = pd.DataFrame(candidate_rows)

if candidates.empty:
    raise SystemExit(
        "ERROR: no physically compatible "
        "TA component pairs were found."
    )

candidates = candidates.sort_values(
    [
        "accession",
        "strict_adjacent_module",
        "combined_bitscore",
        "minimum_identity_percent",
        "minimum_query_coverage_percent",
        "genomic_gap_bp",
        "protein_a",
        "protein_b",
    ],
    ascending=[
        True,
        False,
        False,
        False,
        False,
        True,
        True,
        True,
    ],
)

candidates.to_csv(
    ALL_CANDIDATES_FILE,
    sep="\t",
    index=False,
)

strict = candidates.loc[
    candidates["strict_adjacent_module"]
].copy()

strict.to_csv(
    STRICT_CANDIDATES_FILE,
    sep="\t",
    index=False,
)



selected_rows = []
used_proteins_by_plasmid = defaultdict(set)

strict_for_selection = strict.sort_values(
    [
        "accession",
        "combined_bitscore",
        "minimum_identity_percent",
        "minimum_query_coverage_percent",
        "genomic_gap_bp",
        "compatible_pair_id_count",
    ],
    ascending=[
        True,
        False,
        False,
        False,
        True,
        True,
    ],
)

for row in strict_for_selection.to_dict(
    "records"
):
    accession = row["accession"]

    protein_a = row["protein_a"]
    protein_b = row["protein_b"]

    used = used_proteins_by_plasmid[
        accession
    ]

    if protein_a in used or protein_b in used:
        continue

    selected_rows.append(row)
    used.add(protein_a)
    used.add(protein_b)


selected = pd.DataFrame(selected_rows)

if selected.empty:
    selected = strict.iloc[0:0].copy()

selected.to_csv(
    SELECTED_MODULES_FILE,
    sep="\t",
    index=False,
)


# Proteins participating in more than one strict candidate module.
protein_candidate_counts = Counter()

for row in strict.itertuples(index=False):
    protein_candidate_counts[row.protein_a] += 1
    protein_candidate_counts[row.protein_b] += 1

reused_rows = []

for protein_id, count in sorted(
    protein_candidate_counts.items()
):
    if count <= 1:
        continue

    metadata = protein_metadata[protein_id]

    reused_rows.append({
        "target_protein_id": protein_id,
        "accession": metadata["accession"],
        "strict_candidate_module_count": count,
        "gene": metadata["gene"],
        "product": metadata["product"],
    })

pd.DataFrame(
    reused_rows,
    columns=[
        "target_protein_id",
        "accession",
        "strict_candidate_module_count",
        "gene",
        "product",
    ],
).to_csv(
    REUSED_PROTEIN_FILE,
    sep="\t",
    index=False,
)


if not selected.empty:
    pair_summary = (
        selected
        .groupby(
            ["best_pair_id"],
            as_index=False,
        )
        .agg(
            selected_modules=(
                "accession",
                "size",
            ),
            positive_plasmids=(
                "accession",
                "nunique",
            ),
            median_combined_bitscore=(
                "combined_bitscore",
                "median",
            ),
            median_identity_percent=(
                "minimum_identity_percent",
                "median",
            ),
            median_query_coverage_percent=(
                "minimum_query_coverage_percent",
                "median",
            ),
            median_gap_bp=(
                "genomic_gap_bp",
                "median",
            ),
            same_strand_modules=(
                "same_strand",
                "sum",
            ),
            example_toxin_description=(
                "protein_a_reference_description",
                "first",
            ),
            example_antitoxin_description=(
                "protein_b_reference_description",
                "first",
            ),
        )
        .sort_values(
            [
                "positive_plasmids",
                "selected_modules",
                "best_pair_id",
            ],
            ascending=[
                False,
                False,
                True,
            ],
        )
    )
else:
    pair_summary = pd.DataFrame()

pair_summary.to_csv(
    PAIR_SUMMARY_FILE,
    sep="\t",
    index=False,
)


if not selected.empty:
    plasmid_module_summary = (
        selected
        .groupby(
            "accession",
            as_index=False,
        )
        .agg(
            selected_TypeII_modules=(
                "best_pair_id",
                "size",
            ),
            distinct_TADB_pair_ids=(
                "best_pair_id",
                "nunique",
            ),
            pair_ids=(
                "best_pair_id",
                lambda values: ",".join(
                    map(
                        str,
                        sorted(set(values)),
                    )
                ),
            ),
        )
    )
else:
    plasmid_module_summary = pd.DataFrame(
        columns=[
            "accession",
            "selected_TypeII_modules",
            "distinct_TADB_pair_ids",
            "pair_ids",
        ]
    )

all_plasmids = plasmids[
    ["accession"]
].drop_duplicates()

plasmid_module_summary = (
    all_plasmids
    .merge(
        plasmid_module_summary,
        on="accession",
        how="left",
        validate="one_to_one",
    )
)

for column in [
    "selected_TypeII_modules",
    "distinct_TADB_pair_ids",
]:
    plasmid_module_summary[column] = (
        plasmid_module_summary[column]
        .fillna(0)
        .astype(int)
    )

plasmid_module_summary[
    "pair_ids"
] = (
    plasmid_module_summary["pair_ids"]
    .fillna("")
)

plasmid_module_summary.to_csv(
    PLASMID_MODULE_SUMMARY_FILE,
    sep="\t",
    index=False,
)


selected_plasmids = (
    selected["accession"].nunique()
    if not selected.empty
    else 0
)

multi_module_plasmids = int(
    (
        plasmid_module_summary[
            "selected_TypeII_modules"
        ] >= 2
    ).sum()
)

strict_same_strand = int(
    strict["same_strand"].sum()
)

selected_same_strand = (
    int(selected["same_strand"].sum())
    if not selected.empty
    else 0
)

summary_lines = [
    "===== TADB2 TYPE II MODULE VALIDATION =====",
    (
        "TA-component-positive target proteins: "
        f"{len(positive_protein_ids)}"
    ),
    (
        "Protein pairs sharing at least one TADB "
        f"pair ID: {len(candidates)}"
    ),
    (
        "Strict adjacent candidates "
        f"(gap <= {MAX_INTERGENIC_GAP_BP} bp "
        "and no translated CDS between): "
        f"{len(strict)}"
    ),
    (
        "Strict candidates on the same strand: "
        f"{strict_same_strand}/{len(strict)}"
    ),
    (
        "Selected non-overlapping Type II "
        f"candidate modules: {len(selected)}"
    ),
    (
        "Selected modules on the same strand: "
        f"{selected_same_strand}/{len(selected)}"
    ),
    (
        "Plasmids with at least one selected "
        f"candidate Type II module: "
        f"{selected_plasmids}/937"
    ),
    (
        "Plasmids with two or more selected "
        f"candidate Type II modules: "
        f"{multi_module_plasmids}/937"
    ),
    (
        "Distinct selected TADB pair IDs: "
        f"{selected['best_pair_id'].nunique() if not selected.empty else 0}"
    ),
    (
        "Proteins participating in multiple "
        "strict candidate modules: "
        f"{len(reused_rows)}"
    ),
    "",
    (
        "Primary module definition: two distinct "
        "proteins on the same complete plasmid, "
        "with no translated CDS between them, "
        f"an intergenic/overlap distance <= "
        f"{MAX_INTERGENIC_GAP_BP} bp, and at least "
        "one shared numeric TADB toxin-antitoxin "
        "reference pair identifier."
    ),
    "",
    (
        "Selected modules were greedily resolved "
        "by combined bitscore so that one protein "
        "was not counted in more than one module."
    ),
    "",
    (
        "These remain sequence-supported candidate "
        "Type II maintenance modules"
    ),
    "",
    "STATUS: PASSED",
]

SUMMARY_FILE.write_text(
    "\n".join(summary_lines) + "\n",
    encoding="utf-8",
)

print("\n".join(summary_lines))
