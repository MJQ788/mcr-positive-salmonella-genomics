#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path

def categorize(raw):
    text = str(raw).strip().lower()
    if "human" in text or text in ("homo sapiens", "clinical"):
        return "human"
    if any(x in text for x in ("chicken", "egg", "poultry", "broiler")):
        return "chicken_related"
    return "other_or_unknown"

def load_metadata(path, id_col, source_col):
    if path is None:
        return {}
    path = Path(path)
    dialect = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter=dialect)
        if id_col not in reader.fieldnames or source_col not in reader.fieldnames:
            raise SystemExit(
                f"Metadata requires columns '{id_col}' and '{source_col}'. "
                f"Found: {reader.fieldnames}"
            )
        result = {}
        for row in reader:
            isolate = row[id_col].strip()
            if isolate in result:
                raise SystemExit(f"Duplicate metadata ID: {isolate}")
            result[isolate] = categorize(row[source_col])
        return result

class DSU:
    def __init__(self, ids):
        self.parent = {i: i for i in ids}
        self.size = {i: 1 for i in ids}
    def find(self, i):
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i
    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return
        if self.size[a] < self.size[b]:
            a, b = b, a
        self.parent[b] = a
        self.size[a] += self.size[b]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--distances", required=True,
                    help="Tab-delimited square SNP-dists output")
    ap.add_argument("--metadata", help="Optional CSV/TSV with isolate source")
    ap.add_argument("--id-column", default="Assembly")
    ap.add_argument("--source-column", default="Isolation source")
    ap.add_argument("--threshold", type=int, default=10)
    ap.add_argument("--exclude-id", action="append", default=["Reference"],
                    help="Sequence ID to exclude; default: Reference")
    ap.add_argument("--outdir", default="snp_network")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    metadata = load_metadata(args.metadata, args.id_column, args.source_column)
    excluded = set(args.exclude_id)

    with open(args.distances, newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        if len(header) < 2:
            raise SystemExit("Distance file has no sample column header")
        col_ids = [x.strip() for x in header[1:]]
        data = {}
        for line in reader:
            if not line or all(not c for c in line):
                continue
            name = line[0].strip()
            if len(line) != len(col_ids) + 1:
                raise SystemExit(f"Malformed distance row for {name}")
            if name in data:
                raise SystemExit(f"Duplicate SNP-distance row: {name}")
            data[name] = line[1:]

    if set(data) != set(col_ids):
        raise SystemExit("SNP-dists row IDs and column IDs differ")
    ids = sorted(set(data) - excluded)
    if args.metadata:
        missing = set(ids) - set(metadata)
        if missing:
            raise SystemExit(
                f"{len(missing)} sample IDs missing in metadata (example: {sorted(missing)[:5]})"
            )
    parents = DSU(ids)
    indexes = {name: i for i, name in enumerate(col_ids)}
    edges = []

    for index, a in enumerate(ids):
        row = data[a]
        for b in ids[index + 1:]:
            d = int(float(row[indexes[b]]))
            if d <= args.threshold:
                parents.union(a, b)
                edges.append((a, b, d))

    components = defaultdict(list)
    for node in ids:
        components[parents.find(node)].append(node)
    sorted_components = sorted(
        (sorted(nodes) for nodes in components.values()),
        key=lambda nodes: (nodes[0], len(nodes))
    )
    component_for_node = {
        node: f"C{idx:04d}"
        for idx, nodes in enumerate(sorted_components, start=1)
        for node in nodes
    }

    with (out / "edges_leq10.tsv").open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_a", "sample_b", "core_snp_distance", "component_id"])
        writer.writerows(
            [a, b, d, component_for_node[a]] for a, b, d in edges
        )

    mixed_count = 0
    mixed_n = 0
    with (out / "components.tsv").open("w", newline="") as f, \
         (out / "mixed_source_nodes.tsv").open("w", newline="") as g:
        summary = csv.writer(f, delimiter="\t")
        mixed = csv.writer(g, delimiter="\t")
        summary.writerow([
            "component_id", "n_genomes", "n_human", "n_chicken_related",
            "n_other_or_unknown", "mixed_source"
        ])
        mixed.writerow(["component_id", "sample_id", "source_group"])
        for nodes in sorted_components:
            groups = [metadata.get(x, "other_or_unknown") for x in nodes]
            n_h = groups.count("human")
            n_c = groups.count("chicken_related")
            is_mixed = n_h > 0 and n_c > 0
            cid = component_for_node[nodes[0]]
            summary.writerow([
                cid, len(nodes), n_h, n_c,
                len(nodes) - n_h - n_c, int(is_mixed)
            ])
            if is_mixed:
                mixed_count += 1
                mixed_n += len(nodes)
                for node in nodes:
                    mixed.writerow([cid, node, metadata.get(node, "other_or_unknown")])

    print(f"Genomes analyzed: {len(ids)}")
    print(f"Edges with distance <= {args.threshold}: {len(edges)}")
    print(f"All connected components: {len(sorted_components)}")
    if args.metadata:
        print(f"Mixed-source components: {mixed_count}; genomes in these: {mixed_n}")
    else:
        print("No metadata supplied: mixed-source status is not evaluated.")
    print(f"Outputs: {out}")

if __name__ == "__main__":
    main()
