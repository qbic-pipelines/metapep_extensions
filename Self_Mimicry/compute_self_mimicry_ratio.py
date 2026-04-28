#!/usr/bin/env python3
import sys
import pandas as pd
from pathlib import Path

# -------------------------------
# Usage
# -------------------------------
if len(sys.argv) < 2:
    print("Usage: python compute_self_mimicry_ratio.py <host_peptides.tsv>")
    print("Example: python compute_self_mimicry_ratio.py human_peptides.tsv.gz")
    sys.exit(1)

host_path = sys.argv[1]

# -------------------------------
# Load host peptides
# -------------------------------
print("Loading host peptides...")

host_df = pd.read_csv(host_path, sep="\t", compression="infer", dtype=str)
print(f"→ Loaded {len(host_df)} host proteins from {host_path}")

def parse_peptide_cell(s):
    if not isinstance(s, str):
        return []
    parts = [p.strip() for p in s.split(",")]
    return [p for p in parts if p]

peptide_cols = [c for c in host_df.columns if c.lower() == "peptides" or c.lower().startswith("peptides_") or c.lower() == "all_peptides"]

if not peptide_cols:
    sys.exit("No peptide columns found in host file (expected e.g. peptides, peptides_9mer, peptides_11mer, or all_peptides)")

host_peptides = set()
for col in peptide_cols:
    for val in host_df[col].fillna(""):
        host_peptides.update(parse_peptide_cell(val))

print(f"→ Unique host peptides found: {len(host_peptides):,}\n")

# -------------------------------
# Compute self-mimicry ratios per allele
# -------------------------------
print("Computing self-mimicry ratios per allele...")

def compute_ratio(peptide_string: str) -> float:
    if pd.isna(peptide_string) or not peptide_string.strip():
        return 0.0
    peptides = set(p.strip() for p in peptide_string.split(",") if p.strip())
    if not peptides:
        return 0.0
    shared = peptides.intersection(host_peptides)
    return len(shared) / len(peptides)

# Loop über alle entities_peptides.allele_X.tsv Dateien
for file in sorted(Path(".").glob("entities_peptides.allele_*.tsv")):
    allele_id = file.stem.split("_")[-1]
    print(f"→ Processing {file} (allele {allele_id})")

    entities_df = pd.read_csv(file, sep="\t", dtype=str)
    entities_df["self_mimicry_ratio"] = entities_df["peptides"].apply(compute_ratio)

    output_path = f"entities_self_mimicry_ratios.allele_{allele_id}.tsv"
    entities_df[["entity_id", "self_mimicry_ratio"]].to_csv(output_path, sep="\t", index=False)

    print(f"  Saved: {output_path} ({len(entities_df):,} entities)")

print("\nAll alleles processed successfully.")