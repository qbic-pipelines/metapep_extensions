#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import sys
import gzip

# -------------------------------
# Arguments
# -------------------------------
if len(sys.argv) < 5:
    print("Usage: python prepare_self_mimicry_ratios.py <entities_proteins.tsv> <proteins_peptides.tsv> <peptides.tsv.gz> <predictions.tsv.gz>")
    sys.exit(1)

entities_proteins_path = sys.argv[1]
proteins_peptides_path = sys.argv[2]
peptides_path = sys.argv[3]
predictions_path = sys.argv[4]

bindingThreshold = 0.426
mapEnitityPeptideIDs = False
mergeEntities = False
mapPeptideIDsSequences = True


# -------------------------------
# Prepare output
# -------------------------------
output_ids_path = Path("entities_peptides_ids.tsv")
output_final_path = Path("entities_peptides.tsv")

# -------------------------------
# Step 1: Chunkwise entity → peptide_id mapping and filter peptide IDs by binding score
# -------------------------------
if mapEnitityPeptideIDs:
    print("Filtering peptides by prediction score ≥", bindingThreshold, "...")

    pred_reader = pd.read_csv(
        predictions_path,
        sep="\t",
        dtype={"peptide_id": str, "prediction_score": float, "allele_id": str},
        usecols=["peptide_id", "prediction_score", "allele_id"],
        compression="gzip",
        chunksize=5_000_000
    )

    passing_peptides = {}
    for i, chunk in enumerate(pred_reader, start=1):
        chunk_filtered = chunk.loc[chunk["prediction_score"] >= bindingThreshold]
        print(f"  → Chunk {i}: {len(chunk_filtered):,} passing rows")

        for allele_id, group in chunk_filtered.groupby("allele_id"):
            if allele_id not in passing_peptides:
                passing_peptides[allele_id] = set()
            passing_peptides[allele_id].update(group["peptide_id"].astype(str))

    print("\nSummary of passing peptides per allele:")
    for allele_id, peps in passing_peptides.items():
        print(f"  Allele {allele_id}: {len(peps):,} peptides passed threshold")

    chunk_size = 2_000_000
    entities_proteins = pd.read_csv(entities_proteins_path, sep="\t", dtype=str)

    for allele_id, peptide_set in passing_peptides.items():
        print(f"\n→ Processing allele {allele_id} ({len(peptide_set):,} peptides)")
        output_ids_path = Path(f"entities_peptides_ids.allele_{allele_id}.tsv")
        if output_ids_path.exists():
            output_ids_path.unlink()

        reader = pd.read_csv(proteins_peptides_path, sep="\t", dtype=str, usecols=["protein_id", "peptide_id"], chunksize=chunk_size)

        for i, chunk in enumerate(reader, start=1):
            print(f"  Chunk {i}: {len(chunk):,} protein→peptide rows")
            chunk = chunk[chunk["peptide_id"].isin(peptide_set)]
            print(f"    → {len(chunk):,} rows remain after filtering")

            merged = chunk.merge(entities_proteins, on="protein_id", how="left")
            merged.dropna(subset=["entity_id"], inplace=True)

            grouped = merged.groupby("entity_id")["peptide_id"].apply(lambda x: ",".join(sorted(set(x)))).reset_index()
            grouped.to_csv(output_ids_path, sep="\t", index=False, mode="a", header=(i == 1))

        print(f"  Finished allele {allele_id} mapping. Output: {output_ids_path}")
    
    print("\nStep 1 finished.\n")

# -------------------------------
# Step 2: Group duplicate entities (direkt nach Step 1)
# -------------------------------
if mergeEntities:
    print("\nStep 2: Group duplicate entities and merge peptide IDs...")

    for file in Path(".").glob("entities_peptides_ids.allele_*.tsv"):
        allele_id = file.stem.split("_")[-1]
        print(f"  Processing {file} (allele {allele_id})")
        entities_df = pd.read_csv(file, sep="\t", dtype=str)
        entities_df = entities_df.groupby("entity_id")["peptide_id"].apply(lambda lists: ",".join(sorted(set(",".join(lists).split(","))))).reset_index() # Group and joined data
        entities_df.rename(columns={"peptide_id": "peptide_ids"}, inplace=True)
        entities_df.to_csv(file, sep="\t", index=False)
        print(f"  Finished merging for allele {allele_id}")
    
    print("Step 2 finished. Duplicate entities merged.")

# -------------------------------
# Step 3: Map peptide IDs → sequences (chunked peptide DB, inplace)
# -------------------------------
if mapPeptideIDsSequences:
    print("Step 3: Map peptide IDs → sequences (per allele)")

    files = sorted(Path(".").glob("entities_peptides_ids.allele_*.tsv"), key=lambda p: int(p.stem.split("_")[-1]))
    for file in files:
        allele_id = file.stem.split("_")[-1]
        output_final_path = Path(f"entities_peptides.allele_{allele_id}.tsv")
        print(f"  → Mapping peptides for allele {allele_id}")

        entities_df = pd.read_csv(file, sep="\t", dtype=str)

        pep_reader = pd.read_csv(
            peptides_path,
            sep="\t",
            dtype=str,
            usecols=["peptide_id", "peptide_sequence"],
            compression="gzip",
            chunksize=2_000_000,
            engine="python"
        )

        for i, pep_chunk in enumerate(pep_reader, start=1):
            pep_dict = pep_chunk.set_index("peptide_id")["peptide_sequence"].to_dict()

            def map_ids(ids_str):
                return ",".join([pep_dict.get(pid, pid) for pid in ids_str.split(",")])

            entities_df["peptide_ids"] = entities_df["peptide_ids"].apply(map_ids)
            print(f"    Chunk {i:03d}: mapped {len(pep_dict):,} peptides")
            del pep_dict

        entities_df.rename(columns={"peptide_ids": "peptides"}, inplace=True)
        entities_df.to_csv(output_final_path, sep="\t", index=False)
        print(f"  Finished allele {allele_id} → saved {output_final_path}")

    print("\nStep 3 finished.")
