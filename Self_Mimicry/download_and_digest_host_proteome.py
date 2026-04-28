#!/usr/bin/env python3
"""
Download a complete proteome from NCBI and digest it into peptides.
"""

import os
import sys
import time
import gzip
import pandas as pd
from Bio import Entrez, SeqIO
from urllib.error import HTTPError
import argparse


# -------------------------------
# Helper: peptide digestion
# -------------------------------
def generate_peptides(sequence, lengths):
    peptides_by_length = {L: [] for L in lengths}
    seq_len = len(sequence)
    for L in lengths:
        if seq_len < L:
            continue
        peptides_by_length[L] = [sequence[i:i+L] for i in range(seq_len - L + 1)]
    return peptides_by_length


# -------------------------------
# CLI
# -------------------------------
parser = argparse.ArgumentParser(description="Download and digest proteome via NCBI Entrez.")
parser.add_argument("-t", "--taxon", default="9606" required=True, help="NCBI TaxID")
parser.add_argument("-e", "--email", required=True, help="Email for NCBI Entrez")
parser.add_argument("-k", "--key", help="NCBI API key (optional)")
parser.add_argument("-l", "--lengths", default="9", help="Comma-separated peptide lengths (e.g. '9,11,12')")
parser.add_argument("-o", "--output", default="proteome_peptides.tsv", help="Output TSV file")
args = parser.parse_args()

lengths = [int(x) for x in args.lengths.split(",")]
Entrez.email = args.email
if args.key:
    Entrez.api_key = args.key

# -------------------------------
# Step 1: Fetch proteome
# -------------------------------
taxid = args.taxon  # TODO: Maybe also allow Species names and resolve TaxID
query = f"txid{taxid}[Organism:exp] AND refseq[filter]"
print(f"Downloading proteome for {args.taxon}...", file=sys.stderr)

retmax = 10000
retstart = 0
all_ids = []

while True:
    try:
        with Entrez.esearch(db="protein", term=query, retstart=retstart, retmax=retmax) as handle:
            result = Entrez.read(handle)
        batch = result["IdList"]
        if not batch:
            break
        all_ids.extend(batch)
        print(f"  → Collected {len(all_ids)} IDs", file=sys.stderr)
        if len(batch) < retmax:
            break
        retstart += retmax
        time.sleep(0.34)
    except HTTPError as e:
        print(f"HTTPError: {e}, retrying...", file=sys.stderr)
        time.sleep(10)

if not all_ids:
    sys.exit("No protein IDs found!")

# -------------------------------
# Step 2: Fetch protein sequences
# -------------------------------
records = []
chunks = [all_ids[i:i + 5000] for i in range(0, len(all_ids), 5000)]

for i, chunk in enumerate(chunks):
    print(f"Fetching chunk {i+1}/{len(chunks)}...", file=sys.stderr)
    for attempt in range(3):
        try:
            with Entrez.efetch(db="protein", id=",".join(chunk), rettype="fasta", retmode="text") as handle:
                for rec in SeqIO.parse(handle, "fasta"):
                    records.append((rec.id, str(rec.seq)))
            break
        except HTTPError as e:
            print(f"HTTPError: {e}, retrying...", file=sys.stderr)
            time.sleep(10)
    time.sleep(0.5)

print(f"Fetched {len(records)} protein sequences", file=sys.stderr)

# -------------------------------
# Step 3: Digest to peptides
# -------------------------------
out_records = []
for prot_id, seq in records:
    peptides_by_length = generate_peptides(seq, lengths)
    all_peptides = sum(peptides_by_length.values(), [])
    out_records.append({
        "protein_id": prot_id,
        "num_peptides": len(all_peptides),
        "peptides": ",".join(all_peptides)
    })

df = pd.DataFrame(out_records)
print(f"Generated peptides for {len(df)} proteins", file=sys.stderr)

# -------------------------------
# Step 4: Save result
# -------------------------------
df.to_csv(args.output, sep="\t", index=False)
print(f"Done! Saved to {args.output}")