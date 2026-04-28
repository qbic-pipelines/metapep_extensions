#!/usr/bin/env python3
import os
import gzip
import hashlib
import csv
import argparse
import sys
from collections import defaultdict
from typing import List, Optional

# ============================================================
# Tee logger
# ============================================================
class Tee:
    def __init__(self, filename):
        self.file = open(filename, "w")
        self.stdout = sys.stdout

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        self.file.flush()
        self.stdout.flush()


# ============================================================
# File helpers
# ============================================================
def is_gzip(path: str) -> bool:
    return path.endswith(".gz")


def open_maybe_gzip(path: str, mode="rt"):
    return gzip.open(path, mode) if is_gzip(path) else open(path, mode)


def file_size(path: str) -> int:
    return os.path.getsize(path)


def compute_md5(path: str) -> str:
    h = hashlib.md5()
    with open_maybe_gzip(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


# ============================================================
# Streaming diff core
# ============================================================
def stream_keys(
    path: str,
    key_columns: Optional[List[str]],
    max_examples: int = 10,
):
    """
    Returns:
        header
        key_counter
        example_rows (key -> original row)
        line_count
    """
    key_counter = defaultdict(int)
    example_rows = {}
    line_count = 0

    with open_maybe_gzip(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames

        if key_columns:
            for col in key_columns:
                if col not in header:
                    raise ValueError(f"Key column '{col}' not found in {path}. Available: {header}")

        for row in reader:
            line_count += 1

            if key_columns:
                key = tuple(row[col].strip() for col in key_columns)
            else:
                key = tuple(row[h].strip() for h in header)

            key_counter[key] += 1

            if key not in example_rows and len(example_rows) < max_examples:
                example_rows[key] = row

    return header, key_counter, example_rows, line_count


def print_examples(title, rows, header):
    if not rows:
        return
    print(f"[INFO] {title}:")
    print("       " + "\t".join(header))
    for row in rows:
        print("       " + "\t".join(row[h] for h in header))


# ============================================================
# Comparison
# ============================================================
def compare_pair(name, dev_path, orig_path, key_columns=None, numeric=False, max_examples=10):
    print(f"\n=== Comparing {name} ===")

    print(f"[INFO] Dev file size: {file_size(dev_path)} bytes")
    print(f"[INFO] Original file size: {file_size(orig_path)} bytes")

    print("[INFO] Computing MD5 checksums ...")
    print(f"[INFO] Dev MD5: {compute_md5(dev_path)}")
    print(f"[INFO] Original MD5: {compute_md5(orig_path)}")

    if numeric:
        print("[INFO] Numeric comparison ...")
        dev_vals, orig_vals = {}, {}

        with open_maybe_gzip(dev_path) as f:
            r = csv.DictReader(f, delimiter="\t")
            key = r.fieldnames[0]
            for row in r:
                dev_vals[row[key]] = { k: round(float(v), 5) for k, v in row.items() if k != key }

        with open_maybe_gzip(orig_path) as f:
            r = csv.DictReader(f, delimiter="\t")
            key = r.fieldnames[0]
            for row in r:
                orig_vals[row[key]] = { k: round(float(v), 5) for k, v in row.items() if k != key }

        for cond in sorted(set(dev_vals) & set(orig_vals)):
            print(f"[INFO] Condition: {cond}")
            for col in dev_vals[cond]:
                d, o = dev_vals[cond][col], orig_vals[cond][col]
                pct = (min(d, o) / max(d, o) * 100) if max(d, o) else 100
                print(f"       {col}: dev={d}, orig={o}, overlap={pct:.2f}%")
        return

    print("[INFO] Streaming keys ...")
    hdr_d, dev_keys, dev_examples, dev_lines = stream_keys(dev_path, key_columns, max_examples)
    hdr_o, orig_keys, orig_examples, orig_lines = stream_keys(orig_path, key_columns, max_examples)

    dev_set = set(dev_keys)
    orig_set = set(orig_keys)

    shared = dev_set & orig_set
    only_dev = dev_set - orig_set
    only_orig = orig_set - dev_set

    print(f"[INFO] Dev line count: {dev_lines}")
    print(f"[INFO] Original line count: {orig_lines}")

    print("[INFO] Line comparison:")
    print(f"       Unique dev keys: {len(dev_set)}")
    print(f"       Unique original keys: {len(orig_set)}")
    print(f"       Shared keys: {len(shared)}")
    print(f"       Only dev: {len(only_dev)}")
    print(f"       Only original: {len(only_orig)}")

    print_examples("Entries present in dev but NOT in original", [dev_examples[k] for k in only_dev if k in dev_examples], hdr_d)
    print_examples("Entries present in original but NOT in dev", [orig_examples[k] for k in only_orig if k in orig_examples], hdr_o)


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dev", required=True)
    parser.add_argument("--original", required=True)
    args = parser.parse_args()

    dev_name = os.path.basename(os.path.normpath(args.dev))
    orig_name = os.path.basename(os.path.normpath(args.original))
    log_file = f"comparison_log_{dev_name}_vs_{orig_name}.txt"

    sys.stdout = Tee(log_file)

    pairs = [
        ("proteins.tsv.gz", ["protein_id"], f"{args.dev}/db_tables/proteins.tsv.gz", f"{args.original}/db_tables/proteins.tsv.gz", False),
        ("peptides.tsv.gz", ["peptide_sequence"], f"{args.dev}/db_tables/peptides.tsv.gz", f"{args.original}/db_tables/peptides.tsv.gz", False),
        ("predictions.tsv.gz", None, f"{args.dev}/db_tables/predictions.tsv.gz", f"{args.original}/db_tables/predictions.tsv.gz", False),
        ("proteins_peptides.tsv", ["peptide_id"], f"{args.dev}/db_tables/proteins_peptides.tsv", f"{args.original}/db_tables/proteins_peptides.tsv", False),
        ("stats.tsv", None, f"{args.dev}/db_tables/stats.tsv", f"{args.original}/db_tables/stats.tsv", True),
    ]

    print(f"[INFO] Starting comparison (log saved to {log_file})")

    for name, key_cols, dev_file, orig_file, numeric in pairs:
        if not os.path.exists(dev_file) or not os.path.exists(orig_file):
            print(f"[ERROR] Missing file for {name}")
            continue

        compare_pair(name, dev_file, orig_file, key_columns=key_cols, numeric=numeric)

    print("\n[INFO] Finished all comparisons.")


if __name__ == "__main__":
    main()
