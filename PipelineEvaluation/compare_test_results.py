#!/usr/bin/env python3
import os
import gzip
import hashlib
import csv
import argparse
import sys
from typing import Optional, List, Tuple

# Run: python3 compare_test_results.py --dev ../output_fork2_small_test --original ../output_fork_small_test --lightweight


# =========================
# Tee logger
# =========================
class Tee:
    """Write stdout to both console and file."""
    def __init__(self, filename):
        self.file = open(filename, "w")
        self.stdout = sys.stdout

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        self.file.flush()
        self.stdout.flush()


# =========================
# File helpers
# =========================
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


# =========================
# Lightweight streaming hash
# =========================
def stream_hash(path: str, columns: Optional[List[int]]):
    h = hashlib.md5()
    line_count = 0

    with open_maybe_gzip(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        for row in reader:
            if columns is None:
                data = "\t".join(row)
            else:
                data = "\t".join(row[i] for i in columns if i < len(row))

            h.update(data.encode("utf-8"))
            h.update(b"\n")
            line_count += 1

    return header, line_count, h.hexdigest()


# =========================
# Heavy comparison helpers (RAM-safe)
# =========================
def collect_keys_and_examples(
    path: str,
    columns: Optional[List[int]],
    max_examples: int = 10
) -> Tuple[list, set, list, int]:
    """
    Returns:
        header
        unique_keys
        example_rows (full rows, limited)
        line_count
    """
    keys = set()
    examples = []
    line_count = 0

    with open_maybe_gzip(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        for row in reader:
            line_count += 1

            if columns is None:
                key = tuple(cell.strip() for cell in row)
            else:
                key = tuple(row[i].strip() for i in columns if i < len(row))

            if key not in keys:
                keys.add(key)
                if len(examples) < max_examples:
                    examples.append(tuple(row))

    return header, keys, examples, line_count


def collect_diff_examples(
    path: str,
    columns: Optional[List[int]],
    diff_keys: set,
    max_examples: int = 10
) -> list:
    """
    Second streaming pass to collect example rows for specific keys.
    """
    examples = []

    with open_maybe_gzip(path) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)  # skip header

        for row in reader:
            if columns is None:
                key = tuple(cell.strip() for cell in row)
            else:
                key = tuple(row[i].strip() for i in columns if i < len(row))

            if key in diff_keys:
                examples.append(tuple(row))
                if len(examples) >= max_examples:
                    break

    return examples


def print_example_lines(title: str, rows: list, header: list, max_examples: int = 10):
    if not rows:
        return

    print(f"[INFO] {title} (showing up to {max_examples} examples):")
    print("       " + "\t".join(header))
    for row in rows[:max_examples]:
        print("       " + "\t".join(row))


# =========================
# Compare functions
# =========================
def compare_pair_lightweight(name, dev_path, original_path, columns: Optional[List[int]]):
    print(f"\n=== Lightweight comparing {name} ===")

    print(f"[INFO] Dev file size: {file_size(dev_path)} bytes")
    print(f"[INFO] Original file size: {file_size(original_path)} bytes")

    print("[INFO] Streaming comparison ...")
    hdr_d, cnt_d, hash_d = stream_hash(dev_path, columns)
    hdr_o, cnt_o, hash_o = stream_hash(original_path, columns)

    print(f"[INFO] Dev line count: {cnt_d}")
    print(f"[INFO] Original line count: {cnt_o}")

    if hdr_d != hdr_o:
        print("[WARNING] Headers differ!")
        print(f"Dev:      {hdr_d}")
        print(f"Original: {hdr_o}")

    print("[INFO] Content hashes:")
    print(f"       Dev:      {hash_d}")
    print(f"       Original: {hash_o}")

    if cnt_d == cnt_o and hash_d == hash_o:
        print("[RESULT] Files are IDENTICAL (lightweight)")
    else:
        print("[RESULT] Files DIFFER (lightweight)")


def compare_pair_heavy(name, dev_path, original_path, columns: Optional[List[int]], numeric):
    print(f"\n=== Comparing {name} ===")

    print(f"[INFO] Dev file size: {file_size(dev_path)} bytes")
    print(f"[INFO] Original file size: {file_size(original_path)} bytes")

    print("[INFO] Computing MD5 checksums ...")
    print(f"[INFO] Dev MD5: {compute_md5(dev_path)}")
    print(f"[INFO] Original MD5: {compute_md5(original_path)}")

    if numeric:
        dev_vals, orig_vals = {}, {}

        with open_maybe_gzip(dev_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            key_col = reader.fieldnames[0]
            for row in reader:
                dev_vals[row[key_col]] = {
                    k: int(v) for k, v in row.items() if k != key_col
                }

        with open_maybe_gzip(original_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            key_col = reader.fieldnames[0]
            for row in reader:
                orig_vals[row[key_col]] = {
                    k: int(v) for k, v in row.items() if k != key_col
                }

        for cond in sorted(set(dev_vals) & set(orig_vals)):
            print(f"[INFO] Condition: {cond}")
            for col in dev_vals[cond]:
                d, o = dev_vals[cond][col], orig_vals[cond][col]
                pct = (min(d, o) / max(d, o) * 100) if max(d, o) else 100
                print(f"       {col}: dev={d}, orig={o}, overlap={pct:.2f}%")
        return

    print("[INFO] Collecting unique keys (streaming) ...")
    header_dev, dev_keys, _, dev_lines = collect_keys_and_examples(dev_path, columns)
    header_orig, orig_keys, _, orig_lines = collect_keys_and_examples(original_path, columns)

    shared = dev_keys & orig_keys
    only_dev = dev_keys - orig_keys
    only_orig = orig_keys - dev_keys

    print(f"[INFO] Dev line count: {dev_lines}")
    print(f"[INFO] Original line count: {orig_lines}")

    print("[INFO] Line comparison:")
    print(f"       Unique dev: {len(dev_keys)}")
    print(f"       Unique original: {len(orig_keys)}")
    print(f"       Shared: {len(shared)}")
    print(f"       Only dev: {len(only_dev)}")
    print(f"       Only original: {len(only_orig)}")

    if only_dev:
        dev_examples = collect_diff_examples(dev_path, columns, only_dev)
        print_example_lines(
            "Entries present in dev but NOT in original",
            dev_examples,
            header_dev,
        )

    if only_orig:
        orig_examples = collect_diff_examples(original_path, columns, only_orig)
        print_example_lines(
            "Entries present in original but NOT in dev",
            orig_examples,
            header_orig,
        )


# =========================
# Main
# =========================
def main():
    parser = argparse.ArgumentParser(description="Compare output files from two pipelines")
    parser.add_argument("--dev", required=True)
    parser.add_argument("--original", required=True)
    parser.add_argument("--lightweight", action="store_true")
    args = parser.parse_args()

    dev_name = os.path.basename(os.path.normpath(args.dev))
    orig_name = os.path.basename(os.path.normpath(args.original))
    log_file = f"comparison_log_{dev_name}_vs_{orig_name}.txt"

    sys.stdout = Tee(log_file)

    pairs = [
        ("proteins.tsv.gz", [0],
         f"{args.dev}/db_tables/proteins.tsv.gz",
         f"{args.original}/db_tables/proteins.tsv.gz", False),

        ("peptides.tsv.gz", [0],
         f"{args.dev}/db_tables/peptides.tsv.gz",
         f"{args.original}/db_tables/peptides.tsv.gz", False),

        ("predictions.tsv.gz", [0],
         f"{args.dev}/db_tables/predictions.tsv.gz",
         f"{args.original}/db_tables/predictions.tsv.gz", False),

        ("proteins_peptides.tsv", None,
         f"{args.dev}/db_tables/proteins_peptides.tsv",
         f"{args.original}/db_tables/proteins_peptides.tsv", False),

        ("stats.tsv", None,
         f"{args.dev}/db_tables/stats.tsv",
         f"{args.original}/db_tables/stats.tsv", True),
    ]

    print(f"[INFO] Starting comparison (log saved to {log_file})")

    for name, columns, dev_file, orig_file, numeric in pairs:
        if not os.path.exists(dev_file):
            print(f"[ERROR] Missing dev file: {dev_file}")
            continue
        if not os.path.exists(orig_file):
            print(f"[ERROR] Missing original file: {orig_file}")
            continue

        if args.lightweight:
            compare_pair_lightweight(name, dev_file, orig_file, columns)
        else:
            compare_pair_heavy(name, dev_file, orig_file, columns, numeric)

    print("\n[INFO] Finished all comparisons.")


if __name__ == "__main__":
    main()
