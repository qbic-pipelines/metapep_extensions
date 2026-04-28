# Metapep – Supplementary Work & Analyses

This folder contains analysis scripts, data processing pipelines, and intermediate results that were developed alongside the main [Metapep pipeline](https://github.com/nf-core/metapep) but are **not yet integrated** into the official repository. It serves as a reference and starting point for future contributors.

The work here covers four main areas: microbiome data preprocessing, pipeline input preparation, self-mimicry analysis, and pipeline validation tooling.


---

## Folder Overview

```
metapep_extensions/
├── DataFactory/          # Preprocessing: raw ASV data → microbiome TSV files
├── PipelineInput/        # Ready-to-use input files for pipeline runs
├── Debugging/            # Core pipeline scripts and intermediate outputs
├── Self_Mimicry/         # Self-mimicry analysis (host proteome cross-reactivity)
└── PipelineEvaluation/   # Tools for comparing pipeline output across versions
```

---

## DataFactory

**Purpose:** Convert raw amplicon sequence variant (ASV) tables into the `taxon_id | abundance` format required by the Metapep pipeline.

**Key file:** [`DataFactory/preprocessDataTables.ipynb`](DataFactory/preprocessDataTables.ipynb)

**What it does:**
1. Loads human gut microbiome ASVs (`humanASV.xlsx`) and soil microbiome ASVs (`soilASV.tsv`)
2. Filters entries with <95% taxonomic confidence at species level
3. Uses `taxonkit` and the NCBI Entrez API to map taxonomy strings to NCBI TaxIDs at strain level
4. Validates that strains have protein/nucleotide sequences available in NCBI
5. Outputs `humanMicrobiome.tsv` and `soilMicrobiome.tsv`

**Output format:**
```tsv
taxon_id    abundance
1226752     0.00437
1210088     0.00098
```

**Dependencies:** `taxonkit` (must be in PATH), BioPython, pandas, NCBI Entrez API access (email + API key)

**Potential for reuse:** Any new microbiome dataset (from 16S/ITS surveys, qPCR, etc.) can be fed through this notebook to produce pipeline-ready input files.

---

## PipelineInput

**Purpose:** Curated, ready-to-use input files for running the Metapep pipeline. This is the "sample library" of test and real conditions.

**Microbiome definition files** (`taxon_id | abundance` TSV):
| File | Description |
|------|-------------|
| `humanMicrobiome.tsv` | 121 strains from human gut microbiome (from DataFactory) |
| `soilMicrobiome.tsv` | 44 strains from soil microbiome (from DataFactory) |
| `commensals.tsv` | Small validation set: known commensal bacteria |
| `pathogens.tsv` | Small validation set: known pathogens |
| `strongBinders.tsv` | Validation: organisms expected to produce strong HLA binders |
| `weakBinders.tsv` | Validation: organisms expected to produce weak HLA binders |
| `debugCondition.tsv` | Minimal test case for debugging |
| `parabacteroidesDistasonis.tsv` | Single-organism test case |
| `testAssemblySelection.tsv` | Tests assembly selection logic |

**Samplesheet CSVs** (pipeline entry point):
| File | Description |
|------|-------------|
| `samplesheet_small_test.csv` | Minimal test run |
| `samplesheet_big_test.csv` | Full-scale test with multiple conditions |
| `samplesheet_proteins.csv` | Direct protein input mode (bypasses NCBI download) |
| `samplesheet_sophie.csv` | Specific configuration used in collaboration |

**Samplesheet format:**
```csv
condition,type,microbiome_path,alleles,weights_path
humanMicrobiome,taxa,PipelineInput/humanMicrobiome.tsv,"A*01:01;A*02:01;B*07:02",
```

**FASTA protein files:** Several bacterial proteomes (Bacteroides sp., Parabacteroides distasonis) for direct protein input mode.

---

## Debugging

**Purpose:** Houses the core pipeline step scripts in standalone form, plus intermediate data files generated during development. Useful for testing individual pipeline steps in isolation without running the full Nextflow pipeline.

**Scripts:**

### `download_proteins_entrez.py`
Downloads bacterial proteomes from NCBI given a list of TaxIDs.
- Selects the largest available assembly per strain
- Traverses NCBI links: Assembly → Nucleotide → Protein
- Output: `proteins.entrez.tsv.gz`, `taxa_assemblies.tsv`, entity/microbiome association tables
- Written by Sabrina Krakau, Leon Kuchenbecker, Till Englert
- CLI: `python download_proteins_entrez.py --microbiomes microbiomes.tsv --email ... --api-key ...`

### `generate_peptides.py`
Digests protein sequences into fixed-length k-mers (peptides).
- Chunks processing by first amino acid for RAM efficiency
- Validates sequences against standard amino acid codes
- Output: `peptides.tsv.gz`, `proteins_peptides.tsv`, `proteins_lengths.tsv`

### `generate_protein_and_entity_ids.py`
Merges predicted proteins (from genome assemblies) and downloaded NCBI proteins into a unified ID scheme.
- Assigns global protein IDs and entity IDs
- Output: `proteins.tsv.gz`, `entities_proteins.tsv`, `entities.tsv`, `microbiomes_entities.tsv`

**Intermediate data files** (can be used to re-run downstream steps without redownloading):
- `proteins.tsv.gz` / `proteins.entrez.tsv.gz` – protein sequences
- `peptides.tsv.gz` – generated peptides
- `proteins_peptides.tsv` – protein-to-peptide mapping
- `entities.tsv` / `entities_proteins.tsv` – entity/protein ID mappings

---

## Self_Mimicry

**Purpose:** Analyze cross-reactivity between microbiome-derived peptides and the human proteome (self-mimicry). This is a complete standalone analysis that is **not part of the main Metapep pipeline**.

**Core idea:** For each microbial entity and HLA allele, compute the fraction of HLA-binding peptides that also appear in the human proteome. A high self-mimicry ratio may indicate potential for autoimmune cross-reactivity or tolerance induction.

### Scripts

#### `download_and_digest_host_proteome.py`
Downloads the complete human proteome (TaxID 9606) from NCBI and digests it into peptides of specified length.
- Output: `proteome_peptides.tsv`
- CLI: `python download_and_digest_host_proteome.py -t 9606 -e your@email.com -k API_KEY -l 9 -o proteome_peptides.tsv`

#### `prepare_self_mimicry_ratios.py`
Prepares entity-peptide tables for self-mimicry calculation.
- Filters peptides by HLA binding prediction score (default threshold: 0.426)
- Maps filtered peptide IDs → entity IDs → actual sequences
- Chunked processing (2–5M rows/chunk) for large datasets
- Output per allele: `entities_peptides.allele_N.tsv`

#### `compute_self_mimicry_ratio.py`
Calculates the self-mimicry ratio for each entity.
- For each entity: `self_mimicry_ratio = shared_peptides_with_host / total_predicted_binding_peptides`
- Output per allele: `entities_self_mimicry_ratios.allele_N.tsv`

#### `plot_self_mimicry.R`
Produces scatter plots of binding ratio vs. self-mimicry ratio per HLA allele with Pearson correlation.
- Input: `entity_binding_ratios.allele_N.tsv` + `entities_self_mimicry_ratios.allele_N.tsv`

#### `plot_entity_binding_ratios.R`
Box plots of entity-wise HLA binding rates grouped by condition (e.g., soil vs. human microbiome).
- Optional Mann-Whitney U test p-values between conditions
- Output: PDFs per allele

#### `plot_entity_binding_ratios.nf`
Nextflow process wrapper for `plot_entity_binding_ratios.R` – ready to be plugged into the main pipeline.

**Pre-computed results** (HLA alleles A\*01:01, A\*02:01, B\*07:02, B\*08:01; conditions: humanMicrobiome, soilMicrobiome):
- `entity_binding_ratios.allele_[0-3].tsv`
- `entities_self_mimicry_ratios.allele_[0-3].tsv`
- PNG plots for each allele

**Potential for integration into main pipeline:** The Nextflow process in `plot_entity_binding_ratios.nf` and the self-mimicry scripts are candidates for inclusion in the official Metapep pipeline as optional output modules.

---

## PipelineEvaluation

**Purpose:** Validate that two versions of the Metapep pipeline produce identical (or near-identical) outputs. Used to verify refactoring, fork merges, or parameter changes don't introduce regressions.

### `compare_test_results.py`
Compares two pipeline output directories across multiple test conditions.
- File size and MD5 checksum comparison
- Order-independent line-by-line content hashing (works for TSV/gzip files where row order may differ)
- Streaming/chunked comparison (RAM-safe for large compressed files)
- Produces a detailed report of unique vs. shared entries

```bash
python3 compare_test_results.py \
    --dev output_fork2 \
    --original output_fork \
    --lightweight
```

**Existing comparison logs** in the folder show validation results for several test scenarios (default, bins-only, co-assembly, mouse), with ~99.99–100% identity between pipeline versions.

---

## Recommended Next Steps / Integration Opportunities

| Component | Status | Notes |
|-----------|--------|-------|
| DataFactory notebook | Standalone | Could be wrapped as a preprocessing script for the pipeline's documentation |
| Self-mimicry scripts | Standalone, complete | Strong candidate for an optional pipeline module |
| `plot_entity_binding_ratios.nf` | Nextflow process, ready | Can be directly added to `modules/local/` in the main pipeline |
| `compare_test_results.py` | Standalone test utility | Could be added to the pipeline's CI test suite |
| Validation datasets in PipelineInput | Ready | Already usable as CI test conditions |

---

## Dependencies

| Tool | Used in | Install |
|------|---------|---------|
| Python ≥3.8 | All scripts | — |
| pandas | All scripts | `pip install pandas` |
| BioPython | Protein download/digestion | `pip install biopython` |
| taxonkit | DataFactory | [taxonkit.bioinf.shenwei.me](https://bioinf.shenwei.me/taxonkit/) |
| R ≥4.0 | Self_Mimicry plots | — |
| ggplot2, dplyr, ggpubr | R scripts | `install.packages(c("ggplot2","dplyr","ggpubr"))` |
| Nextflow | `.nf` process | [nextflow.io](https://www.nextflow.io/) |
| NCBI Entrez API key | Download scripts | Register at [ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov/account/) |

---

## Author

Florian Kramer  
Contact: flo_kramer@t-online.de
