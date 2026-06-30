<p align="center">
  <h1 align="center">AsmSV</h1>
  <p align="center">
    <strong>Assembly-based Structural Variant Caller for Haplotype-resolved Genomes</strong>
  </p>
  <p align="center">
    <a href="#installation">Installation</a> &nbsp;·&nbsp;
    <a href="#quick-start">Quick Start</a> &nbsp;·&nbsp;
    <a href="#method-overview">Method</a> &nbsp;·&nbsp;
    <a href="#output-format">Output</a> &nbsp;·&nbsp;
    <a href="#parameters">Parameters</a> &nbsp;·&nbsp;
    <a href="#examples">Examples</a> &nbsp;·&nbsp;
    <a href="#citation">Citation</a>
  </p>
</p>

---

## Introduction

**AsmSV** detects structural variants (SVs) from haplotype-resolved *de novo* assemblies aligned to a reference genome. It takes two BAM files — one per haplotype — and produces a standard VCFv4.2 with phased genotypes (`1|0`, `0|1`, `1|1`).

AsmSV is built on a **Backbone–Branch strategy**: a dynamic-programming chain finds the optimal collinear alignment backbone, and off-backbone ("branch" / "orphan") segments are analyzed to resolve complex rearrangements — INS+DEL, composite inversions, and multi-breakpoint compound events — providing a comprehensive, accurate, genome-wide view of structural variation.

### Key Features

- **Comprehensive SV spectrum** — DEL, INS, DUP (tandem / proximal / dispersed), INV, BND/TRA, and compound CSV (complex sequence variant) types including INV-INS, INV-DEL, InvDUP, and INS_DEL replacements.
- **Phased genotyping && Fine‑grained traceable annotation fields** — Cross‑haplotype clustering generates phased genotypes directly from diploid assemblies, alongside standard VCF‑v4.2 output enriched with traceable INFO‑field annotations to facilitate in‑depth interrogation of individual complex structural variants (CSVs).
- **Complex event decomposition** — multi-branch events are decomposed into individual SV components linked by shared `EVENT` IDs, plus a macro `*_SUMMARY` record per complex locus. Use the bundled `example/plot.py` to visualize CRR (complex rearrangement region) topology.
- **Fast & lightweight** — multi-threaded with safe IPC serialization of alignment records. On a 3950X (16-core / 128 GB) at default 8 threads, one human genome takes ~74 s (index reused, 6.2 GB) or ~207 s (index rebuilt, 8.0 GB).
- **Flexible filtering** — BED-based region masking, sliding-window density filtering, and configurable size/quality thresholds.

---

## Installation

### Requirements

- Python ≥ 3.8
- [pysam](https://github.com/pysam-developers/pysam) ≥ 0.19 (tested on 0.23.3)
- [edlib](https://github.com/Martus901/edlib) (tested on 1.3.9)
- [samtools](https://github.com/samtools/samtools) (for BAM sorting, tested on 1.17+)

### Install from source

```bash
git clone https://github.com/Milia1368/AsmSV.git
cd AsmSV
pip install -e .
```

This exposes the `asmSV` command-line entry point.

---

## Quick Start

### Input preparation

Align each haplotype assembly to the reference with a long-read aligner such as [minimap2](https://github.com/lh3/minimap2):

```bash
minimap2 -a --eqx -x asm5 reference.fa hap1.fa | samtools sort -o hap1_to_ref.bam
minimap2 -a --eqx -x asm5 reference.fa hap2.fa | samtools sort -o hap2_to_ref.bam
samtools index hap1_to_ref.bam
samtools index hap2_to_ref.bam
```

### Basic usage

```bash
asmSV \
  --bam1 hap1_to_ref.bam \
  --bam2 hap2_to_ref.bam \
  --ref  reference.fa \
  -o     results.vcf \
  -t     8
```

### Reuse pre-sorted BAMs (skip sorting)

```bash
asmSV \
  --bam1 hap1.bam --sorted-bam1 hap1.qsort.bam \
  --bam2 hap2.bam --sorted-bam2 hap2.qsort.bam \
  --ref  reference.fa \
  -o     results.vcf
```
### Resort BAMs (force sorting)

```bash
asmSV \
  --bam1 hap1.bam \
  --bam2 hap2.bam  \
  --ref  reference.fa \
  -o     results.vcf --force-resort
```

### With region masking

```bash
asmSV \
  --bam1 hap1.bam --bam2 hap2.bam \
  --ref  reference.fa -o results.vcf \
  --mask repeats.bed
```

### Tune size thresholds

```bash
asmSV \
  --bam1 hap1.bam --bam2 hap2.bam \
  --ref  reference.fa -o results.vcf \
  --min-sv 50 --max-sv 2000000 --min-csv 500
```

---

## Method Overview

AsmSV processes each haplotype independently through the pipeline below, then merges results across haplotypes for genotyping.

```
BAM (Hap1 / Hap2)
     │
     ▼
Sort by query name (or use --sorted-bam)
     │
     ▼
Group alignments by read name
  · MAPQ ≥ 10 filtering
  · Extract forward-strand query sequence
     │
     ▼
Per-read SV extraction (parallel across threads)
  ├─ split_segments           # resolve query overlap (optional)
  ├─ solve_chain_dp           # DP optimal alignment chain
  ├─ find_collinear_backbone  # identify uni-directional backbone
  ├─ resolve_chain_svs        # backbone–branch SV extraction
  │   ├─ Scenario B: no branches → linear SV (DEL/INS/DUP/INS_DEL)
  │   └─ Scenario A: branches present → complex SV resolution
  │       ├─ [A1] single inv. branch → INV / CSV(INV-INS/DEL) / InvDUP
  │       ├─ [A2] nested linear SVs across branch network
  │       └─ [A3] macro topology → BND pairs + span records
  └─ parse_intra_alignment_sv # CIGAR-internal INDELs
     │
     ▼
BED mask filtering (optional)
Density filtering (optional)
     │
     ▼
Cross-haplotype clustering & genotyping
     │
     ▼
VCFv4.2 output
```

### Backbone–Branch Strategy

1. **DP Chaining** (`solve_chain_dp`) — Segments are sorted by query coordinate and chained via a penalized longest-path DP (similar to minimap2 chaining). Penalties account for query gaps/overlaps, reference distance, strand flips, and cross-chromosome jumps.

2. **Backbone Identification** (`find_collinear_backbone`) — A second DP pass selects the dominant uni-directional chain. Segments whose strand disagrees with the global majority are forced into branch status. An adaptive overlap filter (dual thresholds: 5 kb ceiling, 50% ratio) prevents alignment artifacts from corrupting the backbone.

3. **SV Extraction** (`resolve_chain_svs`) — Adjacent backbone pairs define gaps. Branch segments within these gaps reveal complex rearrangements:
   - **Scenario B** (no branches): linear SVs — DEL, INS, DUP:TANDEM, or CSV (INS_DEL replacement).
   - **Scenario A** (branches present): single inverted branches trigger local inversion analysis with flanking INDEL decomposition (**[A1]**); all branch–backbone transitions are scanned for nested linear SVs (**[A2]**); branches are condensed into typed clusters and emitted as BND pairs plus span records (**[A3]**).

### Edge Refinement

Each alignment undergoes automatic edge refinement at construction time: CIGAR ends are scanned for stable anchors (≥10 bp match with positive neighborhood score), and unstable flanking regions are trimmed to prevent false-positive INDELs.

---

## Output Format

AsmSV writes standard VCFv4.2. All records include `SVTYPE`, `SVLEN`, `END`, and `GT`.

### Structural Variant Types

| SVTYPE | Description |
|--------|-------------|
| DEL    | Deletion |
| INS    | Insertion (includes short tandem duplications < 1 kb) |
| DUP    | Tandem duplication ≥ 1 kb |
| INV    | Inversion |
| BND    | Breakend / novel adjacency |
| CSV    | Complex Sequence Variant (see below) |

### CSV — Complex Sequence Variant

CSV is a custom SVTYPE for compound events. The `SUBTYPE` INFO field specifies the mechanism:

| SUBTYPE | Description | Source |
|---------|-------------|--------|
| INS_DEL | Simultaneous insertion + deletion (sequence replacement) | Backbone linear |
| INV-INS | Inversion with flanking insertion | Branch [A1] |
| INV-DEL | Inversion with flanking deletion | Branch [A1] |
| INV-INS-DEL | Inversion with flanking insertion and deletion | Branch [A1] |
| InvDUP:TANDEM | Local inverted tandem duplication | Branch [A1] |
| InvDUP:PROXIMAL | Proximal inverted dispersed duplication | Branch [A3] |
| InvDUP:DISP:INTRA | Distal intra-chromosomal inverted duplication | Branch [A1]+[A3] |
| DUP:PROXIMAL | Proximal same-strand dispersed duplication | Branch [A3] |
| DUP:DISP:INTRA | Distal intra-chromosomal same-strand duplication | Branch [A3] |
| DUP:DISP:INTER / InvDUP:DISP:INTER | Inter-chromosomal dispersed duplication breakend pair | Branch [A3] |
### Branch Classification Logic

Branch blocks are classified by distance to the backbone gap and strand orientation:

| Priority | Condition | Same strand | Reverse strand | Cross-chr |
|----------|-----------|-------------|----------------|-----------|
| 0 (highest) | Fully inside gap | LOCAL_RETAIN | InvDUP | — |
| 1 | ≤ 500 bp to gap | DUP:TANDEM | InvDUP:TANDEM | — |
| 2 | ≤ 10 kb to gap | DUP:PROXIMAL | InvDUP:PROXIMAL | — |
| 3 (default) | > 10 kb to gap | DUP:DISP:INTRA | InvDUP:DISP:INTRA | [Inv]DUP:DISP:INTER |

> **Note:** `LOCAL_RETAIN` (priority 0, same strand) produces no span record — the cluster is tracked only via the flanking BND pair. `InvDUP` (priority 0, reverse strand) produces a `CSV` record (`SUBTYPE=InvDUP`) anchored at the orphan's actual reference coordinates.

### Core INFO Fields

| Field | Description |
|-------|-------------|
| SVTYPE | Structural variant type |
| SVLEN | SV length (negative for DEL) |
| END | End position on reference |
| SUBTYPE | Variant subclass (see tables above) |
| MATEID | Mate BND record ID for breakend pairing |
| SUPPORT_H1 / SUPPORT_H2 | Per-haplotype support counts |
| CHR2 / POS2 / STRANDS | BND partner coordinates and strand orientation |
| DEL_LEN / INS_LEN | Deleted and inserted lengths for CSV events |
| SEQ | Inserted sequence |
| EVENT | Shared event ID linking related sub-records |
| SRCREG | Source genomic region of variant sequence |
| INV_SPAN | Span of inner inverted region |
| SOURCE | Detection layer (backbone / CIGAR / COMPLEX / COMPLEX_EDGE) |
| RNAMES | Supporting read names |
| SUB_SVS | Pipe-separated topology trace on COMPLEX `*_SUMMARY` records |

### FILTER Tags

| Tag | Description |
|-----|-------------|
| PASS | High-confidence variant |
| COMPLEX | Multi-cluster (≥2 branch clusters) complex SV summary record |
| LOCAL_BND | Short-range intra-chromosomal BND (not a true translocation; may be upgraded to PASS when supported by both haplotypes) |
| OVERSIZE | SV exceeds `--max-sv` threshold |
| EDGE_COMPLEX | Complex rearrangement at read edge (orphans before/after the backbone) |

### Genotypes

| GT | Meaning |
|----|---------|
| `1\|0` | Variant on Haplotype 1 only |
| `0\|1` | Variant on Haplotype 2 only |
| `1\|1` | Variant on both haplotypes (homozygous) |
| `./.` | Unknown |

---

## Parameters

### Input / Output

| Option | Default | Description |
|--------|---------|-------------|
| `--bam1` *(required)* | — | Input BAM for Hap1 (sorted & indexed) |
| `--bam2` *(required)* | — | Input BAM for Hap2 (sorted & indexed) |
| `--ref` *(required)* | — | Reference genome FASTA (indexed) |
| `-o, --output` | `./results.vcf` | Output VCF path (directory → `results.vcf` inside) |
| `--mask` | none | BED file of regions to exclude from output |
| `--sorted-bam1` | none | Pre-sorted (queryname) BAM for Hap1; skips sorting |
| `--sorted-bam2` | none | Pre-sorted (queryname) BAM for Hap2; skips sorting |

### VCF Format

| Option | Default | Description |
|--------|---------|-------------|
| `--enable-large-sequence` | off | Always write full ALT/REF sequences (otherwise symbolic alleles above 10 kb) |

### Filtering

| Option | Default | Description |
|--------|---------|-------------|
| `--min-sv` | 30 | Minimum SV span (bp) to report |
| `--max-sv` | 1500000 | Maximum SV span; larger events tagged `OVERSIZE` (`-1` to disable) |
| `--min-mapq` | 10 | Minimum mapping quality |
| `--insdel-tolerance` | 10 | INS vs DEL length-difference tolerance (bp) |
| `--min-tra-size` | 150000 | Min ref gap to classify as TRA/BND |
| `--min-csv` | 400 | Min size (bp) for INS_DEL replacement events |
| `--max-csv` | 40000 | Max size (bp) for INS_DEL replacement events |
| `--split-segments` | off | Split overlapping alignment segments before DP chaining |

### Clustering

| Option | Default | Description |
|--------|---------|-------------|
| `--merge-indels-tol` | 3 | Distance tolerance to merge adjacent identical INDELs |
| `--cluster-pos-tol` | 3 | Max breakpoint distance for merging |
| `--cluster-len-tol` | 0 | Length tolerance adjustment |
| `--min-seq-sim` | 1.0 | Min sequence similarity for INS merging |

### Advanced

| Option | Default | Description |
|--------|---------|-------------|
| `-t, --threads` | 8 | Parallel processing threads |
| `--enable-density-filter` | off | Suppress dense clusters of candidates |
| `--density-window` | 1000 | Window size (bp) for density filtering |
| `--max-sv-in-window` | 3 | Max SVs allowed per density window |
| `--no-retain-sorted` | off | Delete intermediate `.qsort.bam` after processing |
| `--force-resort` | off | Force re-sort even if a valid `.qsort.bam` exists |
| `--debug` | off | Verbose DEBUG logging |

### Algorithm Details

**DP Chain Scoring** — For each segment *i*, the algorithm considers the *K* = 32 nearest predecessors *j*:

```
score[i] = max { score[j] + AS[i] - penalty(i, j) }
```

| Condition | Penalty |
|-----------|---------|
| Query gap > 0 | gap × 0.25 |
| Query overlap | overlap × density × 1.0 |
| Same-strand ref gap < 10 kb | distance × 0.1 |
| Same-strand ref gap ≥ 10 kb | +500 |
| Strand flip (INV) | +200 |
| Cross-chromosome | +500–600 |

The global score additionally penalizes unmapped query regions at both ends.

**Clustering Rules**

| SV Type | Matching Criteria |
|---------|-------------------|
| DEL / DUP / INV | Δpos ≤ 3 bp, Δlen < 2 bp |
| INS | Δpos ≤ 3 bp, Δlen ≤ 0 bp, seq similarity ≥ 1.0 |
| BND / TRA | Both Δpos ≤ 3 bp, same CHR2 |
| CSV (INS_DEL / INV-*) | Δpos ≤ 3 bp, ΔEND ≤ 3 bp, INS_LEN + DEL_LEN within tolerance or 10% |
| CSV (DUP:* / InvDUP:*) | Δpos ≤ 3 bp, Δlen < 2 bp |

> CSV detection involves several empirical parameters — adjust them to suit your data. Ideas, requirements, and algorithm improvements are very welcome!

---

## Examples

The `example/` directory provides helper scripts for filtering VCF output and benchmarking. All accept `-h` for full help; outputs default to a derived name next to the input when `-o` is omitted.

### Filter VCF for benchmarking

Keep PASS single-type SVs; drop CSV by default.

```bash
python example/filter_sv.py -i results.vcf
python example/filter_sv.py -i results.vcf -o clean.vcf --keep-csv
```

### Decompose CSV into simple SVs

```bash
# INV-related CSV → INV / DEL / INS
python example/trans_inv.py -i results.vcf --min-sv 50

# INS_DEL CSV → DEL + INS
python example/trans_CINS_DEL.py -i results.vcf
```


### Visualize a COMPLEX region

Plot the 2D topology (Ref vs Assembly) of a `*_SUMMARY` record from its `SUB_SVS` field. Requires `matplotlib`.

```bash
# List plottable COMPLEX events
python example/plot.py -i results.vcf --list

# Plot one event
python example/plot.py -i results.vcf -e haplotype1-0000097_CPLX_2 --hap Hap1 -o out.png
```

---

## Project Structure

```
AsmSV/
├── src/AsmSV/
│   ├── asmSV.py        # Entry point & per-haplotype pipeline
│   ├── bam_info.py     # Backbone–branch SV extraction & complex event resolution
│   ├── structures.py   # Core data structures (alignment nodes, clusters, candidates)
│   ├── genotype.py     # Cross-haplotype clustering & phased genotyping
│   ├── utils.py        # CLI parsing, VCF I/O, helpers
│   └── config.py       # Global parameters & defaults
├── example/            # Plotting, filtering & benchmarking helpers
├── setup.py
├── LICENSE
└── README.md
```

---

## Citation

If you use AsmSV in your research, please cite:

> *Citation information will be added upon publication.*

---

## License

This project is licensed under the [MIT License](LICENSE).
