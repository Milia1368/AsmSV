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
    <a href="#citation">Citation</a>
  </p>
</p>

---

## Introduction

**AsmSV** detects structural variants (SVs) from haplotype-resolved *de novo* assemblies aligned to a reference genome. It takes two BAM files — one per haplotype — and produces a standard VCFv4.2 with phased genotypes (`1|0`, `0|1`, `1|1`).

AsmSV is designed around a **Backbone–Branch strategy**: a dynamic-programming chain finds the optimal collinear alignment backbone, and off-backbone ("branch" or "orphan") segments are analyzed to resolve complex rearrangements including inversions, inverted/dispersed duplications, and multi-breakpoint compound events.

### Key Features

- **Comprehensive SV spectrum** — DEL, INS, DUP (tandem / proximal / dispersed), INV, BND/TRA, and compound CSV (complex sequence variant) types including INV-INS, INV-DEL, InvDUP, and INS_DEL replacements.
- **Phased genotyping** — cross-haplotype clustering produces phased genotypes directly from diploid assemblies.
- **Complex event decomposition** — multi-branch-cluster events are decomposed into individual SV components linked by shared `EVENT` IDs, with a macro summary record for each complex locus.
- **Parallel execution** — multi-threaded processing with safe IPC serialization of alignment records.
- **Flexible filtering** — BED-based region masking, sliding-window density filtering, and configurable size/quality thresholds.

---

## Installation

### Requirements

- Python ≥ 3.8
- [pysam](https://github.com/pysam-developers/pysam)
- [edlib](https://github.com/Martinsos/edlib) (Python bindings)
- [samtools](https://github.com/samtools/samtools) (for BAM sorting)

### Install from source

```bash
git clone https://github.com/Milia1368/AsmSV.git
cd AsmSV
pip install -e .
```

### Install dependencies only

```bash
pip install pysam edlib
```

---

## Quick Start

### Basic usage

```bash
asmSV \
  --bam1 hap1_to_ref.bam \
  --bam2 hap2_to_ref.bam \
  --ref  reference.fa \
  -o     output.vcf \
  -t     16
```

### With pre-sorted BAMs (skip sorting step)

```bash
asmSV \
  --bam1 hap1.bam --sorted-bam1 hap1.qsort.bam \
  --bam2 hap2.bam --sorted-bam2 hap2.qsort.bam \
  --ref  reference.fa \
  -o     output.vcf
```

### With region masking

```bash
asmSV \
  --bam1 hap1.bam --bam2 hap2.bam \
  --ref reference.fa -o output.vcf \
  --mask repeats.bed
```

### Input preparation

AsmSV expects each haplotype assembly aligned to the reference genome with a long-read aligner such as [minimap2](https://github.com/lh3/minimap2):

```bash
minimap2 -a --eqx -x asm5 reference.fa hap1.fa | samtools sort -o hap1_to_ref.bam
minimap2 -a --eqx -x asm5 reference.fa hap2.fa | samtools sort -o hap2_to_ref.bam
samtools index hap1_to_ref.bam
samtools index hap2_to_ref.bam
```

---

## Method Overview

AsmSV processes each haplotype independently through the following pipeline, then merges results across haplotypes for genotyping.

```
BAM (Hap1/Hap2)
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
  │   ├─ Scenario B: no branches → linear SV (DEL/INS/DUP/CSV)
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

1. **DP Chaining** (`solve_chain_dp`) — Alignment segments are sorted by query coordinate and chained via a penalized longest-path DP (similar to minimap2 chaining). Penalties account for query gaps/overlaps, reference distance, strand flips, and cross-chromosome jumps.

2. **Backbone Identification** (`find_collinear_backbone`) — A second DP pass selects the dominant uni-directional chain. All segments whose strand disagrees with the global majority are forced into branch status. An adaptive overlap filter with dual thresholds (absolute ceiling of 5 kb and ratio cutoff of 50%) prevents alignment artifacts from corrupting the backbone.

3. **SV Extraction** (`resolve_chain_svs`) — Adjacent backbone pairs define gaps. Branch segments within these gaps reveal complex rearrangements:
   - **Scenario B** (no branches): linear SVs — DEL, INS, DUP:TANDEM, or CSV (INS_DEL replacement).
   - **Scenario A** (branches present): single inverted branches trigger local inversion analysis with flanking INDEL decomposition (**[A1]**); all branch–backbone transitions are scanned for nested linear SVs (**[A2]**); branches are condensed into typed clusters and emitted as BND pairs plus span records (**[A3]**).

### Edge Refinement

Each alignment undergoes automatic edge refinement at construction time: CIGAR ends are scanned for stable anchors (≥10 bp match with positive neighborhood score), and unstable flanking regions are trimmed to prevent false-positive INDELs.

---

## Output Format

AsmSV writes standard VCFv4.2. All records include `SVTYPE`, `SVLEN`, `END`, and `GT` fields.

### Structural Variant Types

| SVTYPE | Description |
|--------|-------------|
| DEL | Deletion |
| INS | Insertion (includes short tandem duplications < 1 kb) |
| DUP | Tandem duplication ≥ 1 kb |
| INV | Inversion |
| BND | Breakend / novel adjacency |
| CSV | Complex Sequence Variant (see below) |

### CSV — Complex Sequence Variant

CSV is a custom SVTYPE defined by AsmSV for compound events. The `SUBTYPE` INFO field specifies the mechanism:

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
| LOCAL_RETAIN | In-gap retained sequence (same strand, not a duplication) | Branch [A3] |

### Branch Classification Logic

Branch blocks are classified by distance to the backbone gap and strand orientation:

| Priority | Condition | Same strand | Reverse strand | Cross-chr |
|----------|-----------|-------------|----------------|-----------|
| 0 (highest) | Fully inside gap | LOCAL_RETAIN | INV | — |
| 1 | ≤ 500 bp to gap | DUP:TANDEM | InvDUP:TANDEM | — |
| 2 | ≤ 10 kb to gap | DUP:PROXIMAL | InvDUP:PROXIMAL | — |
| 3 (default) | > 10 kb to gap | DUP:DISP:INTRA | InvDUP:DISP:INTRA | [Inv]DUP:DISP:INTER |

### Key INFO Fields

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

### FILTER Tags

| Tag | Description |
|-----|-------------|
| PASS | High-confidence variant |
| COMPLEX | Multi-cluster (≥2 branch clusters) complex SV summary record |
| LOCAL_BND | Short-range intra-chromosomal BND (not a true translocation) |
| LOCAL_RETAIN | In-gap same-strand branch (retained sequence, not duplication) |
| OVERSIZE | SV exceeds MAX_SV_SIZE threshold |
| EDGE_COMPLEX | Complex rearrangement at read edge |

### Genotype

Phased genotypes are output in the `GT` FORMAT field:

| Genotype | Meaning |
|----------|---------|
| `1\|0` | Variant on Hap1 only |
| `0\|1` | Variant on Hap2 only |
| `1\|1` | Variant on both haplotypes (homozygous) |
| `./.` | Unknown |

---

## Parameters

### Input / Output

| Option | Description |
|--------|-------------|
| `--bam1` | BAM for Haplotype 1 *(required)* |
| `--bam2` | BAM for Haplotype 2 *(required)* |
| `--ref` | Reference genome FASTA *(required)* |
| `-o, --output` | Output VCF path \[default: `./results.vcf`\] |
| `--mask` | BED file of regions to exclude |
| `--sorted-bam1` | Pre-sorted (queryname) BAM for Hap1 |
| `--sorted-bam2` | Pre-sorted (queryname) BAM for Hap2 |

### Variant Filtering

| Option | Default | Description |
|--------|---------|-------------|
| `--min-sv` | 30 | Minimum SV size (bp) |
| `--max-sv` | 1,500,000 | Maximum SV size; -1 = no limit |
| `--min-mapq` | 10 | Minimum mapping quality |
| `--min-csv` | 400 | Minimum CSV (INS_DEL) size (bp) |
| `--max-csv` | 40,000 | Maximum CSV (INS_DEL) size (bp) |
| `--min-tra-size` | 150,000 | Minimum gap for BND/TRA classification |
| `--insdel-tolerance` | 10 | Length tolerance for split-read INS/DEL (bp) |

### Clustering

| Option | Default | Description |
|--------|---------|-------------|
| `--cluster-pos-tol` | 3 | Breakpoint position tolerance (bp) |
| `--cluster-len-tol` | 0 | Length tolerance offset (bp) |
| `--merge-indels-tol` | 3 | Max distance for merging adjacent CIGAR INDELs (bp) |
| `--min-seq-sim` | 1.0 | Minimum sequence similarity for INS clustering |

### Output Format

| Option | Default | Description |
|--------|---------|-------------|
| `--max-alleles-len` | 5,000 | Max sequence length in ALT column |
| `--no-symbolic` | off | Disable symbolic alleles; force full sequences |

### Advanced

| Option | Default | Description |
|--------|---------|-------------|
| `-t, --threads` | 8 | Number of parallel threads |
| `--split-segments` | off | Enable query overlap splitting before DP |
| `--enable-density-filter` | off | Enable sliding-window density noise filter |
| `--density-window` | 1,000 | Density filter window size (bp) |
| `--max-sv-in-window` | 3 | Max SVs per window before suppression |
| `--debug` | off | Enable verbose DEBUG logging |

---

## Project Structure

```
AsmSV/
├── example/
├──src/
  ├──__init__.py  
  ├── asmSV.py        # Main pipeline: parallel dispatch, BAM I/O, serialization
  ├── bam_info.py     # Core: AlignedSegment, DP chaining, backbone, SV extraction
  ├── structures.py   # SVCandidate data structure, BedMaskFilter
  ├── genotype.py     # Cross-haplotype clustering and phased genotyping
  ├── config.py       # Global parameters (SVConfig)
  └── utils.py        # CLI parsing, VCF writing, utility functions
```

---

## Algorithm Details

### DP Chain Scoring

For each segment *i*, the algorithm considers the *K* = 32 nearest predecessors *j*:

```
score[i] = max{ score[j] + AS[i] - penalty(i, j) }
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

### Clustering Rules

| SV Type | Matching Criteria |
|---------|-------------------|
| DEL / DUP / INV | Δpos ≤ 3 bp, Δlen < 2 bp |
| INS | Δpos ≤ 3 bp, Δlen ≤ 0 bp, seq similarity ≥ 1.0 |
| BND / TRA | Both Δpos ≤ 3 bp, same CHR2 |
| CSV (INS_DEL / INV-\*) | Δpos ≤ 3 bp, ΔEND ≤ 3 bp, INS\_LEN + DEL\_LEN within tolerance or 10% |
| CSV (DUP:\* / InvDUP:\*) | Δpos ≤ 3 bp, Δlen < 2 bp |

---

## Examples

#图片

---

## Citation

If you use AsmSV in your research, please cite:

> *Citation information will be added upon publication.*

---

## License

This project is licensed under the [MIT License](https://www.google.com/search?q=LICENSE). to be added.
