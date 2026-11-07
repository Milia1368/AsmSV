# AsmSV

```text
\ \/ /    █████╗ ███████╗███╗   ███╗███████╗██╗   ██╗
 \  /    ██╔══██╗██╔════╝████╗ ████║██╔════╝██║   ██║
 /  \    ███████║███████╗██╔████╔██║███████╗██║   ██║
/ /\ \   ██╔══██║╚════██║██║╚██╔╝██║╚════██║╚██╗ ██╔╝
        ██║  ██║███████║██║ ╚═╝ ██║███████║ ╚████╔╝
        ╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝╚══════╝  ╚═══╝
```

**AsmSV: Assembly-based Haplotype-Resolved Structural Variant Caller**

[Installation](#installation) · [Quick start](#quick-start) · [User guide](docs/USAGE.md) · [Output format](docs/OUTPUT_FORMAT.md) · [CSV examples](docs/CSV_EXAMPLES_EN.md) / [中文](docs/CSV_EXAMPLES_ZH.md) · [Citation](#citation)

## Introduction

AsmSV detects structural variants from two haplotype-resolved genome assemblies aligned to the same reference. It takes one BAM per haplotype and writes VCF 4.2 records with genotype and structural annotations.

This branch provides **1.4.0.dev0**, a development version with CIGAR-first evidence scoring, exact flanked inversion recovery, complete known insertion payloads and independent size-reporting policies. Its historical branch name remains **`release/asmsv-v1.3.0`**; the previous 1.3.0 snapshot is commit `eb64bc2`. The caller uses one common AQ/DP scoring template across supported assembly aligners and retains incomplete-region evidence. BAM sorting caches remain reusable between runs. The package name remains **AsmSV**, and the command remains **`asmSV`**. See [release changes and compatibility](CHANGELOG.md).

### Features

- INS, DEL, INV and DUP detection, with DELINS and supported composite structures represented explicitly.
- Main events, linked components and a reference-relative INS/DEL output view.
- Unresolved-region annotations with explicit boundaries, members and available connections.
- Genotype evidence filters and phase-set annotations where supported by the assemblies.
- Configurable inputs, multithreaded execution and reusable BAM sorting caches.

## Installation

Requirements: Python **3.10+**, `pysam>=0.19` and `edlib>=1.3.9`. Linux is the tested platform. A pysam build with **libdeflate** is recommended for fast BAM compression.

### Install the current branch

Create an environment with a libdeflate-enabled pysam before installing AsmSV. The following Linux setup uses the validated `defaults` build; pip installs edlib with AsmSV:

```bash
conda create -y -n asmsv_v1.4_dev --override-channels \
  -c defaults python=3.10 pysam=0.23.3 pip
conda activate asmsv_v1.4_dev
git clone --branch release/asmsv-v1.3.0 https://github.com/Milia1368/AsmSV.git
cd AsmSV
python -m pip install .
python -c 'import pysam, pysam.config; print("pysam", pysam.__version__, "libdeflate", pysam.config.HAVE_LIBDEFLATE)'
asmSV --version
asmSV --help
```

The dependency check should report `libdeflate 1` and `asmSV --version` should report `AsmSV 1.4.0.dev0`. The startup log also reports libdeflate availability. Builds without it can sort large BAMs substantially more slowly.

### Python virtual environment

```bash
git clone --branch release/asmsv-v1.3.0 https://github.com/Milia1368/AsmSV.git
cd AsmSV
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
asmSV --help
```

Python dependencies are installed automatically; compression performance depends on the pysam build. See [temporary storage and installation checks](docs/USAGE.md#memory-and-temporary-files).

The GitHub commands above select the historical release branch, now containing 1.4.0.dev0. For an editable installation, use `python -m pip install -e .`. The source-checkout launcher is `bash run.sh ...`.

### Update an existing installation

In a checkout with its runtime dependencies already installed, select the release branch and reinstall the local source while retaining the conda dependency builds. Commit or otherwise preserve any local edits before switching branches:

```bash
git fetch origin
git switch release/asmsv-v1.3.0
git pull --ff-only
python -m pip install --force-reinstall --no-deps .
asmSV --version
git rev-parse HEAD
```

Record the Git commit as well as the version when reproducing results. The branch name retains 1.3.0 for continuity, while the package and CLI report 1.4.0.dev0. Use the Git commit and saved scoring configuration to distinguish updates; see [compatibility notes](CHANGELOG.md) before comparing with earlier builds.


## Quick start

Prepare two assembly-to-reference BAMs and an indexed reference FASTA. Preserve supplementary alignments, real CIGARs and alignment-difference tags. For example, using minimap2 and samtools:

```bash
samtools faidx reference.fa
minimap2 -ax asm5 --eqx -t 16 reference.fa hap1.fa | samtools view -b -o hap1.bam
minimap2 -ax asm5 --eqx -t 16 reference.fa hap2.fa | samtools view -b -o hap2.bam
```

Call variants:

```bash
asmSV \
  --aligner minimap2 \
  --bam1 hap1.bam --bam2 hap2.bam \
  --ref reference.fa --sample SAMPLE \
  -t 4 -o calls.vcf
```

AsmSV uses one fixed AQ/DP scoring template for the supported assembly aligners. `--aligner`, BAM source modes and upstream scoring flags retain their provenance role; they do not change caller scores. An explicit `--preset-file` can override the common template for controlled experiments.

If BAM hard clipping prevents recovery of complete contig sequences, add `--query-fasta1 hap1.fa --query-fasta2 hap2.fa` with indexed original assemblies. BAM and reference sequence names and lengths must agree. By default, the caller creates `hap1.asmsv.qname.bam` beside `hap1.bam` and the corresponding cache beside `hap2.bam`. Each run validates the cache and its source receipt before reuse; changing the output directory does not trigger another sort. Use `--work-dir` to choose a different cache directory, or `--sorted-bam1/2` to supply existing grouped BAMs.

See the [user guide](docs/USAGE.md) for parameters, input checks, examples and troubleshooting.

## Output

| File | Purpose |
| --- | --- |
| `calls.vcf` | Main events and independent variants |
| `calls.vcf.components.vcf` | Linked event components |
| `calls.vcf.indels.vcf` | Reference-relative INS/DEL view |
| `calls.vcf.events.jsonl` | Additional structural and provenance information |
| `calls.vcf.regions.jsonl`, `.regions.bed` | Incomplete-region annotations, members, boundaries and connections |
| `calls.vcf.genotypes.json` | Genotype summary |
| `calls.vcf.run.json`, `.preset.json`, `.log` | Run status, configuration and log |
| `calls.vcf.junctions.vcf` | Optional BND projection with `--junction-vcf` |

Use `.indels.vcf` when comparing with an INS/DEL truth set. These files provide different views of the same calls; do not concatenate them and count every record independently. Field definitions and coordinate conventions are in the [output specification](docs/OUTPUT_FORMAT.md).

Incomplete structures are retained as region annotations even when no complete event can be called. `REGION_ID` links existing VCF records to these annotations. Region membership and the optional CRR observation label do not add variants or establish one shared mutation. See [unresolved regions](docs/OUTPUT_FORMAT.md#unresolved-region-annotations).

To browse composite SVs and export a VERIX breakpoint view after calling:

```bash
asmSV-cpx --vcf calls.vcf --scope csv --out-prefix calls.csv
```

Open `calls.csv.html` in a browser. The tool also writes `calls.csv.jsonl` and `calls.csv.verix.vcf`, preserving original SVTYPE and genotypes. Default scope includes CPX, inverted/dispersed DUP and local DELINS replacements, classified by the export-only ASMSV_CSV_CLASS field. Use `--scope cpx` for the previous CPX-only population. It does not rerun detection. See the [selection rules, format and VERIX usage](docs/OUTPUT_FORMAT.md#csv-browser-and-verix-export).

Worked CSV examples with a small real-output subset: [English](docs/CSV_EXAMPLES_EN.md), [中文](docs/CSV_EXAMPLES_ZH.md). They cover parent/component links, copy-local edits, zero-net replacements and incomplete regions.

## Installation example

```bash
python examples/make_toy_inputs.py --output toy
asmSV --aligner minimap2 \
  --bam1 toy/hap1.bam --bam2 toy/hap2.bam \
  --ref toy/reference.fa -t 2 -o toy/calls.vcf
```

The example generates synthetic inputs with a 200 bp deletion and a 60 bp insertion on Hap1. It is an installation check, not a biological accuracy benchmark.

## Scope

AsmSV reports structures supported by the input assemblies and alignments. Assembly gaps, repeats and incomplete alignments can leave structures partially resolved. CPX calls do not establish mutation chronology, and BND records do not necessarily represent complete translocations. Review genotype filters and partial-call annotations before downstream analysis.

## Installation checks

CI checks package installation, both command-line entry points and the synthetic calling example above. Development regression tests and evaluation scripts are maintained locally and are excluded from the public repository and source distributions.

## Citation

Citation information will be added upon publication.

## License and support

AsmSV is distributed under the [MIT license](LICENSE). Report problems through [GitHub Issues](https://github.com/Milia1368/AsmSV/issues), including the version, command, relevant input header and a small reproducing example when possible.
