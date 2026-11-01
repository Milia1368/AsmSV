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

[Installation](#installation) · [Quick start](#quick-start) · [User guide](docs/USAGE.md) · [Output format](docs/OUTPUT_FORMAT.md) · [Citation](#citation)

## Introduction

AsmSV detects structural variants from two haplotype-resolved genome assemblies aligned to the same reference. It takes one BAM per haplotype and writes VCF 4.2 records with genotype and structural annotations.

Version **1.3.0** updates the implementation and output annotations of the existing AsmSV project. The package name remains **AsmSV**, and the command remains **`asmSV`**.

### Features

- INS, DEL, INV and DUP detection, with DELINS and supported composite structures represented explicitly.
- Main events, linked components and a reference-relative INS/DEL output view.
- Genotype evidence filters and phase-set annotations where supported by the assemblies.
- Configurable inputs, multithreaded execution and reusable BAM sorting caches.

## Installation

Requirements: Python **3.10+**, `pysam>=0.19` and `edlib>=1.3.9`. Linux is the tested platform. The Python dependencies are installed automatically.

```bash
git clone --branch release/asmsv-v1.3.0 https://github.com/Milia1368/AsmSV.git
cd AsmSV
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
asmSV --help
```

For an editable installation, use `python -m pip install -e .`. The source-checkout launcher is `bash run.sh ...`.

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

If BAM hard clipping prevents recovery of complete contig sequences, add `--query-fasta1 hap1.fa --query-fasta2 hap2.fa` with indexed original assemblies. BAM and reference sequence names and lengths must agree. The caller creates reusable name-sorted BAM caches automatically.

See the [user guide](docs/USAGE.md) for parameters, input checks, examples and troubleshooting.

## Output

| File | Purpose |
| --- | --- |
| `calls.vcf` | Main events and independent variants |
| `calls.vcf.components.vcf` | Linked event components |
| `calls.vcf.indels.vcf` | Reference-relative INS/DEL view |
| `calls.vcf.events.jsonl` | Additional structural and provenance information |
| `calls.vcf.genotypes.json` | Genotype summary |
| `calls.vcf.run.json`, `.preset.json`, `.log` | Run status, configuration and log |
| `calls.vcf.junctions.vcf` | Optional BND projection with `--junction-vcf` |

Use `.indels.vcf` when comparing with an INS/DEL truth set. These files provide different views of the same calls; do not concatenate them and count every record independently. Field definitions and coordinate conventions are in the [output specification](docs/OUTPUT_FORMAT.md).

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

## Tests

```bash
python -B -m unittest discover -s tests -v
```

## Citation

Citation information will be added upon publication.

## License and support

AsmSV is distributed under the [MIT license](LICENSE). Report problems through [GitHub Issues](https://github.com/Milia1368/AsmSV/issues), including the version, command, relevant input header and a small reproducing example when possible.
