# AsmSV user guide

## Inputs

AsmSV requires `--bam1`, `--bam2`, `--ref` and `--aligner`. Both BAMs must describe assemblies aligned to the same reference. Reference identifiers and lengths must match the FASTA dictionary.

Supported source declarations are `minimap2`, `winnowmap` and `vacmap`. Input handling is tested for all three; whole-genome validation currently uses minimap2 asm5. Do not assume identical performance across aligners or presets.

Retain primary and supplementary alignments with real CIGARs. EQX, MD or cs provide alignment-difference information; valid NM tags can also be used when supported. Text `CG:Z` replacement CIGARs are unsupported. For hard-clipped inputs without a complete BAM sequence, supply indexed original assemblies through `--query-fasta1/2`.

## Basic command

```bash
asmSV --aligner minimap2 \
  --bam1 hap1.bam --bam2 hap2.bam \
  --ref reference.fa --sample SAMPLE \
  -t 4 -o calls.vcf
```

The source mode is inferred from consistent BAM program headers. Use `--aligner-preset` when an explicit declaration is needed. Conflicting source declarations are reported.

## Common options

| Option | Default | Description |
| --- | --- | --- |
| `--bam1`, `--bam2` | Required | Haplotype assembly alignments |
| `--ref` | Required | Indexed reference FASTA |
| `--aligner` | Required | Source aligner |
| `--sample` | SAMPLE | Output VCF sample name |
| `-o`, `--output` | results.vcf | Main output path |
| `-t`, `--threads` | 4 | Workers per haplotype |
| `--min-sv` | 30 | Minimum variant-size setting, in bp |
| `--min-mapq` | 10 | Minimum alignment MAPQ |
| `--query-fasta1/2` | None | Indexed original assembly FASTAs |
| `--sorted-bam1/2` | Automatic | Existing name-sorted BAMs |
| `--work-dir` | Output-local cache | Reusable sorting cache |
| `--genotyping` | assembly | assembly, evidence or presence |
| `--vcf-info` | minimal | minimal or full annotation profile |
| `--junction-vcf` | Disabled | Write a separate BND projection |
| `--mask` | None | Exclude candidates overlapping a reference BED |
| `--max-sv` | 1,500,000 | Duplication-gain OVERSIZE setting; -1 disables it; not a universal length filter |

Use `asmSV --help` for all options. Changing the minimum size can affect detection and grouping; it is not equivalent to applying a post-calling VCF size filter.

## Reuse sorted BAMs

```bash
samtools sort -n -@ 4 -o hap1.qname.bam hap1.bam
samtools sort -n -@ 4 -o hap2.qname.bam hap2.bam
asmSV --aligner minimap2 \
  --bam1 hap1.bam --bam2 hap2.bam \
  --sorted-bam1 hap1.qname.bam --sorted-bam2 hap2.qname.bam \
  --ref reference.fa -t 4 -o calls.vcf
```

Supplied query-name groups must be contiguous. Coordinate-sorted BAM indexes are not required by the caller.

## Source annotations and optional realignment

`--dup-source-search` accepts `off`, `adjacent` (default) or `local`. Local search additionally uses a bounded same-chromosome region controlled by `--dup-source-radius`. A source annotation alone does not prove retained duplication. Wider searches can take longer and may expose additional ambiguity.

Optional replacement realignment is disabled by default:

```bash
python -m pip install '.[realign]'
asmSV --aligner minimap2 --bam1 hap1.bam --bam2 hap2.bam \
  --ref reference.fa --replacement-realign fallback -o calls.vcf
```

The optional feature has execution limits. Failure to resolve a structure is retained in diagnostic output. Existing `--copy-source-*` spellings remain aliases for `--dup-source-*`.

## Memory and temporary files

Increasing the worker count can increase peak memory. Use `--work-dir` for reusable sorting storage and `--scratch-dir` for temporary files. `--sequence-transport copy` selects process copies instead of automatic shared-memory transport. Original assembly FASTAs can provide sequence access without relying on a complete BAM SEQ.

## Upgrading from 1.2.x

The distribution name remains AsmSV and the executable remains `asmSV`. This version requires Python 3.10 or newer. Declare the source aligner explicitly and review `asmSV --help` before reusing older commands.

The current output uses CPX and DELINS annotations and separate event/component/reference-change views. Older CSV subtype names and script assumptions should not be applied without checking the [output specification](OUTPUT_FORMAT.md). Historical development-only `--detector` and `--min-csv` switches are not supported.

## Troubleshooting

| Symptom | Check |
| --- | --- |
| BAM/FASTA dictionary mismatch | Reference names, sequence lengths and reference version |
| Missing complete query sequence | Supply indexed matching original assembly FASTAs |
| Query-name groups are not contiguous | Recreate the name-sorted BAM cache |
| Source/preset conflict | BAM `@PG` commands and `--aligner` / `--aligner-preset` |
| Empty or partial result | Log, MAPQ, assembly coverage, run completion and partial-call annotations |
| Unexpected genotype | FORMAT/FT, phase sets and completeness of both assemblies |
| Excess memory | Worker count, input fragmentation, temporary storage and sequence transport |

`--mask` removes calls; it is not a repeat-region annotation option. Avoid applying an exclusion BED when the intention is only to label regions.

## Optional evaluation

```bash
python -m pip install '.[benchmark]'
python scripts/evaluate_benchmark.py calls.vcf.indels.vcf \
  --ref reference.fa --truth truth.vcf.gz --bed benchmark.bed
python scripts/benchmark_genotypes.py calls.vcf.indels.vcf \
  --output genotype_benchmark \
  --ref reference.fa --truth truth.vcf.gz --bed benchmark.bed
```

The scripts require bcftools on PATH. Use fresh benchmark output directories and consistent reference/truth/region definitions. The benchmark extra pins Truvari to the evaluation-script version. Detection and genotype matching policies differ; inspect the saved commands before comparing metrics.
