# AsmSV user guide

This checkout describes local **1.4.0.dev0**. Published stable 1.3.0 has different calling behavior; retain the version, source hashes and input provenance when comparing runs. The run manifest records `software_version`, and VCF headers record `asmsvVersion`.

## Inputs

AsmSV requires `--bam1`, `--bam2`, `--ref` and `--aligner`. Both BAMs must describe assemblies aligned to the same reference. Reference identifiers and lengths must match the FASTA dictionary.

Supported source declarations are `minimap2`, `winnowmap` and `vacmap` (`winnowmap2` is an alias). HG002–CHM13 reference INS/DEL evaluation covers minimap2 asm5/asm20, Winnowmap asm20 and VACmap asm inputs. Alignment evidence and accuracy still vary across these inputs.

Retain primary and supplementary alignments with real CIGARs. EQX, MD or cs provide alignment-difference information; valid NM tags can also be used when supported. Text `CG:Z` replacement CIGARs are unsupported. For hard-clipped inputs without a complete BAM sequence, supply indexed original assemblies through `--query-fasta1/2`.

## Basic command

```bash
asmSV --aligner minimap2 \
  --bam1 hap1.bam --bam2 hap2.bam \
  --ref reference.fa --sample SAMPLE \
  -t 4 -o calls.vcf
```

The source mode is inferred from consistent BAM program headers. Use `--aligner-preset` when an explicit declaration is needed. Conflicting source declarations are reported. Source modes describe the BAM's origin; they do not select numerical caller parameters.

All supported sources use the fixed `common-assembly-hq-v5` template: AQ mismatch 19, gap open (39,81), gap extend (3,1), structural threshold/scale 50/10; DP query gap 0.25, excess overlap 0.5, chromosome change 500 and inversion 200. AQ still scores each alignment's actual differences and applies MAPQ policy. Native AS and upstream `-A/-B/-O/-E` flags do not supply caller coefficients. Different scoring flags between the two BAMs are allowed and recorded.

For an explicit experiment, pass `--preset-file presets/common-assembly-hq-v5.json` or a custom schema-5 JSON. Partial JSON overrides inherit the fixed common template. Source-specific defaults from earlier builds must be supplied explicitly to reproduce their numerical scoring; selecting `--aligner-preset asm20` alone now only declares provenance. The supplied example is the same template used without a preset file. High-divergence or error-rich assemblies require separate calibration; the HG002 evaluation does not establish an optimum for those inputs.

## Common options

| Option | Default | Description |
| --- | --- | --- |
| `--bam1`, `--bam2` | Required | Haplotype assembly alignments |
| `--ref` | Required | Indexed reference FASTA |
| `--aligner` | Required | Source aligner |
| `--aligner-preset` | BAM header or declared-tool mode | Source mode for provenance and consistency checks |
| `--preset-file` | Fixed common template | Explicit schema-5 AQ/DP override |
| `--sample` | SAMPLE | Output VCF sample name |
| `-o`, `--output` | results.vcf | Main output path |
| `-t`, `--threads` | 4 | Workers per haplotype |
| `--min-sv` | 30 | Minimum variant-size setting, in bp |
| --min-inv / --min-inv-size | follows --min-sv | Minimum standalone INV span; independent of INS/DEL and CPX inversion frames |
| --max-indel-size | -1 | Inclusive INS/DEL reporting maximum; larger records retain type/sequence with OVERSIZE |
| --max-ref-gap | 150,000 | Reference-gap bound for heuristic fallback; not a universal cap or TRA classification threshold |
| `--min-mapq` | 10 | Minimum alignment MAPQ |
| `--query-fasta1/2` | None | Indexed original assembly FASTAs |
| `--sorted-bam1/2` | Automatic | Existing name-sorted BAMs |
| `--work-dir` | Beside each input BAM | Override the reusable sorting cache directory |
| `--scratch-dir` | System temporary directory | BAM sort spill files and temporary structural evidence |
| `--genotyping` | assembly | assembly, evidence or presence |
| `--vcf-info` | minimal | minimal or full annotation profile |
| `--junction-vcf` | Disabled | Write a separate BND projection, including retained connections from incomplete regions |
| `--mask` | None | Exclude candidates overlapping a reference BED |
| `--max-sv` | 1,500,000 | Duplication-gain OVERSIZE setting; -1 disables it; not a universal length filter |

Use `asmSV --help` for all options. Changing the minimum size can affect detection and grouping; it is not equivalent to applying a post-calling VCF size filter.

## Reuse sorted BAMs

Without `--sorted-bam1/2` or `--work-dir`, each input gets a cache in its own directory:

```text
/data/hap1.bam  -> /data/hap1.asmsv.qname.bam
/other/hap2.bam -> /other/hap2.asmsv.qname.bam
```

Each cache has a `.bam.json` receipt recording the input's resolved path, size and modification time, and the cached BAM's identity. Before reuse, AsmSV checks this receipt, BAM integrity with quickcheck, and its queryname sort header. A missing/stale receipt, changed input or invalid cache triggers regeneration. Complete cache files are published atomically; original inputs are not modified. Automatically managed `*.asmsv.qname.bam` and their receipts may be replaced when validation fails.

Changing the VCF output directory reuses the same input-local caches. For an input directory that is not writable, or to place caches on a local disk, set `--work-dir /path/to/writable/cache`. An explicit work directory overrides input-local caches and retains the previous hash-based cache layout, including reuse of valid older receipts. Version 1.3.0 output-local `cache/` files remain available via `--work-dir /previous/output/cache`.

You can also prepare grouped BAMs yourself. Explicit `--sorted-bam1/2` takes precedence over automatic caching for that haplotype and leaves supplied files unchanged:

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

Name sorting uses the system temporary directory (or `--scratch-dir`) for spill files and publishes each completed BAM beside its input, or in an explicit `--work-dir`. Prefer a local disk with sufficient free space for temporary files. Input-local caches on network storage can still be slower to create; use a local `--work-dir` when appropriate. A first run includes both BAM sorts; a cache hit skips them. Cached-input timings therefore do not describe a first run.

AsmSV logs each preparation/calling stage, cache hits and a heartbeat every 30 seconds while a stage is running. The `.run.json` records the current stage, per-stage wall times, total elapsed seconds, and the pysam/samtools versions and libdeflate availability. A long sort does not mean variant detection has started.

BAM compression speed also depends on how pysam was built. The startup log reports `libdeflate=False` when that acceleration is unavailable. Check the active interpreter with:

```bash
python -c 'import pysam, pysam.config; print(pysam.__version__, pysam.config.HAVE_LIBDEFLATE)'
```

For a conda installation, install a libdeflate-enabled pysam before installing AsmSV. One validated version is 0.23.3; the version alone does not guarantee build features. Full environment creation and installation steps are in the [README](../README.md#conda-installation-recommended):

```bash
conda install --override-channels -c defaults pysam=0.23.3 pip
python -m pip install .
```

## Upgrading from 1.2.x

The distribution name remains AsmSV and the executable remains `asmSV`. This version requires Python 3.10 or newer. Declare the source aligner explicitly and review `asmSV --help` before reusing older commands.

The current output uses CPX and DELINS annotations and separate event/component/reference-change views. Older CSV subtype names and script assumptions should not be applied without checking the [output specification](OUTPUT_FORMAT.md). Historical development-only `--detector` and `--min-csv` switches are not supported.

## Troubleshooting

| Symptom | Check |
| --- | --- |
| BAM/FASTA dictionary mismatch | Reference names, sequence lengths and reference version |
| Missing complete query sequence | Supply indexed matching original assembly FASTAs |
| Query-name groups are not contiguous | Recreate the name-sorted BAM cache |
| Cannot create a sorting cache beside the input | Supply a writable `--work-dir` or existing `--sorted-bam1/2` |
| Source/preset conflict | BAM `@PG` commands and `--aligner` / `--aligner-preset` |
| Empty or partial result | Log, MAPQ, assembly coverage, run completion and partial-call annotations |
| Unexpected genotype | FORMAT/FT, phase sets and completeness of both assemblies |
| Excess memory | Worker count, input fragmentation, temporary storage and sequence transport |
| Slow first run / long BAM sort | Stage log, libdeflate-enabled pysam, local `--scratch-dir`, and reusable `--work-dir`; compare timings with the same cache conditions |

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


## Complete insertion sequences and source links

Current1.3.0 builds retain all available INS payloads in ALT and .sequences.fa with length/SHA256 links. Inspect .indels.vcf for reference changes and .junctions.vcf (enabled by --junction-vcf) for paired source/target topology. RELOCATION_ID identifies same-haplotype sequence-compatible DEL/INS links; inspect the corresponding RELOCATION record in .events.jsonl for source-search scope, original atoms and exact boundary-equivalence proof.

Use --tandem-target-indels only when you need the separate .target_indels.vcf representation at declared tandem-copy targets. Default reference indels and this view describe equivalent alleles and must not be pooled. Declare the chosen view in a benchmark; a distance window can distinguish equivalent tandem boundaries. See [output semantics](OUTPUT_FORMAT.md).


## Independent inversion and indel size policies

Use --min-sv 50 --min-inv 30 to retain supported standalone inversions down to 30 bp while keeping the INS/DEL minimum at 50 bp. A lower reporting threshold does not create missing alignment evidence; the pure exact-gap model also needs its intrinsic 16-bp core and usable flanks.

Known INS payloads remain complete literal ALT at every length, including OVERSIZE records. --max-indel-size 10000 marks INS/DEL records longer than 10,000 bp; it does not convert them to TRA. --max-ref-gap controls only heuristic fallback, so a directly observed CIGAR deletion or a certified pure deletion can be larger. Duplication --max-sv remains separate. INS/DEL projections of a DUP follow the indel report cap; the DUP parent and BND view have their own semantics and must not be summed as independent events.
