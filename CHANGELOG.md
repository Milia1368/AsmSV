# Release changes

## 1.4.0.dev0 — development branch update, 2026-11-07

The historical `release/asmsv-v1.3.0` branch now contains this development version. The previous 1.3.0 snapshot remains available at commit `eb64bc2`. The package and CLI report `1.4.0.dev0`; this branch update does not create a version tag or a GitHub Release.

- Use complete EQX CIGAR directly; resolve M with synchronized MD/cs, and use NM only as an aggregate fallback.
- Use a continuous, capped-gap evaluation length for AQ density while retaining DP rewards and the common coefficients.
- Recover pure inversions with exact forward flanks before general sequence-search budgets. Optimize the same bounded core model without changing accepted alleles.
- Keep all known INS sequence in ALT and content-addressed FASTA, independent of length.
- Add independent --min-inv, optional --max-indel-size and heuristic --max-ref-gap policies. Size never relabels a local INS/DEL as a translocation.
- Retain certified source/target joins and scoped observed-DEL/INS relations without asserting unique genome-wide origin.
- Keep bounded allele normalization and the explicitly optional tandem-target output; evaluate representation differences separately from new detection.
- Remove the unused full-sequence NN verifier, its dormant scoring-view API and unused VCF sequence helpers. Historical source/tests remain in local research archives.

Scores, output interpretation and inferred calls can differ from stable 1.3.0. Record version, source hashes, input alignment provenance and evaluation protocol.


## 1.3.0 maintenance update — 2026-11-03

Supported assembly alignment inputs now use the same `common-assembly-hq-v5` AQ/DP scoring template. Minimap2 asm5/asm10/asm20, Winnowmap and VACmap source declarations no longer select different caller coefficients. Native AS and upstream scoring flags are recorded as provenance and do not supply caller scores.

`--aligner` is still required. `--aligner-preset` still declares and validates the upstream mode. `--preset-file` remains available for explicit schema-5 AQ/DP overrides. Its omission uses the fixed common template; changing only `--aligner-preset` does not tune detection. Source-specific JSON examples have been replaced by the common example.

The run manifest records the numerical `scoring_sha256`, the source-inclusive `preset_sha256`, and whether scoring came from the fixed template or explicit JSON. Record these fields and the Git commit when comparing experiments.

Incomplete contexts retain region boundaries, path members, observed connections and reasons in `.regions.jsonl` and `.regions.bed`. `REGION_ID` connects existing VCF records to these observations. Regions and the optional CRR label do not count as additional variants or complete reconstructed events. Supported partial INV and BND records retain their partial status.

`asmSV-cpx` exports a CSV browser and a VERIX-compatible breakpoint view without rerunning detection. It preserves caller genotypes and distinguishes export categories from the original SVTYPE.

The production structural defaults remain unchanged, including GAP_TOL=100 bp, INV boundary cost=200 and chromosome-change cost=500. This release does not add new rescue rules for dense inversions. Retained sorting caches and preparation/calling progress reports remain available.

### Compatibility

Earlier 1.3.0 builds could derive AQ costs from the source aligner, mode and command-line scoring flags. Reproduce those runs with their frozen source and saved configuration, or explicitly supply the earlier coefficients using `--preset-file`. Default output for those sources can change under this maintenance update. Region and browser annotations do not establish complete-event accuracy or eliminate known dense-rearrangement limitations.
