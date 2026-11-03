# AsmSV output format

AsmSV writes VCF 4.2 with declared caller-specific INFO fields. The main file, component file, reference INS/DEL file and optional junction file provide different views of the same evidence.

## Variant types

| Type | Meaning |
| --- | --- |
| INS | Inserted sequence |
| DEL | Deleted reference sequence |
| INV | Inverted interval |
| DUP | Additional copy with retained-source support |
| DELINS | Sequence replacement, with insertion and deletion lengths reported separately |
| CPX | Supported composite allele structure |
| BND | Non-reference adjacency endpoint |

Direct tandem copies use `<DUP:TANDEM>`. Other copy annotations may include `DUP_TYPE` and `DUP_ORIENT`. CPX and DELINS are explicitly declared caller extensions. A single structural allele does not establish the number or order of historical mutations. A BND pair represents one adjacency, not necessarily a complete translocation.

## Files

| File suffix | Contents |
| --- | --- |
| Main `.vcf` | Canonical events and independent calls |
| `.components.vcf` | Components associated with canonical events |
| `.indels.vcf` | Reference-relative insertion/deletion projection |
| `.junctions.vcf` | Optional adjacency/BND projection |
| `.events.jsonl` | Supplementary per-contig structures and VCF-record links |
| `.regions.jsonl` | All incomplete contexts, boundaries, members and retained connections |
| `.regions.bed` | Reference interval/port projections of those region annotations |
| `.genotypes.json` | Genotype summary by output view |
| `.run.json` | Completion status, input identities, configuration and output hashes |
| `.preset.json` | Resolved configuration |
| `.log` | Run messages |

The run manifest separates `preset_sha256` (configuration including source identity) from `scoring_sha256` (numerical AQ/DP evaluation). `scoring_policy` is `fixed-common` for defaults or `explicit-json` for an override. Source modes and literal scoring flags remain in the provenance fields; `upstream_scoring_inherited=false` records that they did not supply caller coefficients.

For an output name `calls.vcf`, the component filename is `calls.vcf.components.vcf`. Do not add counts across the VCF views. A deletion within an added copy is not an original-locus deletion; the reference `.indels.vcf` view handles that distinction for INS/DEL comparisons.

## CSV browser and VERIX export

For this research export, CSV is an umbrella set of reported composite alleles: coupled multi-edit/source structures, inverted/dispersed copies with explicit source/target semantics, and local replacements reported separately as REPLACEMENT. It is broader than `SVTYPE=CPX`. Proximity alone does not group independent basic variants, and a composite report does not establish a single historical mutation or a unique decomposition.

The executable selection policy is `asmsv-csv-selection-v1`. Default `--scope csv` includes CPX, valid DELINS, native DUP with `DUP_TYPE=dDUP`, and native inverted DUP (`piDUP_FR/piDUP_RF` or explicit `DUP_ORIENT=INVERTED`). Plain INS/DEL/INV/TRA/BND, ordinary tandem DUP and a single templated INS with unknown source retention are excluded. A generic DUP with only subthreshold filler is not promoted by the exporter. Current INVdel/delINV/delINVdel and inversion/copy-associated edits are already encoded as CPX; inverted copies can remain native DUP. Old unused header labels are not evidence of additional supported production templates.

Run the independent postprocessor after calling:

```bash
asmSV-cpx --vcf calls.vcf --scope csv --out-prefix calls.csv
# From a source checkout: python scripts/export_cpx.py --vcf calls.vcf --out-prefix calls.csv
# Reproduce the previous CPX-only population:
asmSV-cpx --vcf calls.vcf --scope cpx --out-prefix calls.cpx
```

It reads the companion `.components.vcf` and, when present, `.events.jsonl`. Components are required for CPX summaries; standalone complex DUP and DELINS can be exported from their own fields. It does not read `.indels.vcf`, change calls, infer new associations or reconstruct complete ALT sequences. Existing outputs require `--force` to replace. Plain and gzip-compressed inputs are accepted; use `--components` and `--events` for nonstandard companion filenames, or `--no-events` to skip structural evidence. `--event ID` is repeatable; `--pass-only` selects PASS parent/report records.

| Export | Contents |
| --- | --- |
| `calls.csv.jsonl` | One `EXPORT_METADATA` line, then one `CSV_EVENT` line per selected report unit |
| `calls.csv.html` | Self-contained, offline browser with search, REF/ALT schematics and available evidence |
| `calls.csv.verix.vcf` | One original report row per selected composite, with VERIX `single`-format breakpoints |

The JSONL schema is `asmsv-csv-view-v1`. Metadata records `scope` and `selection_policy`. Each event contains `id`, `csv_class`, `main`, `components`, `structures`, `warnings` and `bkps`. Original types and parent/component genotypes remain separate. Each available structure preserves one supporting haplotype's event, context and atoms; different haplotypes are not concatenated into one path. A VCF-only layout is explicitly flagged. HTML embeds the same event data and uses no external network requests. Output prefixes remain user-defined; an old `.cpx` filename does not imply CPX-only selection unless `--scope cpx` is set.

The main HTML view compares REF and ALT using matching segment letters, arrows, retained-source/added-copy groups, insertion blocks and occurrence-specific deletion marks. Dispersed copies show source and target separately. Existing context insertion projections are separate from components. Supported layouts are normalized to reference orientation; uncertain layouts explicitly leave ALT unresolved. VCF-only copy/replacement cartoons describe the reported fields and are not independent biological validation. Diagrams are not to scale and omit small differences. Derived layout models are stored only in HTML.

The export VCF retains the parent's ID, REF/ALT, FILTER, coordinates and sample columns. It adds a self-referencing `EVENT` and these INFO fields:

| Field | Encoding |
| --- | --- |
| `ASMSV_CSV_CLASS` | One controlled structural-family value; present only in this research export, not original caller VCFs |
| `BKPS` | Comma-separated `TYPE:length-chrom:cut[-chrom:cut]` entries for VERIX `single` mode |
| `ASMSV_COMPONENTS` | Comma-separated `id\|type\|chrom\|VCF_POS\|END\|SVLEN\|EDIT_CONTEXT` entries |
| `ASMSV_COMPONENT_GT` | Comma-separated `id\|sample\|GT\|PS\|FT` entries, retaining original component values |
| `ASMSV_EXPORT_WARN` | Missing evidence or shared-component warnings; separate from caller FILTER |

`ASMSV_CSV_CLASS` has `Number=1,Type=String`. It does not replace SVTYPE, CPX_TYPE, CPX_COMPONENTS, DUP_TYPE, DUP_ORIENT or EVENT, and is not a confidence score or mechanism claim. No additional IS_CSV flag is required. Its values are:

| Value | Meaning |
| --- | --- |
| `INV_WITH_EDITS` | Inversion framework with associated edits, including INVdel layouts |
| `DUP_WITH_EDITS` | Added-copy or duplication-junction edits in a copy framework |
| `MULTISOURCE_INSERTION` | Explicit supported multisource insertion framework |
| `INVERTED_DUPLICATION` | Native inverted copy, including piDUP_FR/piDUP_RF |
| `DISPERSED_DUPLICATION` | Native dispersed copy without an explicit inverted annotation |
| `INVERTED_DISPERSED_DUPLICATION` | Native dDUP with inverted orientation |
| `REPLACEMENT` | Local DELINS replacement; no claim of two independent events or unique decomposition |
| `COMPOSITE` | Caller-reported CPX without sufficient fields to safely classify its structural family |

Classification uses explicit frame/source evidence conservatively. A DUP,INS type list alone does not prove whether it is copy editing or multisource insertion. piDUP subtypes imply inversion even when the minimal profile omits redundant DUP_ORIENT. A source interval and target are retained for complex native copies. For example, the caller row `SVTYPE=DUP;DUP_TYPE=dDUP;DUP_ORIENT=INVERTED` remains DUP and gains `ASMSV_CSV_CLASS=INVERTED_DISPERSED_DUPLICATION` only in the export.

Subfields in the component summaries are percent-encoded; split on `|` first, then percent-decode each subfield. Component rows are associated through `EVENT` or `PARENT_EVENT`; a shared component remains associated with every reported parent. The exporter does not regroup records by proximity or CONTEXT_ID.

`BKPS` describes **component-coordinate endpoints**, not directed junctions. CPX_INTERVALS `[first,last]` (1-based inclusive) becomes cuts `first-1,last`; native DUP contributes its original source POS/END and target; DEL uses POS/END; INS uses anchoring POS/SVLEN. DELINS contributes DEL and INS with their separate DEL_LEN/INS_LEN, even when net SVLEN is near zero. The target is `TARGET:0-chrom:cut`. Duplicate entries are removed. Added-copy DEL retains its source-template coordinates and ADDED_DUPLICATE context. No whole-allele INS projection is added; a native DUP's linked INS_ORIGIN=DUPLICATION basis is not counted again as novel INS. Missing required geometry causes an error. The tested parser cannot encode `-` or `:` within contig names.

For two compatible `BKPS` files, a VERIX breakpoint comparison is:

```bash
verix bench -q calls.csv.verix.vcf -t truth.csv.verix.vcf \
  --formats single single --csv_links BKPS BKPS \
  -d 500 -b 0 -s 0 --passonly -o bench_csv_breakpoints
```

This contract was checked against VERIX commit `7cb8524a941ce0d81dd04b0de51e009199386614`. VERIX compares supplied records; it does not detect or define CSVs from assemblies. Truth and query must use the same selection scope and compatible coordinates. Ordinary INS/DEL truth does not become CSV truth through this export. `-s/-S` filter breakpoint intervals, not event length; a dispersed source/target gap can be much larger than the copied segment. Type/genotype/insertion enforcement require compatible semantics. If both inputs have the same normalized class definitions, `--types ASMSV_CSV_CLASS ASMSV_CSV_CLASS --enforce_type` can compare those class strings. It still does not validate copy orientation: ordinary single mode ignores direction. `--breakend_mode` requires suitable BND inputs and cannot be added to this single view. Full endpoint coverage does not establish correct wiring, copy assignment, sequence or mutational history. Report CSV strata separately, especially replacements.

## Main INFO fields

| Field | Description |
| --- | --- |
| SVTYPE | Variant type |
| END | Reference extent where applicable |
| SVLEN | Signed length: negative for DEL, span for INV, net length for composite/replacement calls |
| RNAMES | Supporting assembly contig names, not read depth |
| CPX_COMPONENTS | Component types, such as INV,DEL or DUP,INS |
| CPX_INTERVALS | Associated reference/source intervals |
| CPX_TYPE | Named composite subtype when applicable |
| DUP_TYPE, DUP_ORIENT | Copy layout and orientation |
| TARGET_CHR, TARGET_POS | Copy insertion target |
| SOURCE_CHR, SOURCE_START, SOURCE_END, SOURCE_ORIENT | Insertion-template information where established |
| SOURCE_STATE | Retention state when needed; unknown does not mean deleted |
| INS_ORIGIN | Insertion-origin annotation, including TEMPLATED |
| INS_LEN, DEL_LEN | Both replacement lengths |
| EDIT_INS_LEN, EDIT_DEL_LEN, JUNCTION_INS_LEN | Internal edits or junction insertion length |
| EVENT | Canonical parent ID for a component with one parent |
| PARENT_EVENT | Multiple distinct parent IDs when necessary |
| EDIT_CONTEXT | REFERENCE, ADDED_DUPLICATE or DUPLICATION_JUNCTION |
| MATEID | Paired BND endpoint record |
| REGION_ID | One or more unresolved-region annotation IDs; association does not establish a shared event |

Fields are conditional. The default `--vcf-info minimal` omits redundant labels and normal default states; `--vcf-info full` adds diagnostics. The VCF header is the field dictionary for a particular run.

## Coordinates

VCF POS follows VCF anchoring. SOURCE_START/END and JSON intervals are 0-based half-open; TARGET_POS is a 0-based interbase port; CPX_INTERVALS uses 1-based inclusive intervals. Do not convert all fields with the same offset.

## Genotypes

| Field/value | Interpretation |
| --- | --- |
| GT=1/1 | Both input haplotypes support the same ALT; not genomic copy number |
| GT=1\|0 or 0\|1 with PS | Input haplotype assignment within that phase set |
| GT=0/1 with FT | Dosage estimate with limited counterpart evidence |
| GT=1/. | Counterpart allele unresolved or not directly comparable at a reference locus |
| PS | Phase set associated with compatible physical contig support |
| FT=NO_REFERENCE_SUPPORT | Insufficient counterpart reference traversal |
| FT=AMBIGUOUS_COVERAGE | Competing coverage |
| FT=NONREFERENCE_COMPONENT | Component belongs to an added copy or duplication junction |

Default `--genotyping assembly` retains dosage estimates with evidence filters. `evidence` reports missing counterparts more conservatively. Input haplotype order is not parental origin, and phasing should not be extended beyond PS.

## Filters and partial structures

PASS concerns the variant record; FORMAT/FT concerns genotype evidence. They can differ. PARTIAL denotes an incompletely resolved structure or adjacency. OVERSIZE marks duplication gains meeting the configured size setting.

Unresolved contexts can contain reliable individual calls. Their diagnostic records should not be treated as additional independent CPX events. CRR is a region annotation, not a variant type. CPX_INTERVALS and BND records alone are not a complete reconstructed allele sequence.

## Unresolved-region annotations

Every context remaining PARTIAL or UNRESOLVED is exported by default, including ordinary backbone double gaps with no branch. An output `calls.vcf` produces `calls.vcf.regions.jsonl` and `calls.vcf.regions.bed`, even when the main VCF has no calls. An empty result has metadata/header only. The companion `.events.jsonl` also preserves these regions within each CONTIG bundle.

A region is exactly one incomplete detector context on one physical contig of one input haplotype, after the detector's certified separators. Regions are not merged by distance, reference overlap or graph connected components, and are not merged across contigs or haplotypes. A region may contain zero, one or several reported variants; membership never changes atom ownership, variant counts, FILTER or GT. Region records have `counts_as_variant=false` and `genotype=null`.

The standalone JSONL starts with `record_kind=REGION_METADATA`, `schema=asmsv-regions-v1` and `policy=asmsv-unresolved-regions-v1`, followed by `record_kind=UNRESOLVED_REGION` records:

| Field | Meaning |
| --- | --- |
| region_id, context_id | Stable annotation ID and original detector context ID |
| haplotype, query_name | Physical supporting assembly contig; input haplotype is not parental origin |
| resolution, reasons | PARTIAL/UNRESOLVED and specific original failure reasons |
| labels | Always UNRESOLVED_REGION; optional observational CRR flag |
| reason_categories | ASSEMBLY_GAP, REFERENCE_GAP, COMPUTE_LIMIT, OPEN_BOUNDARY, ALIGNMENT_AMBIGUITY, SEQUENCE_AMBIGUITY, SOURCE_UNCERTAINTY or STRUCTURE_UNRESOLVED |
| query_span | Extent on the forward assembly query; includes an existing framework span when relevant |
| boundaries.left/right | Backbone block, query cut, reference port, exact-anchor support and MAPPED/UNMAPPED/OPEN state |
| boundary_order | FORWARD or OVERLAP for the outer query cuts; an overlap is not silently normalized into a closed interval |
| reference_intervals | Separate TARGET, BOUNDARY_PORT, MEMBER and SOURCE intervals; spans describe evidence, not a claim that every base is altered |
| members | Branch block indices, internal atom/event IDs and retained junction IDs |
| observed_path | Query-ordered alignment occurrences, including support anchors, with distinct occurrence IDs and entry/exit ports |
| junctions | Retained template connections and independently certified selected-path connections |
| unretained_transitions | Transitions not certified, with the explicit reason; absence of an edge does not mean reference structure |
| separators, sequence_checks | Existing decomposition and residual sequence-check evidence |
| canonical_event_ids | Associated main-parent IDs resolved through the event alias table; internal atom IDs are not converted into fictitious VCF IDs |
| vcf_records | Exact `(view, record_id)` links to records actually written, including optional junction records |

`members.block_indices` contains branch members; backbone support blocks are separately identified in boundaries and observed_path. Internal atom/event IDs refer to the companion CONTIG bundle. The VCF ID index is authoritative for published records. Main, component, indel and junction records retain `REGION_ID` in both minimal and full profiles, including all associated region IDs after cross-haplotype clustering. Membership can overlap and implies association, not one biological event.

All JSON intervals are 0-based half-open; reference ports are exact interbase cuts and can have `start=end`. A TARGET interval is bounded by two mapped backbone ports on the same chromosome. Otherwise each mapped side is exported separately as BOUNDARY_PORT. Source/member intervals stay separate across loci and chromosomes; an open side does not invent a missing reference coordinate. Query extent reaches query position 0 or contig length on an open side. MAPPED means a CIGAR projection exists; `exact_anchor` separately states whether the boundary has the required exact support. These are observable context boundaries, not confidence intervals for a unique historical event.

The BED file is BED6 plus `interval_role,haplotype,query_name,resolution,labels,reasons`. Query names are percent-encoded. It has one row per reference projection, so several rows can share a region ID. Zero-length ports are projected to `[cut,cut+1)` for viewing, or the last base at chromosome end; exact cuts remain in JSON. Do not sum BED rows as events or treat overlapping role intervals as disjoint affected sequence.

New selected-path connections require mapped ports and an exact terminal anchor of the existing `min_anchor` length on both sides. Conflicting query overlap is rejected; exact shared-query overlap can be normalized. Junction sequence must be ACGT and within the existing seed-length budget. A same-axis transition with SV-sized residuals on both query and reference is left unresolved. Normal collinear steps do not become BNDs. Each rejected transition remains visible with its reason. Repeated visits to the same reference ports retain separate physical occurrences in JSON.

Connections with `evidence=SELECTED_PATH_EXACT_PORTS,adjacency_status=OBSERVED` passed these local checks. Existing template connections without that additional certificate are labeled `TEMPLATE_RECORDED`. Neither label evaluates alternative mappings or assembly correctness. Known connecting sequence is preserved in query order in JSON; the optional BND view includes certified connecting sequence in ALT. BNDs from incomplete contexts retain FILTER=PARTIAL, separately from whether an individual connection was observed. Shared port geometry with different known junction sequences is not collapsed.

CRR is optional: an incomplete context must contain at least two locally certified noncollinear connections and consecutive member occurrences with incompatible chromosome, orientation or wholly reversed reference order. Fragmentation of a single monotone inversion source, an open boundary, an assembly/reference gap, budget exhaustion or seed ambiguity alone is insufficient. This is a conservative observation flag, not a CSV definition, confidence score, event type or complete rearrangement classification. Regions without CRR are still exported.

Region diagnostics are retained when `--mask` removes variant candidates; metadata declares `region_mask_policy=RETAIN_ALL_DIAGNOSTICS`. `vcf_records` lists only surviving records. The run manifest records region counts, reasons, filenames and hashes. Region localization, connection recovery and complete-event resolution must be evaluated separately. All annotations have SELECTED_PATH scope: rejected or DP-discarded alignments are not recovered by this layer.

The base format is defined by the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf). AsmSV-specific fields and interpretations are declared in the output header.
