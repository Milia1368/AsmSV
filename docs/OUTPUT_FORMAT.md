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
| `.genotypes.json` | Genotype summary by output view |
| `.run.json` | Completion status, input identities, configuration and output hashes |
| `.preset.json` | Resolved configuration |
| `.log` | Run messages |

For an output name `calls.vcf`, the component filename is `calls.vcf.components.vcf`. Do not add counts across the VCF views. A deletion within an added copy is not an original-locus deletion; the reference `.indels.vcf` view handles that distinction for INS/DEL comparisons.

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

The base format is defined by the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf). AsmSV-specific fields and interpretations are declared in the output header.
