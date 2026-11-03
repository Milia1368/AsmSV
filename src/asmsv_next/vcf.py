"""Legacy VCF schema; candidate interpretation is retained."""
import os
import time
import logging
from pathlib import Path
import uuid
import pysam
from typing import List, Dict, Tuple, Any
from collections import defaultdict
from asmsv_next.structures import SVCandidate
from asmsv_next.config import SVConfig
from asmsv_next.vcf_annotations import compact_record,header_definitions,record_evidence
logger = logging.getLogger(__name__)
COMPLEMENT_MAP = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')

def reverse_complement(seq: str) -> str:
    """Computes the reverse complement of a nucleotide sequence."""
    if not seq:
        return ''
    return seq.translate(COMPLEMENT_MAP)[::-1]

def get_target_ref_seq(fasta: pysam.FastaFile, chrom: str, start: int, end: int, strand: str) -> str:
    """
    Safely retrieves the reference sequence for a specified region.
    Returns the reverse complement if the strand is negative.
    """
    try:
        seq = fasta.fetch(chrom, max(0, start), end).upper()
        if not seq:
            return ''
        return reverse_complement(seq) if strand == '-' else seq
    except (ValueError, KeyError, IndexError):
        return ''

def write_vcf(sv_list: List[SVCandidate], output_path: str, header_contigs: Dict[str, int], fasta: pysam.FastaFile,
              view='main',info_level=None,record_sink=None) -> None:
    target=Path(output_path);tmp=target.with_name(target.name+'.'+uuid.uuid4().hex+'.tmp')
    try:
        _write_vcf_body(sv_list,str(tmp),header_contigs,fasta,view,info_level or SVConfig.VCF_INFO,record_sink)
        with pysam.VariantFile(str(tmp)) as stream:
            for _ in stream:pass
        os.replace(tmp,target)
    finally:
        if tmp.exists():tmp.unlink()

def _write_vcf_body(sv_list, output_path, header_contigs, fasta,view='main',info_level='minimal',record_sink=None):
    """
    Generates a standard VCFv4.2 file, optimizing memory and disk I/O operations.
    Sorts SVs in biological order (chr1-22, X, Y, M).
    """

    def get_chrom_rank(chrom: str) -> Tuple[int, Any]:
        """Provides sorting weights for chromosome contigs."""
        c = str(chrom).upper()
        if c.startswith('CHR'):
            c = c[3:]
        if c.isdigit():
            return (0, int(c))
        elif c == 'X':
            return (1, 23)
        elif c == 'Y':
            return (1, 24)
        elif c in ('M', 'MT'):
            return (1, 25)
        else:
            return (2, chrom)
    chrom_ranks = {chrom: get_chrom_rank(chrom) for chrom in header_contigs.keys()}
    missing_chroms = {sv.chrom for sv in sv_list if sv.chrom not in chrom_ranks}
    for c in missing_chroms:
        chrom_ranks[c] = get_chrom_rank(c)
    sv_list.sort(key=lambda x: (chrom_ranks[x.chrom], x.pos))
    try:
        with open(output_path, 'w', buffering=10 * 1024 * 1024) as f:
            f.write('##fileformat=VCFv4.2\n')
            f.write(f"##fileDate={time.strftime('%Y%m%d')}\n")
            f.write('##source=AsmSV\n')
            f.write('##asmsvAnnotationSchema=asmsv-sv-annotation-v3\n')
            f.write(f'##asmsvAnnotationProfile={info_level}\n##asmsvView={view}\n')
            f.write('##asmsvEvidenceScope=SELECTED_PATH\n##asmsvAlternativeAssessment=NOT_AVAILABLE\n')
            f.write(f'##asmsvDupSourceSearch={SVConfig.COPY_SOURCE_SEARCH}\n')
            f.write(f'##asmsvGenotyping={SVConfig.GENOTYPING}\n')
            if SVConfig.GENOTYPING != 'presence':f.write('##phasing=partial\n')
            if SVConfig.CALLER_PRESET is not None:
                f.write(f'##asmsvPresetSHA256={SVConfig.CALLER_PRESET.fingerprint}\n')
            sorted_contigs = sorted(header_contigs.items(), key=lambda item: chrom_ranks[item[0]])
            contig_lines = [f'##contig=<ID={name},length={length}>\n' for (name, length) in sorted_contigs]
            f.writelines(contig_lines)
            alt_defs = [('INS', 'Insertion'), ('DEL', 'Deletion'), ('DUP', 'Duplication'), ('INV', 'Inversion'), ('BND', 'Breakend'), ('CPX', 'Complex structural variant; declared GATK-SV-style extension')]
            alt_defs = [item for item in alt_defs if item[0] in ('INS', 'DEL', 'DUP', 'INV', 'BND', 'CPX')]
            alt_lines = [f'##ALT=<ID={t},Description="{desc}">\n' for (t, desc) in alt_defs]
            alt_lines.append('##ALT=<ID=DUP:TANDEM,Description="Direct tandem duplication">\n')
            alt_lines.append('##ALT=<ID=DELINS,Description="Sequence replacement; REF and ALT lengths are separate">\n')
            f.writelines(alt_lines)
            info_defs = [('SVTYPE', 1, 'String', 'Type of structural variant (DEL, INS, DUP, INV, BND, TRA, CSV)'), ('SVLEN', 1, 'Integer', 'Signed length difference between REF and ALT alleles; negative for deletions'), ('END', 1, 'Integer', 'End position of the variant on the reference (1-based, inclusive); equals POS for breakends'), ('MATEID', 1, 'String', 'ID of the partner breakend record in a BND pair (VCF 4.2)'), ('SUPPORT_H1', 1, 'Integer', 'Number of supporting assembly read groups from haplotype 1 after cross-haplotype clustering'), ('SUPPORT_H2', 1, 'Integer', 'Number of supporting assembly read groups from haplotype 2 after cross-haplotype clustering'), ('RNAMES', 1, 'String', 'Comma-separated names of assembly reads (query sequences) supporting this variant'), ('SOURCE', 1, 'String', 'Detection layer: backbone (collinear split-read gap), CIGAR (single-segment CIGAR INDEL), CIGAR_MERGED (merged adjacent CIGAR INDELs), CIGAR_EDGE (source CIGAR event recovered outside the scoring view), COMPLEX (branch-cluster within backbone gap), COMPLEX_EDGE (edge branch cluster outside backbone)'), ('CHR2', 1, 'String', 'Chromosome of the mate breakend'), ('POS2', 1, 'Integer', '1-based reference position of the mate breakend'), ('STRANDS', 1, 'String', 'Strand orientations at POS and POS2 (first and second character, each + or -), used to construct the BND ALT allele per VCF 4.2'), ('SUBTYPE', 1, 'String', 'Mechanistic subclass. Linear: DUP:TANDEM, MULTI_DUP:TANDEM, INS_DEL. Inversion-related: INV, INV-INS, INV-DEL, INV-INS-DEL, InvDUP, InvDUP:TANDEM, InvDUP:PROXIMAL, InvDUP:DISP:INTRA. Dispersed duplication: DUP:PROXIMAL, DUP:DISP:INTRA, DUP:DISP:INTER. Other: LOCAL_RETAIN, COMPLEX'), ('DEL_LEN', 1, 'Integer', 'Length in base pairs of reference sequence deleted or replaced at this locus'), ('INS_LEN', 1, 'Integer', 'Length in base pairs of novel sequence inserted at this locus (query-derived)'), ('SEQ', 1, 'String', 'Inserted sequence; for INS records placed in the ALT column when below length threshold, for CSV(INS_DEL) records placed in INFO'), ('EVENT', '.', 'String', 'Shared event identifier linking related records from one complex rearrangement on an assembly read (BND pairs, span records, summary); comma-separated if merged across clusters'), ('INV_SPAN', 1, 'Integer', 'Length in base pairs of the inner inverted segment on the reference (inv_end minus inv_start)')]
            retained={'SVTYPE','SVLEN','END','MATEID','SUPPORT_H1','SUPPORT_H2','RNAMES','SOURCE',
                'CHR2','POS2','STRANDS','SUBTYPE','DEL_LEN','INS_LEN','SEQ','EVENT','INV_SPAN'}
            descriptions={
                'SVTYPE':'Variant class; CPX and DELINS are explicitly declared caller extensions',
                'SVLEN':'Signed SV length; INV uses inverted span; NET_LEN records sequence length change in structural summaries',
                'SOURCE':'Evidence layer; not a source chromosome or inferred mutational mechanism',
                'SUBTYPE':'Named duplication subtype: DUP:TANDEM dDUP piDUP_FR or piDUP_RF',
                'INS_LEN':'Inserted edit length; not necessarily novel sequence and not the complete duplicated template',
                'DEL_LEN':'Deleted or replaced edit length; added-duplicate edits do not imply source-locus loss',
                'SEQ':'Inserted sequence when represented in INFO; short INS and DELINS use explicit ALT sequence',
                'EVENT':'Canonical structural event identifier; never a claim about mutational chronology'}
            info_defs=[(k,n,t,descriptions.get(k,d)) for k,n,t,d in info_defs if k in retained]
            info_defs=[(k,'.' if k in ('MATEID','RNAMES') else n,t,d) for k,n,t,d in info_defs]
            info_defs.extend([
                ('CPX_TYPE',1,'String','GATK-SV named template when exact; COMPLEX is the declared unclassified fallback'),
                ('CPX_COMPONENTS','.','String','Ordered structural components using DEL INS INV DUP; repeated types are retained'),
                ('CPX_INTERVALS','.','String','Reference component intervals: TYPE_chr:start-end; 1-based inclusive'),
                ('INS_ORIGIN',1,'String','DUPLICATION TEMPLATED or UNKNOWN; origin annotation does not change an original CIGAR INS'),
                ('TEMPLATE_COUNT',1,'Integer','Number of explicitly represented source segments'),
                ('SOURCE_STATES','.','String','Per-template retained-source states in canonical target order'),
                ('ROLE',1,'String','SIMPLE EVENT COMPONENT JUNCTION or PROJECTION; projections are not independent events'),
                ('ATOM_ID','.','String','Stable input evidence identifiers'),
                ('PARENT_EVENT','.','String','Confirmed canonical parent events'),
                ('CONTEXT_ID','.','String','Candidate contexts; do not imply a shared event'),
                ('REGION_ID','.','String','Unresolved-region annotation IDs in companion regions.jsonl; not variant or event IDs'),
                ('RESOLUTION',1,'String','RESOLVED PARTIAL or UNRESOLVED'),
                ('DECOMPOSITION',1,'String','SINGLE SEPARABLE COUPLED or UNDETERMINED'),
                ('EVIDENCE_SCOPE',1,'String','Scope of evidence evaluated'),
                ('ALT_ASSESSMENT',1,'String','Whether alternative paths were evaluated'),
                ('ALLELE_CONTEXT',1,'String','REFERENCE ADDED_DUPLICATE or DUPLICATION_JUNCTION'),
                ('FRAME',1,'String','Geometric structure framework'),
                ('TOPOLOGY',1,'String','DUP when a retained-source duplication topology is confirmed'),
                ('SOURCE_STATE',1,'String','RETAINED or UNKNOWN; lack of source callability is not excision'),
                ('SOURCE_CHR',1,'String','Reference chromosome of the source'),
                ('SOURCE_START',1,'Integer','0-based source interval start'),
                ('SOURCE_END',1,'Integer','0-based exclusive source interval end'),
                ('TARGET_CHR',1,'String','Chromosome containing the new target'),
                ('TARGET_POS',1,'Integer','0-based target interbase port'),
                ('DUP_ORIENT',1,'String','DIRECT or INVERTED relative to forward target-reference order'),
                ('SOURCE_ORIENT',1,'String','DIRECT or INVERTED source orientation for templated sequence'),
                ('DUP_MODE',1,'String','TANDEM for a direct repeat junction; ADJACENT for other source-boundary placements; otherwise DISPERSED'),
                ('DUP_NORMALIZATION',1,'String','NONE MICROHOMOLOGY RIGHT_ANCHOR_EXTENSION LEFT_ANCHOR_EXTENSION or INS_ANCHOR_SHIFT'),
                ('JUNCTION_INS_LEN',1,'Integer','Known query insertion between duplicated occurrences; remote origin may remain unknown'),
                ('DUP_TEMPLATE_LEN',1,'Integer','Reference template length before copy-local edits'),
                ('DUP_SEQUENCE_LEN',1,'Integer','Length of the observed affected copy traversal'),
                ('DUP_OCCURRENCE',1,'Integer','Affected occurrence in query order; does not imply mutation chronology'),
                ('NET_LEN',1,'Integer','Signed sequence length change'),
                ('SOURCE_SEARCH_SCOPE',1,'String','OFF ADJACENT LOCAL or SELECTED_PATH; does not establish whole-genome source uniqueness'),
                ('SOURCE_ISSUE',1,'String','Unresolved source annotation reason'),
                ('RECHECK_METHOD',1,'String','EXACT KMER EXACT_INVERSE EXACT_SEPARATOR or optional WFA'),
                ('JUNCTION_ID',1,'String','One physical adjacency; BND mates share this identifier')])
            info_defs.append(('BND_INSERTION',1,'String','Known junction sequence in local reference orientation; encoded in BND ALT'))
            info_defs.extend([('QUERY_START',1,'Integer','0-based assembly support start'),('QUERY_END',1,'Integer','0-based exclusive assembly support end')])
            info_defs = header_definitions(info_defs, info_level)
            info_lines = [f'##INFO=<ID={i_id},Number={i_num},Type={i_type},Description="{i_desc}">\n' for (i_id, i_num, i_type, i_desc) in info_defs]
            f.writelines(info_lines)
            filter_defs = [('PASS', 'All filters passed'), ('OVERSIZE', 'Variant size exceeds the user-defined --max-sv threshold')]
            filter_defs = [item for item in filter_defs if item[0] in ('PASS', 'OVERSIZE')]
            for (fid, fdesc) in filter_defs:
                f.write(f'##FILTER=<ID={fid},Description="{fdesc}">\n')
            f.write('##FILTER=<ID=PARTIAL,Description="An observed junction from an incompletely resolved structure">\n')
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype; phased alleles follow Hap1 then Hap2 within PS; a missing allele is not reference">\n')
            if SVConfig.GENOTYPING != 'presence':
                f.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Chromosome-local phase set supported by one physical Hap1/Hap2 contig pair">\n')
                f.write('##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype evidence filter: NO_REFERENCE_SUPPORT, AMBIGUOUS_COVERAGE, NONREFERENCE_COMPONENT, or NO_ALT_SUPPORT; assembly mode retains dosage estimates">\n')
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+SVConfig.SAMPLE_NAME+'\n')
            type_counter = defaultdict(int)
            type_counter_all = defaultdict(int)
            count = 0

            def record_generator():
                nonlocal count
                for sv in sv_list:
                    try:
                        if not sv.id:
                            type_counter[sv.type] += 1
                            current_type_idx = type_counter[sv.type]
                        else:
                            current_type_idx = 0
                        record_tuple = sv.to_vcf_record(current_type_idx, alt_allele=True, reference=fasta)
                        if record_sink is not None:record_sink(record_evidence(sv,record_tuple,view))
                        record_tuple = compact_record(record_tuple, view, info_level)
                        type_counter_all[sv.type] += 1
                        count += 1
                        yield ('\t'.join(map(str, record_tuple)) + '\n')
                    except Exception as e:
                        logger.error(f'Writing VCF record failed for SV at {sv.chrom}:{sv.pos} | Error: {e}')
                        raise
            f.writelines(record_generator())
        logger.info(f'Successfully wrote {count} SVs to output file: {output_path}')
        logger.info(f'SV Counts by Type: {dict(type_counter_all)}')
    except IOError as e:
        logger.error(f"Failed to open or write to output VCF file '{output_path}': {e}")
        raise
