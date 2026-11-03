"""Variant records, VCF coordinate conversion and optional BED exclusion."""
import os
import bisect
import logging
from typing import List, Dict, Tuple, Set, Optional, Union
from collections import defaultdict
import pysam
from asmsv_next.config import SVConfig
logger = logging.getLogger(__name__)

class SVCandidate:
    """
    Data structure representing an identified Structural Variant Candidate.
    """
    __slots__ = ('raw_type', 'type', 'chrom', 'pos', 'sv_len', 'read_names', 'info', 'filter', 'sources', 'genotype', 'end', 'q_start', 'q_end', 'mateid', 'id','reference_gap','format_data')

    def __init__(self, sv_type: str, chrom: str, pos: int, sv_len: int, read_name: Union[str, Set[str], List[str]], info: Dict=None, end: int=None, filter_str: str=None):
        _base = {}
        if info:
            _base.update(info)
        self.info = _base
        self.raw_type = sv_type
        self.type = self._normalize_type(sv_type)
        if self.raw_type != self.type and 'SUBTYPE' not in self.info:
            self.info['SUBTYPE'] = self.raw_type
        self.chrom = chrom
        self.pos = int(pos)
        self.sv_len = abs(int(sv_len))
        self.read_names = set()
        if read_name:
            if isinstance(read_name, str):
                self.read_names.add(read_name)
            else:
                self.read_names.update(read_name)
        self.filter = filter_str if filter_str else 'PASS'
        self.sources: Set[str] = set()
        self.genotype = './.'
        self.format_data = {}
        self.q_start: int = 0
        self.q_end: int = 0
        self.reference_gap = None
        self.mateid: str = ''
        self.id: str = ''
        if end is not None:
            self.end = int(end)
        elif 'END' in self.info:
            self.end = int(self.info['END'])
        elif self.type in ('BND', 'TRA'):
            self.end = self.pos
        elif self.type == 'INS':
            self.end = self.pos + int(self.info.get('REPLACED_REF_LEN', 0))
        else:
            self.end = self.pos + self.sv_len

    @staticmethod
    def _normalize_type(raw_type: str) -> str:
        rt = raw_type.upper()
        if rt == 'DELINS':
            return 'DELINS'
        if rt == 'CPX':return 'CPX'
        if 'CSV' in rt:
            return 'CSV'
        if 'DUP' in rt:
            return 'DUP'
        if 'Inv' in rt:
            return 'INV'
        if 'INV' in rt:
            return 'INV'
        if 'TRA' in rt:
            return 'TRA'
        if 'INS' in rt:
            return 'INS'
        if 'DEL' in rt:
            return 'DEL'
        return 'BND'

    def to_vcf_record(self, sv_id_index: int, alt_allele: bool=True, reference=None) -> Tuple:
        """Transforms candidate data into a VCF-compliant tuple record."""
        sv_id = self.id if self.id else f'SV.{self.type}.{sv_id_index:06d}'
        anchor_base = self._safe_fetch_ref(reference, self.chrom, max(0, self.pos - 1), self.pos) if reference else 'N'
        ref_val = anchor_base
        alt_val = f'<{self.type}>'
        info_list = [f'SVTYPE={self.type}', f'END={self.end}']
        if self.type not in ('BND', 'TRA'):
            val = int(self.info['NET_LEN']) if self.type in ('DELINS','CSV','CPX') and 'NET_LEN' in self.info else -self.sv_len if self.type == 'DEL' else self.sv_len
            info_list.append(f'SVLEN={val}')
        if self.type == 'BND':
            alt_val = self._format_bnd_alt(anchor_base)
        if self.type in ('BND', 'TRA'):
            if 'CHR2' in self.info:
                info_list.append(f"CHR2={self.info['CHR2']}")
            if 'POS2' in self.info:
                info_list.append(f"POS2={self.info['POS2']}")
        elif alt_allele:
            use_full_seq = SVConfig.ENABLE_LARGE_SEQUENCE or self.sv_len < SVConfig.MAX_ALLELE_LEN
            if self.type == 'INS':
                ins_seq = self.info.get('SEQ', '')
                alt_val = anchor_base + ins_seq if ins_seq and use_full_seq else '<INS>'
            elif self.type == 'DEL':
                if reference and use_full_seq:
                    deleted_seq = self._safe_fetch_ref(reference, self.chrom, self.pos, self.end)
                    if len(deleted_seq) == self.sv_len:
                        ref_val = anchor_base + deleted_seq
                        alt_val = anchor_base
                    else:
                        alt_val = '<DEL>'
                else:
                    alt_val = '<DEL>'
            elif self.type == 'DELINS':
                inserted=self.info.get('SEQ','')
                deleted=self._safe_fetch_ref(reference,self.chrom,self.pos,self.end) if reference else ''
                if inserted and len(inserted)<SVConfig.MAX_ALLELE_LEN and len(deleted)==self.end-self.pos and len(deleted)<SVConfig.MAX_ALLELE_LEN:
                    ref_val=anchor_base+deleted;alt_val=anchor_base+inserted
                else:alt_val='<DELINS>'
            elif self.type in ('DUP', 'INV', 'TRA', 'CSV','CPX'):
                alt_val = '<DUP:TANDEM>' if self.type=='DUP' and self.info.get('DUP_MODE')=='TANDEM' else f'<{self.raw_type}>'
        if self.read_names:
            self.info['RNAMES'] = ','.join(sorted(self.read_names))
        if self.type == 'CSV' and self.info.get('SUBTYPE') == 'INS_DEL' and self.info.get('SEQ'):
            info_list.append(f"SEQ={self.info['SEQ']}")
        if self.type=='DELINS' and alt_val=='<DELINS>' and self.info.get('SEQ'):
            info_list.append(f"SEQ={self.info['SEQ']}")
        exclude = {'SEQ', 'CHR2', 'POS2', 'END', 'SVLEN', 'SVTYPE', 'MATEID','STRUCTURE_KEY','SOURCE_SIGNATURE'}
        if self.mateid:
            info_list.append(f'MATEID={self.mateid}')
        for (k, v) in self.info.items():
            if k not in exclude and v is not None and (v != ''):
                if isinstance(v, float):
                    info_list.append(f'{k}={v:.3f}')
                else:
                    info_list.append(f'{k}={v}')
        fields = ['GT'] + [key for key in ('PS', 'FT') if key in self.format_data]
        values = [self.genotype] + [str(self.format_data[key]) for key in fields[1:]]
        return (self.chrom, str(self.pos), sv_id, ref_val, alt_val, '.', self.filter, ';'.join(info_list), ':'.join(fields), ':'.join(values))

    def _format_bnd_alt(self, ref_base: str) -> str:
        """Generates VCF 4.2 compliant BND ALT strings."""
        chr2 = self.info.get('CHR2')
        pos2 = self.info.get('POS2')
        if not chr2 or not pos2:
            return '<TRA>'
        strands = self.info.get('STRANDS', '+-')
        inserted = self.info.get('BND_INSERTION', '')
        p = f'{chr2}:{pos2}'
        if strands == '++':
            return f'{ref_base}{inserted}[{p}['
        elif strands == '--':
            return f']{p}]{inserted}{ref_base}'
        elif strands == '+-':
            return f'{ref_base}{inserted}]{p}]'
        elif strands == '-+':
            return f'[{p}[{inserted}{ref_base}'
        return '<TRA>'

    def _safe_fetch_ref(self, fasta, chrom: str, start: int, end: int) -> str:
        """Safely fetches reference sequences handling absent contigs seamlessly."""
        if not fasta:
            return 'N'
        try:
            seq = fasta.fetch(chrom, start, end).upper()
            return seq or 'N'
        except (ValueError, KeyError):
            alt_chrom = chrom.replace('chr', '') if 'chr' in chrom else f'chr{chrom}'
            try:
                seq = fasta.fetch(alt_chrom, start, end).upper()
                return seq or 'N'
            except (KeyError, ValueError):
                return 'N'

class BedMaskFilter:
    """BED file based region masking filter."""
    __slots__ = ('regions', 'region_starts', 'bed_path', 'is_active')

    def __init__(self, bed_path: Optional[str]):
        self.regions: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        self.region_starts: Dict[str, List[int]] = {}
        self.bed_path = bed_path
        self.is_active = False
        if bed_path and os.path.exists(bed_path):
            self._load_bed()
            self.is_active = True
            logger.info(f'Loaded mask regions from: {bed_path}')
        elif bed_path:
            logger.warning(f'Mask file not found: {bed_path}')

    def _load_bed(self):
        try:
            temp_regions = defaultdict(list)
            with open(self.bed_path, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.strip().split()
                    if len(parts) < 3:
                        continue
                    try:
                        temp_regions[parts[0]].append((int(parts[1]), int(parts[2])))
                    except ValueError:
                        continue
            for (chrom, coords) in temp_regions.items():
                if not coords:
                    continue
                coords.sort()
                merged = []
                (curr_start, curr_end) = coords[0]
                for (next_start, next_end) in coords[1:]:
                    if next_start <= curr_end:
                        curr_end = max(curr_end, next_end)
                    else:
                        merged.append((curr_start, curr_end))
                        (curr_start, curr_end) = (next_start, next_end)
                merged.append((curr_start, curr_end))
                self.regions[chrom] = merged
                self.region_starts[chrom] = [r[0] for r in merged]
        except Exception as e:
            logger.error(f'Failed to load BED mask: {e}')
            self.is_active = False

    def is_overlapped(self, chrom: str, start: int, end: int, threshold_ratio: float=0.5) -> bool:
        """Determines if a given interval significantly overlaps with masked regions."""
        if not self.is_active or chrom not in self.regions:
            return False
        sv_len = end - start
        (query_start, query_end) = (start, start + 1) if sv_len <= 0 else (start, end)
        query_len = max(1, sv_len)
        starts_idx = self.region_starts[chrom]
        idx = bisect.bisect_right(starts_idx, query_end)
        total_overlap = 0
        for i in range(idx - 1, -1, -1):
            (r_start, r_end) = self.regions[chrom][i]
            if r_end <= query_start:
                break
            (ov_s, ov_e) = (max(query_start, r_start), min(query_end, r_end))
            if ov_e > ov_s:
                total_overlap += ov_e - ov_s
        return total_overlap / query_len >= threshold_ratio

def filter_candidates_by_mask(candidates: List[SVCandidate], mask_filter: BedMaskFilter) -> List[SVCandidate]:
    """Filters a list of SVs dropping those falling into masked BED regions."""
    if not mask_filter.is_active:
        return candidates
    filtered_svs = []
    removed_count = 0
    mask_ratio_threshold = 0.5
    for sv in candidates:
        is_masked = mask_filter.is_overlapped(sv.chrom, sv.pos, sv.end, mask_ratio_threshold)
        if not is_masked and sv.type == 'BND':
            chr2 = sv.info.get('CHR2')
            pos2 = sv.info.get('POS2')
            if chr2 and pos2:
                try:
                    p2 = int(pos2)
                    if mask_filter.is_overlapped(chr2, p2 - 10, p2 + 10, 0.01):
                        is_masked = True
                except (ValueError, TypeError):
                    pass
        if not is_masked:
            filtered_svs.append(sv)
        else:
            removed_count += 1
    if removed_count > 0:
        logger.info(f'BedMaskFilter: Removed {removed_count} masked variants.')
    return filtered_svs
