"""Common assembly scoring; upstream identity is provenance, schema 5."""
from dataclasses import asdict, dataclass, fields
import hashlib
import json
import math
from pathlib import Path
import shlex
from asmsv_next.constants import MAPQ_FULL_WEIGHT_MIN, MAPQ_LOW_WEIGHT, REFINE_MIN_ANCHOR, DP_WEAK_TRA_PREMIUM_RATIO


def _nonnegative(name, value):
    if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
        raise ValueError(f'{name} must be finite and nonnegative')


@dataclass(frozen=True)
class AQPreset:
    mismatch: float = 19.0
    gap_open: tuple = (39.0, 81.0)
    gap_extend: tuple = (3.0, 1.0)
    structural_gap_min: int = 50
    structural_scale: float = 10.0

    def __post_init__(self):
        object.__setattr__(self, 'gap_open', tuple(self.gap_open))
        object.__setattr__(self, 'gap_extend', tuple(self.gap_extend))
        for key in ('mismatch', 'structural_scale'):_nonnegative(key,getattr(self,key))
        if len(self.gap_open)!=len(self.gap_extend) or not self.gap_open:
            raise ValueError('Gap open/extend arrays must have equal nonzero length')
        for x in self.gap_open+self.gap_extend:_nonnegative('gap coefficient',x)
        if type(self.structural_gap_min) is not int or self.structural_gap_min<=0:
            raise ValueError('structural_gap_min must be a positive integer')

    @property
    def min_anchor(self):
        # Keep the existing refine procedure unchanged; this is not a tuning key.
        return REFINE_MIN_ANCHOR

    def ordinary_gap_cost(self,length):
        return min(o+e*length for o,e in zip(self.gap_open,self.gap_extend)) if length else 0.0

    def gap_cost(self,length):
        t=self.structural_gap_min
        if length<t:return self.ordinary_gap_cost(length)
        return self.ordinary_gap_cost(t)+self.structural_scale*math.log1p((length-t)/t)

    def mapq_factor(self,quality):
        if quality==255:return None
        return 1.0 if quality>=MAPQ_FULL_WEIGHT_MIN else MAPQ_LOW_WEIGHT


@dataclass(frozen=True)
class DPPreset:
    query_gap: float = 0.25
    excess_overlap_rate: float = 0.5
    cross_chromosome: float = 500.0
    inversion: float = 200.0

    def __post_init__(self):
        for f in fields(self):_nonnegative(f.name,getattr(self,f.name))


@dataclass(frozen=True)
class CallerPreset:
    id: str
    aligner: str
    mode: str
    aq: AQPreset
    dp: DPPreset
    score_unit: str = 'match_equivalent_bp'
    schema_version: int = 5

    def __post_init__(self):
        if self.score_unit!='match_equivalent_bp':raise ValueError('Unsupported score unit')
        if self.schema_version!=5:raise ValueError('Only preset schema 5 is supported by this version')

    def to_dict(self):return asdict(self)

    def derived_parameters(self):
        c=self.dp.cross_chromosome
        return dict(structural_threshold_cost=self.aq.gap_cost(self.aq.structural_gap_min),
            structural_unit=self.aq.structural_gap_min,inversion_cost=self.dp.inversion,
            chromosome_cost=c,maximum_weak_TRA_premium=c*DP_WEAK_TRA_PREMIUM_RATIO,
            information_budgets=[c,2*c])

    @property
    def fingerprint(self):
        return hashlib.sha256(json.dumps(self.to_dict(),sort_keys=True,separators=(',',':')).encode()).hexdigest()

    @property
    def scoring_fingerprint(self):
        """Hash numerical evaluation independently of its source declaration."""
        aq=asdict(self.aq)
        for key in ('mismatch','structural_scale'):aq[key]=float(aq[key])
        for key in ('gap_open','gap_extend'):aq[key]=tuple(float(v) for v in aq[key])
        data=dict(aq=aq,dp={key:float(value) for key,value in asdict(self.dp).items()},score_unit=self.score_unit,
            schema_version=self.schema_version)
        return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def canonical_aligner(name):
    name=name.lower()
    if name=='winnowmap2': name='winnowmap'
    if name not in ('minimap2','winnowmap','vacmap'): raise ValueError('Unsupported aligner: '+name)
    return name


def resolve_preset(aligner,mode=None,*,preset_file=None):
    aligner=canonical_aligner(aligner)
    mode=mode or ('asm' if aligner=='vacmap' else 'asm5')
    modes=('asm',) if aligner=='vacmap' else ('asm5','asm10','asm20')
    if mode not in modes:raise ValueError(f'Unsupported {aligner} assembly source mode: {mode}')
    # Full-genome cross-over evaluation selects the same caller scoring for
    # all supported sources. @PG modes/cost flags describe how BAMs were made;
    # they never supply numerical defaults to this evaluator.
    result=CallerPreset('common-assembly-hq-v5',aligner,mode,AQPreset(),DPPreset())
    if preset_file:
        data=json.loads(Path(preset_file).read_text())
        if type(data.get('schema_version',5)) is not int or data.get('schema_version',5)!=5:
            raise ValueError('Historical preset schema is archived; use schema 5 or replay with its frozen source')
        allowed={f.name for f in fields(CallerPreset)}
        if set(data)-allowed:raise ValueError('Unknown or removed CallerPreset keys: '+str(sorted(set(data)-allowed)))
        if data.get('aligner',aligner)!=aligner or data.get('mode',mode)!=mode:
            raise ValueError('Custom preset source does not match the declared aligner/mode')
        values={}
        for key,cls in (('aq',AQPreset),('dp',DPPreset)):
            custom=data.get(key,{})
            if not isinstance(custom,dict):raise ValueError(key+' must be an object')
            extra=set(custom)-{f.name for f in fields(cls)}
            if extra:raise ValueError('Unknown or removed '+key+' parameters: '+str(sorted(extra)))
            values[key]=cls(**{**asdict(getattr(result,key)),**custom})
        result=CallerPreset(data.get('id','custom-assembly-v5'),aligner,mode,
            values['aq'],values['dp'],data.get('score_unit',result.score_unit))
    return result


def mode_from_command(command):
    argv=shlex.split(command or '')
    mode=None
    for i,token in enumerate(argv):
        if token in ('-x','-ax','-cx','-mode','--mode') and i+1<len(argv):mode=argv[i+1]
        elif token.startswith('--mode='):mode=token.split('=',1)[1]
        else:
            for prefix in ('-ax','-cx','-x'):
                if token.startswith(prefix) and len(token)>len(prefix):mode=token[len(prefix):]
    return mode
