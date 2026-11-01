"""Run-level MD/cs adapters and upstream provenance; no native-score scaling."""
from dataclasses import dataclass
import re
import shlex
from typing import Optional, Tuple
from asmsv_next.presets import mode_from_command


class TagError(ValueError):
    pass


@dataclass(frozen=True, slots=True)
class SourcePreset:
    aligner: str = 'unknown'
    version: Optional[str] = None
    mode: Optional[str] = None
    overrides: Tuple[str, ...] = ()
    command_line: Optional[str] = None
    program_id: Optional[str] = None


def source_from_program(program):
    name = str(program.get('PN', program.get('ID', 'unknown')))
    command = program.get('CL')
    try:
        argv = shlex.split(command or '')
    except ValueError:
        argv = []
    mode, overrides = mode_from_command(command), []
    for i, token in enumerate(argv):
        if token in ('-x', '-ax', '-cx') and i + 1 < len(argv):
            mode = argv[i + 1]
        for prefix in ('-ax', '-cx', '-x'):
            if token.startswith(prefix) and len(token) > len(prefix):
                mode = token[len(prefix):]
                break
        if token in ('-A', '-B', '-O', '-E', '--score-N', '--end-bonus'):
            overrides.extend(argv[i:i + 2])
        elif re.match(r'^-[ABOE].+', token) or token.startswith(('--score-N=', '--end-bonus=')):
            overrides.append(token)
    return SourcePreset(name, program.get('VN'), mode, tuple(overrides), command, program.get('ID'))


def source_from_header(header, program_id=None):
    programs = {p['ID']: p for p in header.get('PG', []) if p.get('ID')}
    if program_id in programs:
        seen = set()
        current = programs[program_id]
        while current.get('ID') not in seen:
            seen.add(current.get('ID'))
            name = str(current.get('PN', current.get('ID', ''))).lower()
            if any(tool in name for tool in ('minimap', 'winnowmap', 'vacmap', 'lra')):
                return source_from_program(current)
            if current.get('PP') not in programs:
                break
            current = programs[current['PP']]
    candidates = [p for p in programs.values() if any(tool in str(p.get('PN', p['ID'])).lower()
                  for tool in ('minimap', 'winnowmap', 'vacmap', 'lra'))]
    if len(candidates) == 1:
        return source_from_program(candidates[0])
    # A merged BAM with several aligners has no unambiguous global source preset.
    return SourcePreset()


def md_tokens(value):
    value = str(value)
    if not re.fullmatch(r'[0-9]+(?:(?:[A-Za-z]|\^[A-Za-z]+)[0-9]+)*', value):
        raise TagError('Malformed MD')
    result = []
    for token in re.findall(r'[0-9]+|\^[A-Za-z]+|[A-Za-z]', value):
        if token.isdigit():
            if int(token): result.append(('=', int(token)))
        elif token.startswith('^'):
            result.append(('D', len(token)-1))
        else:
            result.append(('X', 1))
    return result


def cs_tokens(value):
    value = str(value)
    pattern = re.compile(r':[0-9]+|=[A-Za-z]+|\*[A-Za-z]{2}|[+-][A-Za-z]+|~[A-Za-z]{2}[0-9]+[A-Za-z]{2}')
    cursor, result = 0, []
    for match in pattern.finditer(value):
        if match.start() != cursor:
            raise TagError('Malformed cs')
        token = match.group(); cursor = match.end()
        if token[0] == ':': op, n = '=', int(token[1:])
        elif token[0] == '=': op, n = '=', len(token)-1
        elif token[0] == '*': op, n = 'X', 1
        elif token[0] == '+': op, n = 'I', len(token)-1
        elif token[0] == '-': op, n = 'D', len(token)-1
        else: op, n = 'N', int(token[3:-2])
        if n <= 0:
            raise TagError('Nonpositive cs operation')
        result.append((op, n))
    if cursor != len(value) or not result:
        raise TagError('Empty or malformed cs')
    return result


class TokenCursor:
    def __init__(self, tokens):
        self.tokens = tokens
        self.index = 0
        self.offset = 0

    def take(self, length, allowed):
        result = []
        while length:
            if self.index >= len(self.tokens):
                raise TagError('Difference tag is shorter than CIGAR')
            op, n = self.tokens[self.index]
            if op not in allowed:
                raise TagError('Difference tag disagrees with CIGAR')
            used = min(length, n-self.offset)
            result.append((op, used))
            self.offset += used; length -= used
            if self.offset == n:
                self.index += 1; self.offset = 0
        return result

    def finish(self):
        if self.index != len(self.tokens):
            raise TagError('Difference tag is longer than CIGAR')


def synchronize(cigar, value, kind):
    """Return per-CIGAR aligned sub-runs; MD does not encode I or N."""
    cursor = TokenCursor(md_tokens(value) if kind == 'MD' else cs_tokens(value))
    result = []
    for op, length in cigar:
        if op in (0, 7, 8):
            allowed = {'=', 'X'} if op == 0 else ({'='} if op == 7 else {'X'})
            result.append(cursor.take(length, allowed))
        else:
            if op == 2 or (kind == 'cs' and op in (1, 3)):
                cursor.take(length, { {1:'I', 2:'D', 3:'N'}[op] })
            result.append(None)
    cursor.finish()
    return result
