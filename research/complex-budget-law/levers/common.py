"""expBK shared helpers: load a raw provider word, build the exact audit prefix (precondition + compensation + A + cut + B)
exactly as hpe_audit_and_transport.run does, end-of-B raw frames, and the exact state maps.  Raw coordinates throughout
(tau and the port mapping are linear relabelings, so cone / frame containment is invariant)."""
import gzip, json, importlib.util, sys
from pathlib import Path
sys.dont_write_bytecode = True
from fractions import Fraction as Q

AUDIT = str(Path(__file__).resolve().parents[3] / 'research/complex-571-cold14/code/hpe_audit_and_transport.py')
spec = importlib.util.spec_from_file_location('hpe_audit', AUDIT)
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
basis = A.basis


def load(path):
    return json.load(gzip.open(path))


def regs(e):
    return [e[1]] + ([e[2]] if e[0] == 'add' else [])


def compensation(word):
    v, R = word['report']['v'], word['report']['R']; n = 2*v + R
    main = word['A'] + word['cut'] + word['B']
    m = [{r: Q(1)} for r in range(n)]
    for e in main:
        A.apply(e, m)
    dirty = [{} for _ in range(n)]
    for r in range(v, 2*v):
        dirty[r] = {j: x for j, x in m[r].items() if j >= 2*v}
    for e in word['precondition']:
        A.apply(e, dirty)
    return [['add', r, j, -q.numerator, q.denominator, []] for r in range(v, 2*v) for j, q in sorted(dirty[r].items())]


def frames_trace(word):
    """start frames, per-register list of (part, index, frame) touches in A and B, end-of-B frames (raw basis tuples)."""
    v, R = word['report']['v'], word['report']['R']; n = 2*v + R
    cur = [() for _ in range(n)]
    for i, p in enumerate(word['ports']):
        cur[i] = basis([p])
    for r, f in word.get('start_extra', {}).items():
        cur[int(r)] = basis(f)
    start = list(cur)
    touches = [[] for _ in range(n)]
    for part in ('A', 'B'):
        for i, e in enumerate(word[part]):
            F = basis(e[-1])
            for r in regs(e):
                touches[r].append((part, i, F))
                cur[r] = F
    return start, touches, cur


def inside(a, b):
    return basis(list(a) + list(b)) == tuple(b)


def hyper(y):
    """basis of y^perp in F2^9"""
    return basis([x for x in range(1, 512) if bin(x & y).count('1') % 2 == 0])


def par(x):
    return bin(x).count('1') & 1


def in_yperp(F, y):
    return all(par(x & y) == 0 for x in F)
