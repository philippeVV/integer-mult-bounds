"""Shared geometry/structure for expSL (subspaces of F_2^9 as canonical basis tuples).
Inputs (read-only, repo-relative): ../../helper-resynthesis/expSK/{quads,alpha,instA}.json and
../../helper-bounds/wordmodel.py (basis)."""
import json, os, sys, itertools
from functools import lru_cache
from collections import defaultdict
_HERE = os.path.dirname(os.path.abspath(__file__))
_PKG = os.path.join(_HERE, '..', '..')                      # research/complex-budget-law
_DATA = os.path.join(_PKG, 'helper-resynthesis', 'expSK')
sys.path.insert(0, os.path.join(_PKG, 'helper-bounds')); from wordmodel import basis
from fractions import Fraction as F

d = json.load(open(os.path.join(_DATA, 'quads.json'))); V = d['V']; C = d['C']
rects = [(tuple(A), tuple(B)) for A, B in d['rects']]
al = json.load(open(os.path.join(_DATA, 'alpha.json'))); I = json.load(open(os.path.join(_DATA, 'instA.json'))); n = 112
NOR = V[112:]
dot = lambda a, b: (a & b).bit_count() & 1
cls = {x: i for i, c in enumerate(C) for x in c}


@lru_cache(maxsize=None)
def span(B):
    s = [0]
    for v in B: s = s + [x ^ v for x in s]
    return frozenset(s)


def sub(a, b):  # a subset of b
    sb = span(b); return all(v in sb for v in a)


def join(*Ss): return basis([v for S in Ss for v in S])


@lru_cache(maxsize=None)
def perp(B): return basis([x for x in range(1, 512) if all(dot(x, b) == 0 for b in B)])


def subspaces(B, k):
    """all k-dim subspaces of span(B) (canonical)."""
    out = set()
    vs = [v for v in span(B) if v]
    for S in itertools.combinations(vs, k):
        b = basis(list(S))
        if len(b) == k: out.add(b)
    return out


qk = lambda q: basis([V[x] for x in q])
LINE = lambda a: (V[a],)
TOP = lambda b: perp((V[b],))


def coef(b, s):
    if b < 112: return F(al[str(b)].get(str(s), '0'))
    return F(I['centres'][b - 112]['src'].get(str(s), 0))


byB = defaultdict(list)
for A, B in rects: byB[B].append(A)
crossB = [B for B, As in byB.items() if B[0] < 112 and len(As) == 2]
ownB = [B for B, As in byB.items() if B[0] < 112 and len(As) == 1]
cenB = [B for B, As in byB.items() if B[0] >= 112]
NB = {B: join(*[qk(A) for A in byB[B]]) for B in crossB}
tB = defaultdict(list)          # target -> cross B's containing it
for B in crossB:
    for t in B: tB[t].append(B)
ownQ = {}                         # target -> own A quad
for B in ownB:
    for t in B: ownQ[t] = byB[B][0]
