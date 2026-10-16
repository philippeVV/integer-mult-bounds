"""Port systems for expQS: P120, P112, P104 (#968 / biframe-decoder-ladder RESULTS.json), E7 (expQP construction).
B[a][b] = 1 iff required (or a == b); free iff B = 0 (odd inner product)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import json, os
from pathlib import Path
from e8 import *
REPO = Path(os.environ.get('IMB_REPO') or Path(__file__).resolve().parents[4])
RJ = REPO / 'research/biframe-decoder-ladder/RESULTS.json'

def line_of(v):
    v = tuple(v); nv = tuple(-x for x in v)
    for a in range(NL):
        if ROOT[a] == v or ROOT[a] == nv:
            return a
    raise KeyError(v)

def p104():
    d = json.load(open(RJ))['section4_explicit_P104']
    rm = {line_of(v) for v in d['F0']} | {line_of(v) for v in d['C']}
    assert len(rm) == 16
    return FULL & ~sum(1 << a for a in rm), rm

def e7():
    f = frames()[0]; alpha = f[7]
    return sum(1 << a for a in range(NL) if IP[a][alpha] == 0 and a != alpha)

def Bmat(L):
    idx = list(bits(L))
    return idx, [[(BROW[a] >> b) & 1 for b in idx] for a in idx]

def free_pairs(L):
    idx = list(bits(L)); out = []
    for i, a in enumerate(idx):
        for j in range(i + 1, len(idx)):
            if not (BROW[a] >> idx[j]) & 1:
                out.append((i, j))
    return out
