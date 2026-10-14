"""Lovasz theta of Gamma_E (PSD-rank floor for symmetric PSD decoders) and alpha. usage: theta.py file.json ..."""
import sys, json, numpy as np, cvxpy as cp
from common import *
def theta(F, idx=None):
    idx = list(range(n)) if idx is None else idx; m = len(idx)
    B = cp.Variable((m, m), symmetric=True)
    cons = [B >> 0, cp.trace(B) == 1]
    cons += [B[i, j] == 0 for i in range(m) for j in range(i + 1, m) if F[idx[i], idx[j]]]
    p = cp.Problem(cp.Maximize(cp.sum(B)), cons); p.solve(solver='CLARABEL')
    return p.value
for f in sys.argv[1:]:
    d = json.load(open(f)); d = d[0] if isinstance(d, list) else d
    F = orthmask(d['r'])
    print(f, 'edges', F.sum() // 2, 'theta %.3f' % theta(F), flush=True)
