"""Convert a cell-tree design (pickle) into the request list consumed by build.py (exact Fractions)."""
import json, sys, pickle
from fractions import Fraction as Q
import celltree as CT
n = CT.n; q = CT.q; Gi = CT.D['G']
Cq = [{s: Q(-Gi[t][s], 2) for s in range(n) if s != t and Gi[t][s]} for t in range(n)]
J, zmap, glist = pickle.load(open(sys.argv[1], 'rb'))
tr = CT.Tree(glist, zmap); Jv, users, kids, per, _ = CT.evaluate(tr, detail=True)
print('J', Jv, 'expected L', Jv + n, 'pools', len(users))
def mask(S):
    sp = {0}
    for s in S:
        if q[s] not in sp: sp |= {x ^ q[s] for x in sp}
    m = 0
    for v in sp: m |= 1 << v
    return m
def basis(forms):
    """independent subset (exact) of forms (dicts)"""
    piv = []; out = []
    for f in forms:
        v = dict(f)
        for p, b in piv:
            if v.get(p):
                c = v[p]
                for k, x in b.items():
                    y = v.get(k, 0) - c * x
                    if y: v[k] = y
                    else: v.pop(k, None)
        if v:
            p = min(v); c = v[p]; piv.append((p, {k: x / c for k, x in v.items()})); out.append(f)
    return out
req = []
def emit(child, parent, forms):
    if len(child) == 1:
        s = child[0]
        fs = basis(forms); assert len(fs) == 1
        req.append([str(mask(child)), parent, {str(s): str(fs[0][s])}])
        return
    for f in basis(forms):
        req.append([str(mask(child)), parent, {str(k): str(v) for k, v in f.items()}])
for t in range(n):
    a, b = kids[('H', t)]
    for c in (a, b): emit(c, ['H', t], [{s: Cq[t][s] for s in c}])
for S, us in users.items():
    a, b = kids[S]
    for c in (a, b): emit(c, str(mask(S)), [{s: Cq[t][s] for s in c} for t in sorted(us)])
json.dump(dict(L=Jv + n, req=req), open(sys.argv[2], 'w'))
print('requests', len(req))
