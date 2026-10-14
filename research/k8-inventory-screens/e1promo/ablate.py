"""Which constraint family alone pins each E1 entrance's backward ceiling to dimension 1 (expB).
Families: T = target ends (port hyperplanes), H = helper low ends (helper_final), K = low-end cones at end of B,
Z = A/B boundary pins (centre frames, target levels). Source and other helper ends are full and never bind."""
import sys, json
from collections import Counter
sys.path.insert(0, __import__('os').path.dirname(__file__))
import e1lib as L
w, sha = L.load(sys.argv[1]); W = L.Word(w); W.sha = sha; W.prepare(); W.compute_cones(sys.argv[2])
v = W.v; n = W.n
orig_end, orig_cone, orig_bd = list(W.end), dict(W.cone), dict(W.boundary)
def ceil_with(fams):
    W.end = [orig_end[q] if ((v <= q < 2*v and 'T' in fams) or (q >= 2*v and q in W.lowends and 'H' in fams)) else W.full for q in range(n)]
    W.cone = orig_cone if 'K' in fams else {}
    W.boundary = orig_bd if 'Z' in fams else {}
    # cone ceilings use end[q] of the low end q: keep the true low-end frames for K
    if 'K' in fams:
        save = W.end; W.end = list(save)
        for q in W.lowends: W.end[q] = orig_end[q]
        C, _ = W.ceilings(); W.end = save
        # recompute with the requested ends but keep coneceil from true low ends
        cc = W.coneceil
        W.cone = {}
        C, _ = W.ceilings(); W.coneceil = cc
        # apply cone ceilings manually: redo backward pass with cc
        Cm = list(W.end); gc = [None]*len(W.ev)
        for gi in range(len(W.ev)-1, -1, -1):
            sec, i, rs, U = W.ev[gi]
            if gi == W.nB - 1:
                for k in range(n): Cm[k] &= cc[k]
            if gi == W.nA - 1:
                for q, U0 in W.boundary.items(): Cm[q] &= U0
            c = L.FULLMASK
            for q in rs: c &= Cm[q]
            for q in rs: Cm[q] = c
        return Cm
    C, _ = W.ceilings(); return C
E1 = [q for q in W.entering if L.dim(W.start[q]) == 1]
res = {}
for fams in ['THKZ', 'T', 'H', 'K', 'Z', 'HKZ', 'TKZ', 'THZ', 'THK']:
    C = ceil_with(fams)
    res[fams] = Counter(L.dim(C[q]) for q in E1)
    if fams == 'THKZ': full = {q: L.dim(C[q]) for q in E1}
    if fams == 'T': tonly = {q: L.dim(C[q]) for q in E1}
for k, c in res.items(): print(k, dict(sorted(c.items())))
pinned = [q for q in E1 if full[q] == 1]
print('ceiling-1 helpers', len(pinned), 'of which pinned by target ends alone', sum(tonly[q] == 1 for q in pinned))
