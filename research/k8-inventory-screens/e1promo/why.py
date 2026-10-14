"""For E1 helpers with ceiling dim 1: which hard constraint (boundary / cone / end) blocks each E' (expB)."""
import sys, json
from collections import Counter
sys.path.insert(0, __import__('os').path.dirname(__file__))
import e1lib as L
w, sha = L.load(sys.argv[1]); W = L.Word(w); W.sha = sha; W.prepare(); W.compute_cones(sys.argv[2]); C, gc = W.ceilings()
E1 = [q for q in W.entering if L.dim(W.start[q]) == 1]
stats = Counter(); per = {}
v = W.v
def kind(q):
    if q < v: return 'src'
    if q < 2*v: return 'tgt'
    if q in W.boundary: return 'centre'
    return 'hlp'
for q in E1:
    if L.dim(C[q]) != 1: continue
    E = W.start[q]; seen = set(); reasons = Counter()
    for x in range(1, 512):
        Ep = L.add(E, L.mask_of_basis((x,)))
        if Ep == E or Ep in seen: continue
        seen.add(Ep)
        raised, dprof, legal, bad, touched = W.complete({q: Ep})
        hard = [b for b in bad if b[0] != 'gate_ceiling']
        r = tuple(sorted({(b[0], kind(b[1])) for b in hard}))
        reasons[r] += 1
    per[q] = reasons
    stats[tuple(sorted(reasons))] += 1
    # first event of q
for k, c in stats.items(): print(c, k)
q = next(iter(per)); print('example', q, per[q])
sec, i, rs, U = W.ev[W.revs[q][0]]; print('first event', sec, i, rs, L.basis_of_mask(U), [kind(r) for r in rs])
