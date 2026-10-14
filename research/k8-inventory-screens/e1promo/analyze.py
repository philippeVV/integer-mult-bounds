"""Enumerate E1 entering helpers, their exact backward ceilings and every minimal E1->E2 promotion (expB).
usage: analyze.py WORD.json.gz OUT.json [cone-cache.json]"""
import sys, json, time
from collections import Counter
sys.path.insert(0, __import__('os').path.dirname(__file__))
import e1lib as L

t0 = time.time()
w, sha = L.load(sys.argv[1]); W = L.Word(w); W.sha = sha
paid = W.prepare(); cp = W.centre_paid(); tot = paid + cp
rep = {int(k): c for k, c in w['report']['local_profile'].items()}
print('word', sha[:12], 'profile recomputed', dict(sorted(tot.items())), 'report', rep, 'match', dict(tot) == rep, flush=True)
assert dict(tot) == rep
W.compute_cones(sys.argv[3] if len(sys.argv) > 3 else None)
print('cones', len(W.cone), 'time', round(time.time() - t0, 1), flush=True)
C, gc = W.ceilings()
# sanity: original frames under ceilings
for gi, e in enumerate(W.ev): assert e[3] & ~gc[gi] == 0, gi
for q in range(W.n): assert W.start[q] & ~C[q] == 0, q
E1 = [q for q in W.entering if L.dim(W.start[q]) == 1]
print('entering', len(W.entering), Counter(L.dim(W.start[q]) for q in W.entering), 'E1', len(E1), flush=True)
out = dict(word=sha, profile=dict(sorted(tot.items())), helpers={})
cdims = Counter()
for q in W.entering:
    cdims[(L.dim(W.start[q]), L.dim(C[q]))] += 1
print('(entrance dim, ceiling dim) counts', dict(sorted(cdims.items())), flush=True)
ends = {q: L.dim(W.end[q]) for q in E1}
for q in E1:
    E = W.start[q]; Cq = C[q]
    rec = dict(entrance=list(L.basis_of_mask(E)), end_dim=ends[q], width=ends[q] - 1, ceiling=list(L.basis_of_mask(Cq)),
               ceiling_dim=L.dim(Cq), first_events=[], candidates=[])
    for gi in W.revs[q][:4]:
        sec, i, rs, U = W.ev[gi]
        rec['first_events'].append(dict(sec=sec, idx=i, regs=rs, frame_dim=L.dim(U), ceil_dim=L.dim(gc[gi]), climb=L.dim(U) - L.dim(W.orig_before[(gi, q)])))
    # all 2-dim E' with E < E' <= ceiling
    seen = set()
    for x in L.members(Cq & ~E):
        Ep = L.add(E, L.mask_of_basis((x,)))
        if Ep in seen: continue
        seen.add(Ep)
        raised, dprof, legal, bad, touched = W.complete({q: Ep})
        dw = Counter({rec['width'] - 1: 1, rec['width']: -1})
        drc, parts = L.price(dprof, dw)
        rec['candidates'].append(dict(Ep=list(L.basis_of_mask(Ep)), legal=legal, bad=bad[:3], raised=len(raised),
                                      raised_secs=dict(Counter(W.ev[g][0] for g in raised)), dprofile={str(k): c for k, c in sorted(dprof.items())},
                                      dRC=drc, **parts, touched=len(touched)))
    rec['candidates'].sort(key=lambda c: c['dRC'])
    out['helpers'][str(q)] = rec
json.dump(out, open(sys.argv[2], 'w'), indent=1)
print('done', round(time.time() - t0, 1), 's', flush=True)
