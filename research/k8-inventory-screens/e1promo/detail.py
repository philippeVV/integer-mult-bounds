"""Print the chains of registers touched by one E1->E2 promotion (expB). usage: detail.py WORD cones helper [x]"""
import sys
sys.path.insert(0, __import__('os').path.dirname(__file__))
import e1lib as L
w, sha = L.load(sys.argv[1]); W = L.Word(w); W.sha = sha; W.prepare(); W.compute_cones(sys.argv[2]); C, gc = W.ceilings()
h = int(sys.argv[3]); E = W.start[h]
x = int(sys.argv[4]) if len(sys.argv) > 4 else min(L.members(C[h] & ~E))
Ep = L.add(E, L.mask_of_basis((x,)))
raised, dprof, legal, bad, touched = W.complete({h: Ep})
print('helper', h, 'E', L.basis_of_mask(E), "E'", L.basis_of_mask(Ep), 'ceil', L.basis_of_mask(C[h]), 'legal', legal, 'dprof', dict(dprof))
v = W.v
def kind(q): return 'src' if q < v else 'tgt' if q < 2*v else ('ctr' if q in W.boundary else 'hlp')
for q in sorted(touched):
    s = []
    prev = W.start[q]
    s.append('start d%d' % L.dim(prev))
    for gi in W.revs[q]:
        sec, i, rs, U = W.ev[gi]
        tag = '*' if gi in raised else ''
        d = L.dim(U) - L.dim(W.orig_before[(gi, q)])
        s.append('%s%d%s[%d,+%d%s]' % (sec, i, tag, L.dim(U), d, ('->%d' % L.dim(raised[gi])) if gi in raised else ''))
    s.append('end d%d' % L.dim(W.end[q]))
    print(q, kind(q), ' '.join(s))
for gi in sorted(raised):
    sec, i, rs, U = W.ev[gi]; print('raised', sec, i, w[sec][i] if sec != 'T' else w['early_cleanup'][i], '->', L.basis_of_mask(raised[gi]))
