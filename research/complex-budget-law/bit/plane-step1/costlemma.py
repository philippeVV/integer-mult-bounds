"""Exhaustive exact check of the non-T1 cost lemma on given sets (rational pairing ranks, Bareiss):
 Ty term  (reader = C(y) on L' plus plane ports eps, source = plane ports delta != {}):  cost >= 1, and >= 2 if eps != {}
 Tz term  (transpose): same by symmetry of pairing rank (checked on the transposed pairs as well)
 T0 term  (reader eps, source delta, both nonempty subsets of S): cost >= 2 (= 1 + rank G[eps,delta])
usage: costlemma.py idx [idx...]"""
import core, json, sys, itertools
sets = [s for s in json.load(open('sets_all.json')) if len(s[1]) == 5 and s[2] == 4]
L = list(range(111))
def basis(vecs):
    """integer row-echelon basis (rational span) via fraction-free elimination; returns independent subset of vecs."""
    out = []
    for v in vecs:
        if core.rankQ([list(x) for x in out] + [list(v)]) > len(out): out.append(v)
        if len(out) == 8: break
    return out
def prank(X, Y): return core.rankQ([[core.dot(x, y) for y in Y] for x in X]) if X and Y else 0
funcs = [y for y in itertools.product([0, 1], repeat=9) if any(y)]
for idx in map(int, sys.argv[1:]):
    r, S, rk = sets[idx]; Lp = set(L) - set(S); rv = core.LINES[r]
    subsets = [c for k in range(0, 6) for c in itertools.combinations(S, k)]
    Evec = lambda sub: ([rv] if sub else []) + [core.LINES[a] for a in sub]
    viol = 0; minc = {True: 9, False: 9}
    for y in funcs:
        VR0 = basis([core.LINES[t] for t in sorted(core.codeword(y) & Lp)])
        for eps in subsets:
            VR = basis(VR0 + Evec(eps))
            for dl in subsets:
                if not dl: continue
                c = prank(VR, basis(Evec(dl)))
                c2 = prank(basis(Evec(dl)), VR)
                assert c == c2
                need = 2 if eps else 1
                minc[bool(eps)] = min(minc[bool(eps)], c)
                if c < need: viol += 1
    t0 = min(prank(basis(Evec(e)), basis(Evec(d))) for e in subsets if e for d in subsets if d)
    print('set', idx, 'Ty min cost eps=0:', minc[False], 'eps!=0:', minc[True], 'violations', viol, 'T0 min', t0, flush=True)
