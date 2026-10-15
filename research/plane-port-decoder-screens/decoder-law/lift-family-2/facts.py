import collections, itertools
from s1lib import *
print('|M|', set(len(Mm[b]) for b in N), 'O_hb = E-part + {b} + M_b^-:', all(set(CL[HB[b]]) & set(C) == {b} | set(Mm[b]) for b in N))
ic = collections.Counter((PP[b, c], len(set(Mm[b]) & set(Mm[c]))) for b in N for c in N if b != c)
print('(2 p_b.p_c, |M_b^- cap M_c^-|):', sorted(ic.items()))
ic = collections.Counter((PP[b, c], len(set(Mm[b]) & set(Mp[c]))) for b in N for c in N if b != c)
print('(2 p_b.p_c, |M_b^- cap M_c^+|):', sorted(ic.items()))
sz = collections.Counter((len(cp), s in cp) for cp in CPART.values()); print('C-part sizes', sorted(sz.items()))
# max set Q with pairwise 2pp = 1
best = 0
def ext(Q, cand):
    global best
    best = max(best, len(Q))
    for c in cand:
        if not Q or c > Q[-1]: ext(Q + [c], [d for d in cand if d > c and PP[c, d] == 1])
ext([], N); print('max pairwise-(+1/2) weight set', best)
