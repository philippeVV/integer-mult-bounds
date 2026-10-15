"""(1) Generalised coverage identity priced with the pinned ceiling.py; (2) reductio for multi-total
(matrix-coefficient) term pricing used by PR2126; (3) rank-1-term price of the paired identity fit."""
import json
from fractions import Fraction as Q
from common import *

P = load_ports()
E8 = e8map()
n = len(P)
base = cap(WORD)
w = Q(FAM, ST)
print('baseline a_cap', float(base), 'w per call', float(w))

# (1) uniform u unpaid: break-even (largest integer C with ceiling > base), Theorem-A 2-integral floor 112u
out = {}
for u in (1, 2, 3, 4):
    # first-order break-even: u*w - C = w - 64
    fo = 64 + (u - 1) * w
    # exact: scan C downward from first-order with copies replaced by C rank-1 units? keep word unchanged
    # (favourable bracket used by levers.py: local word unchanged, only centre & unpaid change)
    lo, hi = 0, int(fo) + 5
    while hi - lo > 1:
        mid = (lo + hi) // 2
        try:
            ok = cap(WORD, centre=mid, unpaid=u) > base
        except Exception:
            ok = False
        if ok:
            lo = mid
        else:
            hi = mid
    x2 = None
    for Cc in range(int(fo), -1, -1):
        try:
            if cap(WORD, centre=Cc, unpaid=u) >= 2 * base:
                x2 = Cc; break
        except Exception:
            pass
    floorA = 112 * u
    out[u] = dict(first_order=float(fo), break_even_int=lo, kappa_x2_max_C=x2, thmA_floor=floorA,
                  gain_needed=floorA - lo, ceiling_at_C64=float(cap(WORD, 64, u)) if u < 73 else None)
    print('u', u, out[u])

# (2) reductio: line ports, N = G/2 is a rank-8 form. Under "one query carrying any coefficient matrix,
# cost = support pairing rank", a single query over all 112 ports would cost pair_rank(all, all).
G = [[sum(x * y for x, y in zip(E8[a], E8[b])) // 4 for b in P] for a in P]  # e8map stores doubled roots
assert all(G[i][i] == 2 for i in range(n))
mask_ok = all((G[i][j] == 0) for i in range(n) for j in range(n) if i != j and dot(P[i], P[j]) == 1)
rall = pair_rank(P, P)
print('G/2 is a masked line decoder:', mask_ok, ' F2 pairing rank of span(all labels):', rall, ' span dim', len(span(P)))
loc = {r: c for r, c in WORD.items() if r != 8}
loc[rall] = loc.get(rall, 0) + 1
a_red = cap(loc, centre=rall)
print('line word with one multi-total query: ceiling', float(a_red), 'x', float(a_red / base))

# (3) paired identity fit priced with rank-1 terms (one retained total per copy tree, #446 Thm C)
cert = json.load(open(MAIN + '/research/paired-plane-port-matching-seed6c563a/CERTIFICATE.json'))
cost_matrix_term = cost_rank1 = 0
for i, j in cert['matching']:
    assert dot(P[i], P[j]) == 0
    E = [P[i], P[j]]
    cost_matrix_term += pair_rank(E, E)          # PR2126: one I_2 term on the pair
    cost_rank1 += 2 * pair_rank(E, E)            # e_i e_i^T + e_j e_j^T, each E vs E
print('paired identity: PR2126 matrix-term invoice', cost_matrix_term, ' rank-1-term invoice', cost_rank1)
for Cc in (cost_matrix_term, cost_rank1):
    try:
        print('  ceiling(unpaid=2, C=%d)' % Cc, float(cap(WORD, Cc, 2)), 'x', float(cap(WORD, Cc, 2) / base))
    except Exception as e:
        print('  ceiling(unpaid=2, C=%d) fails:' % Cc, e)
json.dump(dict(uniform=out, reductio=dict(pair_rank_all=rall, ceiling=float(a_red), ratio=float(a_red / base)),
               paired=dict(matrix_term=cost_matrix_term, rank1=cost_rank1)), open('a_pricing.json', 'w'), indent=1)
