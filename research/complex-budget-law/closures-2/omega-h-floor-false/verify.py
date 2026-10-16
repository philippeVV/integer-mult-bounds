"""Exact certificates against the omega*h floor C >= omega*h.
Model (expNZ/expSG, #446 Thm C cost): ports Q with pairing B (F2, unit diagonal); exact decoder N = sum_k p_k q_k^T over Q
with N_aa != 0 and N_ab = 0 whenever B_ab = 1 (a != b); cost = sum_k rank_F2 B[supp p_k, supp q_k].
omega = clique number of the required graph {B_ab = 1}; h = rank_F2 B (the least possible window).
Every decoder here is G/2 restricted to Q (frame gauge), checked entry by entry in Fractions."""
import itertools, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from e8 import *
from punct import omega, rkQ
fr = frames(); F0 = fr[0]; M = [ROOT[a] for a in F0]; T8 = gauge_terms(M)
def exact_check(L):
    D = inv_dual(M); Ls = list(bits(L))
    pv = {a: [Fr(sum(x*y for x, y in zip(ROOT[a], M[k]))) for k in range(8)] for a in Ls}
    qv = {a: [sum(Fr(x)*y for x, y in zip(ROOT[a], D[k])) for k in range(8)] for a in Ls}
    for a in Ls:
        for b in Ls:
            v = sum(pv[a][k]*qv[b][k] for k in range(8)) / 8
            if a == b: assert v == 1
            elif (BROW[a] >> b) & 1: assert v == 0
    # supports used in the cost are exactly the supports of these vectors
    for k, (R, A) in enumerate(T8):
        assert all(((R >> a) & 1) == (pv[a][k] != 0) and ((A >> a) & 1) == (qv[a][k] != 0) for a in Ls)
    return True
def report(name, L):
    exact_check(L)
    c = sum(rho(R & L, A & L) for R, A in T8); V = bin(L).count('1'); w = omega(L); h = rkQ(L)
    print(f'{name:28s} V={V:3d} omega={w} h={h} omega*h={w*h:3d}  C<={c:3d}  ' + ('VIOLATES C>=omega*h' if c < w*h else 'ok') + ('' if V >= w*h else '  (V<omega*h: trivial)'), flush=True)
    VIOL.append((name, c < w * h, V >= w * h))
    return c, w, h
VIOL = []
FULLc = report('P120', FULL)
report('P119 = P120 - 1 line', FULL & ~(1 << F0[0]))
FM = sum(1 << a for a in F0)
report('P112 = P120 - frame (v.t.)', FULL & ~FM)
for fn in ['best_1.json', 'best_2.json']:
    try: d = json.load(open(os.path.join(HERE, fn)))
    except FileNotFoundError: raise SystemExit('missing ' + fn)
    # anneal decoders use their own gauge; re-check with that gauge
    Mx = [[Fr(x) for x in m] for m in d['M']]; L = sum(1 << a for a in d['L'])
    M_save, T_save = M, T8
    M = [[int(x) if x.denominator == 1 else x for x in m] for m in Mx]; T8 = gauge_terms(M)
    report('anneal ' + fn, L); M, T8 = M_save, T8 if False else T_save
# E7 and F7 = E7 + all-ones label
alpha = F0[7]; E7 = sum(1 << a for a in range(120) if IP[a][alpha] == 0)
report('E7 (63 lines)', E7)
c = report('E7 - 1 line', E7 & ~(1 << F0[0]))
# F_7: append a port j whose B-row is all ones on E7 (label 1111111, see e7f7.py)
E7l = list(bits(E7)); n = len(E7l)
rowsF = [(BROW[a] & E7) for a in E7l]
idx = {a: i for i, a in enumerate(E7l)}
rows = [sum(1 << idx[b] for b in bits(r)) | (1 << n) for r in rowsF] + [(1 << (n + 1)) - 1]
hF = rank2(rows)
adj = [rows[i] & ~(1 << i) for i in range(n + 1)]
best = [0]
def ext(sz, cand):
    if sz > best[0]: best[0] = sz
    while cand:
        if sz + bin(cand).count('1') <= best[0]: return
        c = cand.bit_length() - 1; cand &= ~(1 << c); ext(sz + 1, cand & adj[c])
ext(0, (1 << (n + 1)) - 1)
cE7 = sum(rho(R & E7, A & E7) for R, A in T8)
print(f'F_7 = E7 + j (64 ports)          V=64 omega={best[0]} h={hF} omega*h={best[0]*hF}  C<={cE7+1}  ' + ('VIOLATES' if cE7 + 1 < best[0]*hF else 'ok'))
print('F_7 (all 64 odd labels of F2^7, dot form) = E7 + {1111111} (e7f7.py): omega=8, h=7, omega*h=56, '
      'decoder N_E7 (+) [1] is exact (block diagonal; the added port is required-adjacent to all, N=0 there) with C <= 49+1 = 50 < 56: VIOLATES')
assert cE7 + 1 < best[0] * hF                       # F_7 violates
assert [v for _, v, _ in VIOL] == [False, True, True, True, True, False, True], VIOL
assert all(big for name, v, big in VIOL if v)        # every violation has V >= omega*h (non-trivial)
print('PASS_EXPQP_VERIFY')
