"""Diagonal embedding psi: (Sym_m(F2), c) -> (diag matrices over ANY field K, rank), m <= 3.
psi(beta) = diag(lam_1(beta), ..., lam_m(beta)) with F2-linear lam_i, lifted to {0,1} in K.
Requirements (README Sec.3):
  (L) 1-Lipschitz: rank(psi(x) - psi(y)) = #{i: lam_i(x) != lam_i(y)} <= c(x - y) for all x, y;
  (E) endpoint:    lam_i(I) = 1 for all i, so psi(beta + I) - psi(beta) = diag(+-1) is invertible in K.
Then GRS over K with the circuit's own K-realisation (T monomial) gives s >= Wm for every phase frame assignment.
lam_i(v v^T) = q_i(v) where q_i is the indicator of a 2-point coset (m=3) or a point (m=2).
Also checks: no diagonal-linear psi exists at m = 4 (exhaustive over the RM(2,4) support condition)."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
import beta as B


def lam_from_q(qfun, m):
    """linear functional on Sym_m(F2) (bitmask over pairs(m)) with lam(v v^T) = qfun(v)."""
    P = B.pairs(m); L = {}
    for a in range(m):
        e = 1 << a; L[(a, a)] = qfun(e)
    for a, b in P:
        if a < b: L[(a, b)] = (qfun((1 << a) | (1 << b)) - L[(a, a)] - L[(b, b)]) % 2
    lam = lambda mask: sum(L[P[t]] for t in range(len(P)) if mask >> t & 1) % 2
    for v in range(1, 1 << m):                       # verify lam(v v^T) = q(v)
        msk = sum(1 << t for t, (a, b) in enumerate(P) if (v >> a & 1) and (v >> b & 1))
        assert lam(msk) == qfun(v), (v,)
    return lam


SUPP = {2: [{0b01}, {0b10}],
        3: [{0b001, 0b111}, {0b010, 0b011}, {0b100, 0b110}]}
ok_all = True
for m, supps in SUPP.items():
    nb, rk, cost, I = B.tables(m)
    lams = [lam_from_q(lambda v, S=S: int(v in S), m) for S in supps]
    psi = lambda x: tuple(l(x) for l in lams)
    lip = all(sum(a != b for a, b in zip(psi(x), psi(y))) <= cost[x ^ y]
              for x in range(1 << nb) for y in range(1 << nb))
    endp = all(l(I) == 1 for l in lams)
    tight = sum(1 for x in range(1 << nb) for y in range(1 << nb)
                if sum(a != b for a, b in zip(psi(x), psi(y))) < rk[x ^ y])
    print(f"m={m}: Lipschitz(c) {lip}; endpoint lam_i(I)=1 {endp}; pairs where psi-rank < rank_F2: {tight}")
    ok_all &= lip and endp

# m = 4: need 4 functions q_i of degree <= 2, q_i(0) = 0, pairwise disjoint supports, each with an odd number of e_a.
m = 4
deg2 = []
for coeffs in itertools.product((0, 1), repeat=4 + 6):
    lin, quad = coeffs[:4], coeffs[4:]
    pr = [(a, b) for a in range(4) for b in range(a + 1, 4)]
    supp = frozenset(v for v in range(1, 16)
                     if (sum(lin[a] * (v >> a & 1) for a in range(4)) + sum(quad[t] * (v >> a & 1) * (v >> b & 1) for t, (a, b) in enumerate(pr))) % 2)
    if supp and sum(1 for a in range(4) if (1 << a) in supp) % 2 == 1: deg2.append(supp)
deg2 = list(set(deg2)); minsz = min(len(s) for s in deg2)
found = False
small = [s for s in deg2 if len(s) <= 15 - 3 * minsz]
for combo in itertools.combinations(sorted(small, key=len), 4):
    if sum(len(s) for s in combo) <= 15 and all(not (a & b) for a, b in itertools.combinations(combo, 2)):
        found = True; break
print(f"m=4: admissible degree-2 supports {len(deg2)}, min size {minsz}; 4 pairwise disjoint ones exist: {found}")
print("PASS_PSI" if ok_all and not found else "FAIL")
