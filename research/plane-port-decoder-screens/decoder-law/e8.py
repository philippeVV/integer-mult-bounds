"""Self-contained E8 line model. Lines = 120 root lines (doubled coords, ip = dot/4).
Pairing B_ab = (1 + ip_ab) mod 2 (diag: 1) -- equals label pairing q_a.q_b (checked vs pinned labels).
Term cost of support pair (R, A) = rank_F2 B[R, A]."""
import json, itertools
from fractions import Fraction as Fr
E8J = __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', 'lower', 'e8map.json')
_m = {int(k): v for k, v in json.load(open(E8J)).items()}
LAB = sorted(_m); ROOT = [tuple(_m[q]) for q in LAB]; NL = 120
def ip(x, y): return sum(a * b for a, b in zip(x, y)) // 4
def dotq(a, b): return bin(a & b).count('1') & 1
IP = [[ip(ROOT[a], ROOT[b]) for b in range(NL)] for a in range(NL)]
BROW = [sum(1 << b for b in range(NL) if (1 + IP[a][b]) % 2) for a in range(NL)]
assert all(((BROW[a] >> b) & 1) == dotq(LAB[a], LAB[b]) for a in range(NL) for b in range(NL))
def rank2(rows):
    basis = []
    for r in rows:
        for b in basis: r = min(r, r ^ b)
        if r: basis.append(r)
    return len(basis)
def bits(m):
    while m:
        low = m & -m; yield low.bit_length() - 1; m ^= low
def rho(R, A):
    """F2 pairing rank of B[R, A]; R, A bitmasks over the 120 lines."""
    return rank2([BROW[a] & A for a in bits(R)])
def supp(m, L):
    """lines a in L with r_a . m != 0 (m in doubled coords or any rational vector)."""
    return sum(1 << a for a in bits(L) if sum(x * y for x, y in zip(ROOT[a], m)) != 0)
FULL = (1 << NL) - 1
ORTH = [sum(1 << b for b in range(NL) if IP[a][b] == 0) for a in range(NL)]
def frames(L=FULL):
    out = []
    def ext(cl, cand):
        if len(cl) == 8: out.append(tuple(cl)); return
        for c in bits(cand):
            if c > cl[-1]: ext(cl + [c], cand & ORTH[c])
    for a in bits(L): ext([a], ORTH[a] & L)
    return out
def inv_dual(M):
    """M: list of 8 rational vectors (basis). Return dual vectors m'_k with m_j . m'_k = delta/... :
    sum_k m_k m'_k^T = I  (so N = R (1/2 * sum m_k m'_k^T) R^T uses terms (R m_k)(R m'_k)^T / 2 in true coords)."""
    n = 8; A = [[Fr(M[k][i]) for i in range(n)] + [Fr(int(i == k)) for i in range(n)] for k in range(n)]
    # rows = m_k ; solve for X with M X^T = I  -> X rows = dual
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    Minv = [row[n:] for row in A]   # Minv = M^{-1} where M has rows m_k: M Minv = I
    # dual vectors: columns of Minv  (m_j . col_k(Minv) = delta_jk)
    return [[Minv[i][k] for i in range(n)] for k in range(n)]
def gauge_terms(M, L=FULL):
    D = inv_dual(M)
    return [(supp(M[k], L), supp(D[k], L)) for k in range(8)]
def gauge_cost(M, L=FULL):
    return sum(rho(R, A) for R, A in gauge_terms(M, L))
def gauge_check(M, L=FULL):
    """exact check: N_ab = sum_k (r_a.m_k)(r_b.m'_k) * c equals G/2 (c chosen so), i.e. masked target holds."""
    D = inv_dual(M); Ls = list(bits(L))
    for a in Ls:
        for b in Ls:
            v = sum(Fr(sum(x * y for x, y in zip(ROOT[a], M[k]))) * sum(x * y for x, y in zip(ROOT[b], D[k])) for k in range(8)) / 8
            # doubled coords: r.m uses 2r, so (2r_a)^T (sum m m'^T) (2r_b) /8 = r_a.r_b/2 = ip/2 *... check
            if a == b: assert v == 1, v
            elif IP[a][b] == 0: assert v == 0
    return True
