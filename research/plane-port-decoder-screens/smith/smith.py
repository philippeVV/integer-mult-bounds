"""2-adic Smith lemma for plane-port decoders (exact).
r_j(X) := number of Smith invariants of the integer/Z_(2) matrix X with 2-adic valuation <= j
        = min inner dimension c of a factorization X = U V^T mod 2^(j+1)   (subadditive).
LEMMA. N = sum_k d_k t_k^T over Q a masked decoder (N_aa=1, N_ab=0 when E_a not ⊥ E_b), e_k = v2(d_k)+v2(t_k).
       Let Bt = 224x224 0/1 Gram of the plane bases (q_a, r_a) (zero on free blocks, I_2 on diagonal blocks). Then
           224 <= sum_k r_{max(0,-e_k)}( Bt[R_k, A_k] ),
       and r_0(Bt[R,A]) = rho(R,A),  r_d(Bt[R,A]) <= g_d(rho) := sum_{s=1}^{d+1} C(rho, s)  (<= 2^rho - 1).
Proof sketch: m = max(0,-min e_k); A = Bt o (2^m N (x) J) = 2^m I_224 exactly; r_m(A)=224 <= sum r_m(A_k);
A_k = 2^(m+e_k) D_dbar Bt D_tbar => r_m(A_k) <= r_{-e_k}(Bt[R_k,A_k]).  g_d: n mod 2 = sum_{S!=0} ±2^{|S|-1} prod_{l in S} a_l."""
from fractions import Fraction as Fr
from math import comb
def v2(x):
    x = Fr(x)
    if x == 0: return 10**9
    a, b, v = x.numerator, x.denominator, 0
    while a % 2 == 0: a //= 2; v += 1
    while b % 2 == 0: b //= 2; v -= 1
    return v
def smith_vals(M, K):
    """valuations (< K) of Smith invariants of integer matrix M (list of lists) over Z/2^K; returns sorted list."""
    mod = 1 << K
    A = [[x % mod for x in row] for row in M]
    nr = len(A); nc = len(A[0]) if nr else 0
    vals = []; rows = list(range(nr)); cols = list(range(nc))
    def val(x):
        if x == 0: return K
        return (x & -x).bit_length() - 1
    while rows and cols:
        best = None
        for i in rows:
            Ai = A[i]
            for j in cols:
                x = Ai[j]
                if x:
                    v = val(x)
                    if best is None or v < best[0]:
                        best = (v, i, j)
                        if v == 0: break
            if best and best[0] == 0: break
        if best is None: break
        v, i, j = best; vals.append(v)
        p = A[i][j] >> v; pinv = pow(p, -1, mod)   # A[i][j] = 2^v * unit
        for i2 in rows:
            if i2 == i or A[i2][j] == 0: continue
            f = ((A[i2][j] >> v) * pinv) % mod
            Ai2, Ai = A[i2], A[i]
            for j2 in cols: 
                if Ai[j2]: Ai2[j2] = (Ai2[j2] - f * Ai[j2]) % mod
        rows.remove(i); cols.remove(j)
    return sorted(vals)
def r_j(M, j):
    return sum(1 for v in smith_vals(M, j + 1) if v <= j)
def g(d, rho): return sum(comb(rho, s) for s in range(1, d + 2))
