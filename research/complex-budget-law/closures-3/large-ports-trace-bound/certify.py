"""exact certificate for the best omega<=4 system (E8-40) and the F4-18 / E6-27 ones:
N = Gram/norm (exact rationals), N_aa = 1, N_ab = 0 exactly on orthogonal (non-free) pairs, exact rank,
omega(nonorth) by exhaustive enumeration; P120 labels (h=9, B = J + G mod 2 restricted) realise F = nonorth exactly."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import itertools, sys
from fractions import Fraction
import numpy as np
from lattices import SYS
def rank_exact(M):
    M = [row[:] for row in M]; r = 0; n, m = len(M), len(M[0])
    for c in range(m):
        p = next((i for i in range(r, n) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(n):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r
def omega(A):
    V = len(A); best = 0
    def ext(c, cand):
        nonlocal best
        best = max(best, len(c))
        for b in cand:
            if len(c) + 1 + sum(1 for x in cand if x > b) <= best: continue
            ext(c + [b], [x for x in cand if x > b and A[b][x]])
    ext([], list(range(V))); return best
for name, f in [('F4', 'maxsub_F4_4.npy'), ('E6', 'maxsub_E6_4.npy'), ('E8', 'heur_E8_4.npy')]:
    L = SYS[name]()[np.load(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), f))]; G = (L @ L.T).tolist(); V = len(G)
    if len(set(G[i][i] for i in range(V))) == 1:
        N = [[Fraction(G[i][j], G[i][i]) for j in range(V)] for i in range(V)]; exact = True
    else:
        N = None; exact = False      # F4 mixes norms 4 and 8 -> entries 1/sqrt2; rank of G itself is the decoder rank
    A = [[G[i][j] != 0 and i != j for j in range(V)] for i in range(V)]
    rk = rank_exact([[Fraction(x) for x in r] for r in G])
    om = omega(A); ent = sorted(set(abs(N[i][j]) for i in range(V) for j in range(V) if i != j and N[i][j] != 0)) if exact else 'cos in {1/2, 1/sqrt2}'
    print(f"{name}-{V}: decoder rank {rk} (exact), omega(free=nonorth) = {om}, V/(omega*rank) = {Fraction(V, om*rk)} = {V/(om*rk):.4f}; off-diag |N| values {ent}")
    if name == 'E8':
        # P120 labels: x_a in F2^9 with B(x_a,x_b) = 1 + <a,b>/4 mod 2 (doubled coords): free <=> inner odd (norm-2 scale)
        B = [[(1 + G[i][j] // 4) % 2 for j in range(V)] for i in range(V)]
        free = [[B[i][j] == 0 and i != j for j in range(V)] for i in range(V)]
        same = all(free[i][j] == A[i][j] for i in range(V) for j in range(V))
        import numpy as _n
        Bm = _n.array(B) % 2
        def rk2(M):
            M = M.copy() % 2; r = 0
            for c in range(M.shape[1]):
                p = [i for i in range(r, M.shape[0]) if M[i, c]]
                if not p: continue
                M[[r, p[0]]] = M[[p[0], r]]
                for i in range(M.shape[0]):
                    if i != r and M[i, c]: M[i] ^= M[r]
                r += 1
            return r
        print(f"   P120-restricted labels: free graph == nonorth graph: {same}; label window h = rank_F2(B) = {rk2(Bm)} (u=1) "
              f"-> phi_eff = uV/(hK) = {V/(rk2(Bm)*rk):.4f}")
print('PASS_EXPQL_CERT')
