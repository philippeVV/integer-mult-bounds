"""2-adic lower bound (exact, F2 linear algebra).
LEMMA: for any F2 vectors alpha_k subset supp d_k, beta_k subset supp t_k:
   sum_k rho_k >= rank_F2( B o (S (x) J2) ),  S = sum_k alpha_k beta_k^T  (mod 2),
   B = 224x224 F2 Gram of plane bases (block (a,b) = [x.y], x in E_a, y in E_b).
   Proof: B o (alpha beta^T (x) J) is a zero-padded submatrix of B[R_k,A_k] (rank rho_k); rank is subadditive.
THEOREM: if every term d_k t_k^T has 2-integral entries (no even denominators), take alpha=d_k mod 2, beta=t_k mod 2:
   S = N mod 2 = I on all required entries, and B_ab = 0 on free (orthogonal) pairs, so B o (S (x) J) = blockdiag(Gram E_a),
   rank 2*112 = 224.  Hence C >= 224 for every 2-integral decoder (naive N=I is optimal in that class).
This script: (1) randomized check of the lemma; (2) the theorem's block-diagonal rank; (3) naive decoder exact check."""
import random, sys
from base import *
def Bmat_rank(r, S):
    E = planes(r); rows = []
    for a in range(n):
        for x in E[a]:
            v = 0; bit = 0
            for b in range(n):
                for y in E[b]:
                    if S[a][b] and dot(x, y): v |= 1 << bit
                    bit += 1
            rows.append(v)
    return rank2(rows)
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
RC = [rchoices(q) for q in P]
worst = None
for trial in range(40):
    r = [rng.choice(RC[i]) for i in range(n)]
    K = rng.randint(1, 12); C = 0; S = [[0] * n for _ in range(n)]
    for k in range(K):
        R = rng.sample(range(n), rng.randint(1, 20)); A = rng.sample(range(n), rng.randint(1, 20))
        C += term_cost(r, R, A)
        al = [i for i in R if rng.random() < .7]; be = [j for j in A if rng.random() < .7]
        for i in al:
            for j in be: S[i][j] ^= 1
    lb = Bmat_rank(r, S)
    assert C >= lb, (C, lb)
    worst = min(worst or 99, C - lb)
print('PASS lemma randomized: 40 random term families, sum rho_k >= rank(B o S(x)J); min slack', worst)
r = [RC[i][0] for i in range(n)]
I = [[int(i == j) for j in range(n)] for i in range(n)]
print('blockdiag rank', Bmat_rank(r, I), '(= 224 for any planes: each Gram(E_a) nondegenerate)')
terms = [({i: 1}, {i: 1}) for i in range(n)]
print('naive decoder exact (bad, C):', verify(r, terms)[:2])
