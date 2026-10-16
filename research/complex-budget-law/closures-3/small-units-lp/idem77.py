"""Exact witness: the h = 6 column orbit 77 block (beta = 5, 16 sinks x 32 labels, rank_Q 16) carries an idempotent
X = D_r Bt D_a of rank 16 on its FULL support: solve Bt = Bt L Bt with L supported on types, build X, check X^2 = X."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import core6, flint
from rhoprime import system
D = core6.load(6); korb = sorted(set(D['corb'])); i = D['corb'].index(korb[77]); Wm, Qm = D['cols'][i]
W = core6.bits(Wm); Q = core6.bits(Qm); dot = core6.dot
S, A, b = system(W, Q)
M = flint.fmpq_mat(A); rhs = flint.fmpq_mat([[x] for x in b])
# least-norm-free exact solve: pick pivot columns via rref of [A | b]
R, rk = flint.fmpq_mat([r + [x] for r, x in zip(A, b)]).rref()
lam = [flint.fmpq(0)] * len(S); row = 0
piv = []
for c in range(len(S) + 1):
    if row < rk and R[row, c] != 0: piv.append(c); row += 1
assert len(S) not in piv, 'inconsistent'
for rI, c in enumerate(piv): lam[c] = R[rI, len(S)]
# ports: one per type s with r_s = 1, a_s = lam_s  (X = D_r Bt D_a on ports = types S)
X = flint.fmpq_mat([[dot(s[1], t[0]) * lam[j] for j, t in enumerate(S)] for s in S])
print('types used:', len(S), ' rank X =', X.rank(), ' X^2 == X:', X * X == X, ' trace =', sum((X[k, k] for k in range(len(S))), flint.fmpq(0)))
print('lambda values:', sorted({str(x) for x in lam}))
print('PASS_EXPNG_IDEM77' if len(S) == 256 and X.rank() == 16 and X * X == X else 'FAIL_EXPNG_IDEM77')
