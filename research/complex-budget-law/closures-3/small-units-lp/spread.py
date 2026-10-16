"""Linear-algebra check that per-total idempotency defects need not add up: exact rational U (v x K), V with U V^T = I_v,
K = v + 1 (slack s = 1), split into m blocks of rank k_j.  Every X_j = U_j V_j^T is non-idempotent, yet sum rank X_j = v + 1.
So 'v <= sum_j min(rank X_j, rho'_j)' is not implied by I = sum X_j; per-total idempotency is worth O(1) slack only."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import flint, random
random.seed(3); v, m = 12, 4; K = v + 1; ks = [4, 3, 3, 3]
U = flint.fmpq_mat(v, K, [random.randint(-3, 3) for _ in range(v * K)])
# V^T = right inverse of U: V^T = U^T (U U^T)^{-1}  (exact); then U V^T = I
VT = U.transpose() * (U * U.transpose()).inv()
assert U * VT == flint.fmpq_mat(v, v, [int(i == j) for i in range(v) for j in range(v)])
off = 0; tot = 0; allbad = True
for j, k in enumerate(ks):
    Uj = flint.fmpq_mat(v, k, [U[i, off + c] for i in range(v) for c in range(k)])
    Vj = flint.fmpq_mat(k, v, [VT[off + c, i] for c in range(k) for i in range(v)])
    X = Uj * Vj; r = X.rank(); tot += r
    print('block %d: rank X_j = %d, rank(X_j^2 - X_j) = %d, tr X_j = %s' % (j, r, (X * X - X).rank(), sum((X[i, i] for i in range(v)), flint.fmpq(0))))
    allbad &= (X * X - X).rank() == 1
    off += k
print('v = %d, sum rank X_j = %d (slack 1); every block non-idempotent' % (v, tot))
print('PASS_EXPNG_SPREAD' if tot == v + 1 and allbad else 'FAIL_EXPNG_SPREAD')
