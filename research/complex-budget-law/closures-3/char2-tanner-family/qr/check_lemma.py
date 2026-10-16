"""Sanity checks of the algebra in README section 2 (small dense matrices, floats + exact triangular case):
(a) Z(N) = (I - K)A~ with K = A~ N A~^-1 - N  [Z_jk = A_jk (I + N_j - N_k)];
(b) single-session direction N = E_ii: nonzero spec K = +-sqrt(1 - A_ii (A^-1)_ii);
(c) degree-2 part of log det(I - K(D)) at m=1 is -d^T (I - B) d, B = A o A^{-T};
(d) triangular A: det Z(N) = det(A)^m for random non-commuting N."""
import numpy as np
rng = np.random.default_rng(0); n, m = 6, 3
A = rng.normal(size=(n, n)); At = np.kron(A, np.eye(m)); N = [rng.normal(size=(m, m)) for _ in range(n)]
Nb = np.zeros((n*m, n*m))
for k in range(n): Nb[k*m:(k+1)*m, k*m:(k+1)*m] = N[k]
Z = np.block([[A[j, k] * (np.eye(m) + N[j] - N[k]) for k in range(n)] for j in range(n)])
K = At @ Nb @ np.linalg.inv(At) - Nb
ok = np.allclose(Z, (np.eye(n*m) - K) @ At); print('(a)', ok)
Ai = np.linalg.inv(A); i = 2; E = np.zeros((n, n)); E[i, i] = 1
ev = np.linalg.eigvals(A @ E @ Ai - E); ev = ev[np.abs(ev) > 1e-9]
r = np.sqrt(complex(1 - A[i, i] * Ai[i, i])); print('(b)', np.sort_complex(ev), r)
ok &= len(ev) == 2 and np.allclose(np.sort(np.abs(ev)), [abs(r)] * 2)
d = rng.normal(size=n); B = A * Ai.T
f = lambda t: np.log(np.linalg.det(np.eye(n) - (A @ np.diag(t*d) @ Ai - np.diag(t*d))))
h = 1e-3; sec = (f(h) + f(-h) - 2*f(0)) / (2*h*h)
print('(c)', sec, -d @ (np.eye(n) - B) @ d, ' row sums of B:', B.sum(1).round(12))
ok &= abs(sec + d @ (np.eye(n) - B) @ d) < 1e-3 * max(1, abs(sec)) and np.allclose(B.sum(1), 1)
T = np.tril(rng.normal(size=(n, n))); Tt = np.kron(T, np.eye(m))
Z = np.block([[T[j, k] * (np.eye(m) + N[j] - N[k]) for k in range(n)] for j in range(n)])
print('(d)', np.linalg.det(Z), np.linalg.det(T)**m)
ok &= np.isclose(np.linalg.det(Z), np.linalg.det(T)**m)
print('PASS_EXPQR_LEMMA' if ok else 'FAIL_EXPQR_LEMMA')
