"""expL8 small checks (numpy only, one thread, < 10 s).

A. Phase frames C_U = H^{(x)m} diag(i^{Q}) H^{(x)m} (Q: Z/4 quadratic forms on F_2^m) all commute with
   every translation X_v, C_v = aI + bX_v is a 2-term translation convolution, and the commutant of the
   whole phase-frame family is exactly the 2^m-dim translation algebra C[F_2^m]  (Lemma G).
B. Border lemma for the additive rank contract: rank(dA + t dB) for generic t >= rank(dB), so a frame
   family M(t) = t^-1 A + B never beats its finite part B (Lemma B1); kernel dimension is upper
   semicontinuous, so a limit code has dim >= generic dim (Lemma B2).
C. Nonlinear pointwise gates do not commute with a non-monomial frame (Lemma N1 demonstration).
"""
import itertools
import numpy as np

rng = np.random.default_rng(8)
ok = True


def walsh(m):
    H = np.array([[1, 1], [1, -1]], dtype=complex)
    out = np.array([[1]], dtype=complex)
    for _ in range(m):
        out = np.kron(out, H)
    return out


def xshift(m, v):
    n = 1 << m
    P = np.zeros((n, n), dtype=complex)
    for x in range(n):
        P[x ^ v, x] = 1
    return P


# ---- A
m = 3
n = 1 << m
Hm = walsh(m) / np.sqrt(n)
pts = list(range(n))
bits = lambda x: [(x >> j) & 1 for j in range(m)]
frames = []
# all Z/4-valued quadratic forms Q(x) = sum_i c_i x_i + 2 sum_{i<j} d_ij x_i x_j  (mod 4)
pairs = list(itertools.combinations(range(m), 2))
for c in itertools.product(range(4), repeat=m):
    for d in itertools.product(range(2), repeat=len(pairs)):
        Q = []
        for x in pts:
            b = bits(x)
            q = sum(ci * bi for ci, bi in zip(c, b)) + 2 * sum(dd * b[i] * b[j] for dd, (i, j) in zip(d, pairs))
            Q.append(q % 4)
        frames.append(Hm @ np.diag([1j ** q for q in Q]) @ Hm)
maxcomm = max(np.abs(F @ xshift(m, v) - xshift(m, v) @ F).max() for F in frames for v in range(n))
a, b = (1 + 1j) / 2, (1 - 1j) / 2
H2 = walsh(1) / np.sqrt(2)
C = H2 @ np.diag([1, 1j]) @ H2
c_ok = np.allclose(C, a * np.eye(2) + b * xshift(1, 1))
# commutant dimension of the frame family
rows = []
for F in frames[:: max(1, len(frames) // 40)]:
    rows.append(np.kron(np.eye(n), F) - np.kron(F.T, np.eye(n)))
K = np.vstack(rows)
sv = np.linalg.svd(K, compute_uv=False)
comm_dim = int((sv < 1e-9).sum()) + max(0, K.shape[1] - len(sv))
print(f"A: {len(frames)} phase frames, max |[C_U, X_v]| = {maxcomm:.1e}; C = aI+bX: {c_ok}; commutant dim = {comm_dim} (2^m = {n})")
ok &= maxcomm < 1e-12 and c_ok and comm_dim == n

# ---- B1
worst = 0
for trial in range(300):
    k = rng.integers(3, 7)
    rA, rB = rng.integers(0, k + 1), rng.integers(0, k + 1)
    dA = rng.integers(-2, 3, (k, rA)) @ rng.integers(-2, 3, (rA, k)) if rA else np.zeros((k, k))
    dB = rng.integers(-2, 3, (k, rB)) @ rng.integers(-2, 3, (rB, k)) if rB else np.zeros((k, k))
    t = rng.standard_normal()
    gen = np.linalg.matrix_rank(dA + t * dB)  # rank of t^-1 dA + dB up to scale
    worst = min(worst, gen - np.linalg.matrix_rank(dB))
print(f"B1: min over 300 random pencils of rank(dA+t dB) - rank(dB) = {worst} (>= 0 required)")
ok &= worst >= 0
# B1': with a t^-1 AND a t^+1 term the finite part CAN exceed the generic rank (border gain), e.g.
# dM(t) = t^-1 u w^T + (u z^T + v w^T) + t v z^T = t^-1 (u + t v)(w + t z)^T has rank 1, finite part rank 2.
u, v, w, z = (np.array(c, float) for c in ([1, 0], [0, 1], [1, 0], [0, 1]))
tt = 0.37
dM = np.outer(u, w) / tt + (np.outer(u, z) + np.outer(v, w)) + tt * np.outer(v, z)
r_gen, r_fin = np.linalg.matrix_rank(dM), np.linalg.matrix_rank(np.outer(u, z) + np.outer(v, w))
print(f"B1': Laurent edge t^-1 X + Y + t Z: generic rank {r_gen} < finite-part rank {r_fin} (border gain needs t^-1 terms)")
ok &= r_gen == 1 and r_fin == 2

# ---- B2: kernel dim of A + tB generic <= dim at t = 0
viol = 0
for trial in range(300):
    k = rng.integers(3, 7)
    r0 = rng.integers(0, k)
    A = rng.integers(-2, 3, (k, r0)) @ rng.integers(-2, 3, (r0, k)) if r0 else np.zeros((k, k))
    B = rng.integers(-2, 3, (k, k))
    t = rng.standard_normal() * 1e-3
    k_t = k - np.linalg.matrix_rank(A + t * B)
    k_0 = k - np.linalg.matrix_rank(A)
    viol += k_t > k_0
print(f"B2: generic kernel dim > limit kernel dim in {viol}/300 trials (0 required)")
ok &= viol == 0

# ---- C
x = rng.standard_normal(2) + 1j * rng.standard_normal(2)
g = lambda z: z + 0.1 * z * np.abs(z) ** 2  # a pointwise nonlinear 'gate'
err_nl = np.abs(g(C @ x) - C @ g(x)).max()
err_lin = np.abs((2 * (C @ x)) - C @ (2 * x)).max()
print(f"C: nonlinear gate vs C frame commutator = {err_nl:.3e}; linear gate = {err_lin:.1e}")
ok &= err_nl > 1e-3 and err_lin < 1e-12

print("PASS_EXPL8" if ok else "FAIL_EXPL8")
