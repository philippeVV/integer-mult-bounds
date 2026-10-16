"""Degree-4 second-order model at two-sided flat codes.
Flat [4,2] code <-> w in simplex (Bloch vectors r_i in R^3, sum 0, r_i.r_j = 2w_k-1 on matching k).
Codes at points: flat(w_p) scaled by e^u ; at lines: flat(w_l) scaled by e^v.  u = x (gauge), v = x + delta.
Model: Theta_3 - 1/8 = (4/E) K(u,v) + O(3),  K = 1/4 z^T H z, z = ((I-X)u, (I-Y)v),
H = (1/32)[[4Y, 3(I+YX)],[3(I+XY), 4X]]."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
LIM1 = os.path.join(HERE, '..', '..', '..', 'limits-1', 'complex-dual')   # expXK / expXG copies in limits-1
sys.path.insert(0, os.path.join(LIM1, 'xk')); sys.path.insert(0, os.path.join(LIM1, 'xg'))
from xk_sparse import load_any, Struct, theta, orth
PAIRS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]

def flat_code(w, rng=None):
    c = 2 * np.asarray(w) - 1
    G = np.ones((4, 4))
    for k, pr in enumerate(PAIRS):
        for (i, j) in pr: G[i, j] = G[j, i] = c[k]
    ev, U = np.linalg.eigh(G); ev = np.clip(ev, 0, None)
    R = (U * np.sqrt(ev))[:, 1:]          # rows r_i in R^3 (drop the zero eigenvalue)
    if rng is not None:
        Qr = np.linalg.qr(rng.standard_normal((3, 3)))[0]; R = R @ Qr
    Y = np.zeros((4, 2), complex)
    for i, r in enumerate(R):
        r = r / np.linalg.norm(r); th = np.arccos(np.clip(r[2], -1, 1)); ph = np.arctan2(r[1], r[0])
        Y[i] = np.array([np.cos(th / 2), np.exp(1j * ph) * np.sin(th / 2)]) / np.sqrt(2)
    return Y

def XYmats(V, nP, E, W):
    X = np.zeros((E, E)); Y = np.zeros((E, E))
    for i, es in enumerate(V):
        T = X if i < nP else Y
        for k, pr in enumerate(PAIRS):
            for (a, b) in pr: T[es[a], es[b]] = T[es[b], es[a]] = W[i, k]
    return X, Y

def Kmat(X, Y):
    E = len(X); I = np.eye(E)
    H = np.block([[4 * Y, 3 * (I + Y @ X)], [3 * (I + X @ Y), 4 * X]]) / 32
    H = (H + H.T) / 2
    Z = np.block([[I - X, 0 * I], [0 * I, I - Y]])
    return Z.T @ H @ Z / 4          # quadratic form in (u, v)

def xdelta(K):
    E = len(K) // 2; I = np.eye(E); T = np.block([[I, 0 * I], [I, I]])
    Kp = T.T @ K @ T
    return Kp[:E, :E], Kp[:E, E:], Kp[E:, E:]

def schur(K, rc=1e-10):
    G, C, D = xdelta(K)
    ev = np.linalg.eigvalsh(G)
    Sm = D - C.T @ np.linalg.pinv(G, rcond=rc) @ C
    return ev, np.linalg.eigvalsh((Sm + Sm.T) / 2), Sm

if __name__ == '__main__':
    name = sys.argv[1] if len(sys.argv) > 1 else 'W3'
    V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(0)
    for trial in range(3):
        W = rng.dirichlet([2, 2, 2], size=len(V)) if trial else np.full((len(V), 3), 1 / 3)
        Q0 = [flat_code(W[i], rng) * np.exp(2j * np.pi * rng.random(4))[:, None] for i in range(len(V))]
        X, Y = XYmats(V, nP, E, W); K = Kmat(X, Y)
        u = rng.standard_normal(E); v = rng.standard_normal(E)
        for tau in [1e-2, 3e-3]:
            Qs = [Q0[i] * np.exp(tau * (u if i < nP else v)[es])[:, None] for i, es in enumerate(V)]
            th = theta(S, Qs, 3, np.zeros(E)); z = np.r_[u, v] * tau
            print(f'trial {trial} tau {tau}: true {th - 0.125:.4e}  model {4 / E * z @ K @ z:.4e}')
        ev, es_, _ = schur(K)
        print('  Gamma eig min (excl. const)', np.sort(ev)[:3], ' Schur max', es_[-3:])
