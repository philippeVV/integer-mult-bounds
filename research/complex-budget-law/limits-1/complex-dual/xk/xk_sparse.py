"""Sparse (block) version of expXG's Theta_J and its gradients, for large girth-8 graphs (W5: E = 936).
Same conventions as ../xg/xg_core.py (refined kernel-atom bound, REFINE=1)."""
import sys, os, json, numpy as np, scipy.sparse as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'xg'))
from xg_core import orth, _coeffs, refined, markov_ref, rand_codes, fourier_codes

HERE = os.path.dirname(os.path.abspath(__file__))


def load_any(name):
    for f in [os.path.join(HERE, '..', 'xg', 'graphs_small.json'), os.path.join(HERE, 'graphs_xk.json')]:
        if os.path.exists(f):
            d = json.load(open(f))
            if name in d:
                G = d[name]; V = [list(x) for x in G['Pts'] + G['Lns']]
                return V, len(G['Pts']), G['E'], G['girth']
    raise KeyError(name)


class Struct:
    def __init__(self, V, nP, E):
        self.V, self.nP, self.E = V, nP, E
        rP, cP, rL, cL = [], [], [], []
        for i, es in enumerate(V):
            for a in es:
                for b in es:
                    (rP if i < nP else rL).append(a); (cP if i < nP else cL).append(b)
        self.rP, self.cP, self.rL, self.cL = map(np.array, (rP, cP, rL, cL))

    def mats(self, Ys):
        dP, dL = [], []
        for i, (es, Y) in enumerate(zip(self.V, Ys)):
            (dP if i < self.nP else dL).append((Y @ Y.conj().T).ravel())
        A = sp.csr_matrix((np.concatenate(dP), (self.rP, self.cP)), shape=(self.E, self.E))
        B = sp.csr_matrix((np.concatenate(dL), (self.rL, self.cL)), shape=(self.E, self.E))
        return A, B


def tr(Xa, Xb):
    """tr(Xa Xb) for sparse matrices."""
    return complex(Xa.multiply(Xb.T).sum())


def masked(U, W, r, c):
    """entries (U @ W)[r_t, c_t] without forming U @ W."""
    return np.asarray(U[r].multiply(W.T.tocsr()[c]).sum(axis=1)).ravel()


def theta_grad(S, Ys, J):
    """returns th, m, block lists GPb (point blocks) and GLb (line blocks) of the Hermitian gradient (J <= 3)."""
    assert J <= 3
    E = S.E; A, B = S.mats(Ys)
    vA = np.asarray(A[S.rP, S.cP]).ravel(); vB = np.asarray(B[S.rL, S.cL]).ravel()
    X = (A @ B).tocsr(); Xh = (B @ A).tocsr()
    MP = [None, np.asarray(B[S.rP, S.cP]).ravel()]          # masked B (AB)^{k-1} on P pattern
    ML = [None, np.asarray(A[S.rL, S.cL]).ravel()]                                        # masked (AB)^{k-1} A on L pattern
    if J >= 2:
        BX = (B @ X).tocsr(); XA = (X @ A).tocsr()
        MP.append(np.asarray(BX[S.rP, S.cP]).ravel()); ML.append(np.asarray(XA[S.rL, S.cL]).ravel())
    if J >= 3:
        MP.append(masked(BX, X, S.rP, S.cP)); ML.append(masked(XA, Xh, S.rL, S.cL))
    # tr((AB)^k) = sum over P pattern of [B(AB)^{k-1}]_{ij} A_{ji};  A_{ji} = conj(A_{ij})
    m = np.array([1.0] + [float(np.real(np.sum(MP[k] * np.conj(vA)))) / E for k in range(1, J + 1)])
    th, c = _coeffs(refined(m), J); th = th / 2
    gP = sum(c[k] * (k / E) * MP[k] for k in range(1, J + 1)); gL = sum(c[k] * (k / E) * ML[k] for k in range(1, J + 1))
    Gb = []; oP = oL = 0
    for i, es in enumerate(S.V):
        n = len(es)
        if i < S.nP: g = gP[oP:oP + n * n].reshape(n, n); oP += n * n
        else: g = gL[oL:oL + n * n].reshape(n, n); oL += n * n
        Gb.append((g + g.conj().T) / 2)
    return th, m, Gb


def gauge_value_grad(S, Qs, J, logd):
    Ys = [orth(np.exp(logd[es])[:, None] * Q) for es, Q in zip(S.V, Qs)]
    th, m, Gb = theta_grad(S, Ys, J)
    g = np.zeros(S.E)
    for es, Y, G in zip(S.V, Ys, Gb):
        Pi = Y @ Y.conj().T
        g[es] += np.real(np.diag(G @ Pi + Pi @ G - 2 * Pi @ G @ Pi))
    return th, g


def code_grad(S, Ys, J):
    th, m, Gb = theta_grad(S, Ys, J)
    Zs = []
    for es, Y, G in zip(S.V, Ys, Gb):
        Pi = Y @ Y.conj().T
        Zs.append(2 * (np.eye(len(es)) - Pi) @ G @ Y)
    return th, Zs, m


def theta(S, Qs, J, logd=None):
    Ys = Qs if logd is None else [orth(np.exp(logd[es])[:, None] * Q) for es, Q in zip(S.V, Qs)]
    return theta_grad(S, Ys, J)[0]


def leverages(S, Qs, logd=None):
    a = np.zeros(S.E); b = np.zeros(S.E)
    for i, (es, Q) in enumerate(zip(S.V, Qs)):
        Y = orth(Q if logd is None else np.exp(logd[es])[:, None] * Q)
        (a if i < S.nP else b)[es] = np.sum(abs(Y) ** 2, axis=1)
    return a, b


def flat_gauge(S, Qs, side=0, iters=2000, tol=1e-12):
    logd = np.zeros(S.E); rng = range(S.nP) if side == 0 else range(S.nP, len(S.V))
    for it in range(iters):
        err = 0.0
        for i in rng:
            es = S.V[i]; Y = orth(np.exp(logd[es])[:, None] * Qs[i]); a = np.sum(abs(Y) ** 2, axis=1)
            err = max(err, np.max(abs(a - 0.5))); logd[es] += 0.5 * (np.log(0.5) - np.log(a))
        if err < tol: break
    return logd, err


def inner_min(S, Qs, J, starts, box=0.0, maxiter=500):
    from scipy.optimize import minimize
    best = None
    for st in starts:
        st = np.clip(st, -box, box) if box > 0 else st
        r = minimize(lambda x: gauge_value_grad(S, Qs, J, x), st, jac=True, method='L-BFGS-B',
                     bounds=[(-box, box)] * S.E if box > 0 else None, options=dict(maxiter=maxiter))
        if best is None or r.fun < best.fun: best = r
    return float(best.fun), best.x
