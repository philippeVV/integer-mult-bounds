"""expXG core: gauge-invariant moment (Lovasz-type) bound Theta*_J for complex half-rate Tanner codes.

Seed = bipartite graph (points P, lines L, edge set E), local codes A_v <= C^{E(v)} of dim deg(v)/2.
W = U_P cap U_L,  U_P = (+)_p A_p,  U_L = (+)_l B_l.   s = dim W / |E|.
Gauge: positive diagonal D on C^E.  dim(D U_P cap D U_L) = dim W  (D invertible), moments are not invariant.
Theta_J(codes) = max{mu({1}) : mu prob. on [0,1] with moments m_1..m_J of the spectral measure of Pi_P Pi_L Pi_P}
              = 1 / (v^T H^{-1} v)   (Lukacs: J=2j uses H=Hankel(m_0..m_2j), J=2j+1 uses H=Hankel(m_1..m_2j+1), v=1)
Theta*_J = inf_D Theta_J(D codes)  >= s.
"""
import json, os, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    G = json.load(open(os.path.join(HERE, 'graphs_small.json')))[name]
    V = [list(x) for x in G['Pts'] + G['Lns']]
    return V, len(G['Pts']), G['E'], G['girth']


def orth(Y):
    Q, _ = np.linalg.qr(Y)
    return Q


def rand_codes(V, rng, cplx=True):
    out = []
    for es in V:
        n = len(es)
        Y = rng.standard_normal((n, n // 2)) + (1j * rng.standard_normal((n, n // 2)) if cplx else 0)
        out.append(orth(Y))
    return out


def fourier_codes(V, rng=None, shuffle=True):
    """flat-leverage codes: n/2 rows of the n-point DFT (random row subset and column labelling)."""
    out = []
    for es in V:
        n = len(es)
        F = np.exp(2j * np.pi * np.outer(np.arange(n), np.arange(n)) / n) / np.sqrt(n)
        rows = rng.permutation(n)[: n // 2] if (rng is not None and shuffle) else np.arange(n // 2)
        cols = rng.permutation(n) if (rng is not None and shuffle) else np.arange(n)
        out.append(F[np.ix_(cols, rows)])
    return out


def projectors(V, nP, E, Qs, logd=None):
    PP = np.zeros((E, E), complex)
    PL = np.zeros((E, E), complex)
    for i, (es, Q) in enumerate(zip(V, Qs)):
        Y = Q if logd is None else np.exp(logd[es])[:, None] * Q
        Y = orth(Y)
        blk = Y @ Y.conj().T
        if i < nP:
            PP[np.ix_(es, es)] = blk
        else:
            PL[np.ix_(es, es)] = blk
    return PP, PL


def moments(PP, PL, J):
    E = PP.shape[0]
    M = PP @ PL
    out = [1.0]
    X = np.eye(E, dtype=complex)
    for k in range(1, J + 1):
        X = X @ M
        out.append(float(np.real(np.trace(X))) / E)
    return np.array(out)


def markov(m, J):
    """max mass at 1 of a measure on [0,1] with moments m[0..J] (exact Lukacs/Christoffel form)."""
    if J % 2 == 0:
        j = J // 2
        H = np.array([[m[a + b] for b in range(j + 1)] for a in range(j + 1)])
    else:
        j = (J - 1) // 2
        H = np.array([[m[a + b + 1] for b in range(j + 1)] for a in range(j + 1)])
    v = np.ones(j + 1)
    try:
        x = np.linalg.solve(H, v)
    except np.linalg.LinAlgError:
        return 0.0
    val = v @ x
    if val <= 0:
        return 1.0
    return float(1.0 / val)


REFINE = int(os.environ.get('REFINE', '1'))


def refined(m):
    """kernel atom: Pi_P Pi_L Pi_P vanishes on U_P^perp (dim E/2), so nu = 1/2 delta_0 + 1/2 nu'; moments of nu' are 2 m_k."""
    return np.array([1.0] + [2 * v for v in m[1:]])


def markov_ref(m, J):
    return 0.5 * markov(refined(m), J) if REFINE else markov(m, J)


def theta(V, nP, E, Qs, J, logd=None):
    PP, PL = projectors(V, nP, E, Qs, logd)
    return markov_ref(moments(PP, PL, J), J)


def leverages(V, nP, E, Qs, logd=None):
    PP, PL = projectors(V, nP, E, Qs, logd)
    return np.real(np.diag(PP)), np.real(np.diag(PL))


def nullity(V, nP, E, Qs, tol=1e-8, logd=None):
    PP, PL = projectors(V, nP, E, Qs, logd)
    w = np.linalg.eigvalsh(2 * np.eye(E) - PP - PL)
    return int(np.sum(w < tol)), w


def flat_gauge(V, nP, E, Qs, side=0, iters=2000, tol=1e-12):
    """Barthe/Forster scaling: logd making every vertex code on `side` (0=P,1=L) have leverage 1/2 (Sinkhorn-type)."""
    logd = np.zeros(E)
    rng = range(nP) if side == 0 else range(nP, len(V))
    for it in range(iters):
        err = 0.0
        for i in rng:
            es = V[i]
            Y = orth(np.exp(logd[es])[:, None] * Qs[i])
            a = np.sum(abs(Y) ** 2, axis=1)
            target = len(es) / 2 / len(es)
            err = max(err, np.max(abs(a - target)))
            logd[es] += 0.5 * (np.log(target) - np.log(a))
        if err < tol:
            break
    return logd, err


def inner_min(V, nP, E, Qs, J, starts=(None,), maxiter=400):
    """inf over positive diagonal gauge of Theta_J (local search, L-BFGS with finite differences)."""
    from scipy.optimize import minimize
    best = (np.inf, None)
    for x0 in starts:
        x0 = np.zeros(E) if x0 is None else x0
        f = lambda x: theta(V, nP, E, Qs, J, x)
        r = minimize(f, x0, method='L-BFGS-B', options=dict(maxiter=maxiter, eps=1e-6))
        if r.fun < best[0]:
            best = (float(r.fun), r.x)
    return best


# ---------------- analytic gradients ----------------
def _coeffs(m, J):
    """Theta = 1/(v^T H^{-1} v); dTheta = Theta^2 * sum_k c_k dm_k."""
    if J % 2 == 0:
        j = J // 2; off = 0
    else:
        j = (J - 1) // 2; off = 1
    H = np.array([[m[a + b + off] for b in range(j + 1)] for a in range(j + 1)])
    y = np.linalg.solve(H, np.ones(j + 1))
    th = 1.0 / y.sum()
    c = np.zeros(J + 1)
    for a in range(j + 1):
        for b in range(j + 1):
            c[a + b + off] += y[a] * y[b]
    return th, th ** 2 * c


def theta_grad(V, nP, E, Ys, J):
    """Theta_J and the Hermitian gradient matrices G_P, G_L (dTheta = tr(G_P dPi_P) + tr(G_L dPi_L))."""
    A = np.zeros((E, E), complex); B = np.zeros((E, E), complex)
    for i, (es, Y) in enumerate(zip(V, Ys)):
        blk = Y @ Y.conj().T
        if i < nP: A[np.ix_(es, es)] = blk
        else: B[np.ix_(es, es)] = blk
    AB = A @ B
    pw = [np.eye(E, dtype=complex)]
    for k in range(J): pw.append(pw[-1] @ AB)
    m = np.array([1.0] + [float(np.real(np.trace(pw[k]))) / E for k in range(1, J + 1)])
    if REFINE:
        th, c = _coeffs(refined(m), J); th = th / 2   # d(th/2)/dm_k = (1/2) * c'_k * 2 = c'_k
    else:
        th, c = _coeffs(m, J)
    GP = np.zeros((E, E), complex); GL = np.zeros((E, E), complex)
    for k in range(1, J + 1):
        GP += c[k] * (k / E) * (B @ pw[k - 1])
        GL += c[k] * (k / E) * (pw[k - 1] @ A)
    GP = (GP + GP.conj().T) / 2; GL = (GL + GL.conj().T) / 2
    return th, A, B, GP, GL, m


def gauge_value_grad(V, nP, E, Qs, J, logd):
    Ys = [orth(np.exp(logd[es])[:, None] * Q) for es, Q in zip(V, Qs)]
    th, A, B, GP, GL, m = theta_grad(V, nP, E, Ys, J)
    g = np.real(np.diag(GP @ A + A @ GP - 2 * A @ GP @ A) + np.diag(GL @ B + B @ GL - 2 * B @ GL @ B))
    return th, g


def code_grad(V, nP, E, Ys, J):
    """gradient wrt orthonormal local bases Y_v (Euclidean, dTheta = Re tr(Z^* dY))."""
    th, A, B, GP, GL, m = theta_grad(V, nP, E, Ys, J)
    Zs = []
    for i, (es, Y) in enumerate(zip(V, Ys)):
        G = (GP if i < nP else GL)[np.ix_(es, es)]
        Pi = Y @ Y.conj().T
        Zs.append(2 * (np.eye(len(es)) - Pi) @ G @ Y)
    return th, Zs, m


def inner_min_grad(V, nP, E, Qs, J, starts, box=0.0, maxiter=500):
    from scipy.optimize import minimize
    best = None
    for st in starts:
        st = np.clip(st, -box, box) if box > 0 else st
        r = minimize(lambda x: gauge_value_grad(V, nP, E, Qs, J, x), st, jac=True, method='L-BFGS-B',
                     bounds=[(-box, box)] * E if box > 0 else None, options=dict(maxiter=maxiter))
        if best is None or r.fun < best.fun: best = r
    return float(best.fun), best.x
