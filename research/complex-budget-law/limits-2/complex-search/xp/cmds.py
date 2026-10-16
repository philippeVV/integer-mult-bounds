"""Non-equivariant complex MDS Tanner-code search with a gauge-invariant MDS barrier and a canonical gauge.

Codes: A_v = col-span of Q_v (d x k, orthonormal), d-regular bipartite graph, k = d/2.  W = U_P cap U_L,
f_r = sum of the r smallest eigenvalues of K = 2I - Pi_P - Pi_L (= 0 iff dim W >= r).

Protocol (against Lemma D / equivariant degeneracy):
 * canonical gauge, re-imposed every iteration: P-codes Barthe-flat (all leverages k/d), then one scalar per point so that the
   mean L-leverage on each point star is 1/2 (doubly normalised).  Gauge drift cannot fake small eigenvalues of K.
 * MDS barrier: for every code v, its OWN Barthe-flattened gauge D_v (column scaling of that code alone, gauge invariant);
   barrier term -lam * sum_S log( N * |det (D_v Q_v)[S]|^2 / det(Q^H D_v^2 Q) ), annealed lam0 -> lam1 (floor kept > 0).
 * certificate: flattened Pluecker margin  mu_v = N * min_S |p_S|^2  (N = C(d,k)), spectral gap of K in canonical gauge.

usage: cmds.py GRAPHFILE NAME r seed [iters] [lam0] [lam1] [eta]   (env SAVE=1 saves the final codes)"""
import sys, json, itertools, os, time, numpy as np

gf, name, r, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
iters = int(sys.argv[5]) if len(sys.argv) > 5 else 3000
lam0 = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-2
lam1 = float(sys.argv[7]) if len(sys.argv) > 7 else 1e-6
eta = float(sys.argv[8]) if len(sys.argv) > 8 else 1.0
G = json.load(open(gf))[name]
EV = np.array(G['Pts'] + G['Lns']); nP = len(G['Pts']); V, d = EV.shape; k = d // 2; E = G['E']
SUBS = np.array(list(itertools.combinations(range(d), k))); N = len(SUBS)
rng = np.random.default_rng(seed)

def orth(Y):
    Q, R = np.linalg.qr(Y); return Q

def lev(Q):  # rows norms^2 of orthonormalised Q (batched V' x d x k)
    return np.sum(np.abs(orth(Q)) ** 2, axis=2)

def barthe(Q, its=200, tol=1e-10, x0=None):
    """per-code column log-scaling x (V' x d) making leverages k/d.  Returns x, flattened orthonormal Q, residual."""
    x = np.zeros(Q.shape[:2]) if x0 is None else x0.copy()
    for _ in range(its):
        M = orth(np.exp(x)[:, :, None] * Q); l = np.sum(np.abs(M) ** 2, axis=2)
        res = np.max(np.abs(l - k / d))
        if res < tol: break
        x -= 0.5 * np.log(l / (k / d)); x -= x.mean(axis=1, keepdims=True)
    return x, M, res

def pl_margin(Qf):  # Qf orthonormal (flattened); returns log of N*min|p_S|^2 per code
    dets = np.linalg.det(Qf[:, SUBS, :])
    return np.log(N * np.min(np.abs(dets) ** 2, axis=1))

def canon(Q):
    """impose the canonical gauge on all codes (in place on a copy). returns new Q, info."""
    Q = Q.copy()
    xP, MP, resP = barthe(Q[:nP]); Q[:nP] = MP
    # apply the same edge scaling to L codes
    y = np.zeros(E); y[EV[:nP].ravel()] = xP.ravel()
    Q[nP:] = orth(np.exp(y[EV[nP:]])[:, :, None] * Q[nP:])
    # doubly normalise: point scalars s_p so mean_{e in E(p)} b_e = 1/2
    s = np.zeros(nP)
    for it in range(60):
        b = np.zeros(E); b[EV[nP:].ravel()] = lev(Q[nP:]).ravel()
        c = b[EV[:nP]].mean(axis=1); dev = np.max(np.abs(c - 0.5))
        if dev < 1e-9: break
        step = -0.5 * np.log(c / 0.5); step -= step.mean()
        z = np.zeros(E); z[EV[:nP].ravel()] = np.repeat(step, d)
        Q[nP:] = orth(np.exp(z[EV[nP:]])[:, :, None] * Q[nP:])
        Q[:nP] = orth(np.exp(z[EV[:nP]])[:, :, None] * Q[:nP])  # uniform per star: P leverages unchanged
    return Q, dict(resP=float(resP), dn=float(dev))

def buildK(Q):
    K = 2 * np.eye(E, dtype=complex)
    P = Q @ np.conj(np.swapaxes(Q, 1, 2))
    for v in range(V):
        K[np.ix_(EV[v], EV[v])] -= P[v]
    return K

def barrier_grad(Q, x):
    """gradient (ascent direction) of sum_S log(|det Q_S|^2) - N log det(Q^H D^2 Q), D = exp(x) (own flattening)."""
    QS = Q[:, SUBS, :]                                   # V x N x k x k
    inv = np.linalg.inv(QS)                               # V x N x k x k
    g = np.zeros_like(Q)
    gS = np.conj(np.swapaxes(inv, 2, 3))                  # (Q_S^{-1})^H : V x N x k x k
    for t in range(k):
        np.add.at(g, (slice(None), SUBS[:, t]), gS[:, :, t, :])
    D2 = np.exp(2 * x)[:, :, None]
    Mq = np.conj(np.swapaxes(Q, 1, 2)) @ (D2 * Q)
    g -= N * (D2 * Q) @ np.linalg.inv(Mq)
    return g

def run():
    t0 = time.time()
    Q = orth(rng.standard_normal((V, d, k)) + 1j * rng.standard_normal((V, d, k)))
    Q, info = canon(Q)
    xall = None; log = []
    for it in range(iters):
        FR = float(os.environ.get('FRAC', '0.7')); lam = lam0 * (lam1 / lam0) ** min(1.0, it / (FR * iters))
        K = buildK(Q); w, U = np.linalg.eigh(K); X = U[:, :r]
        if it % 10 == 0 or xall is None:
            xall, _, _ = barthe(Q, its=100, x0=xall)
        Xv = X[EV]                                        # V x d x r
        Z = Xv @ np.conj(np.swapaxes(Xv, 1, 2))
        Dg = Z @ Q + lam * barrier_grad(Q, xall)
        Dg -= Q @ (np.conj(np.swapaxes(Q, 1, 2)) @ Dg)
        Q = orth(Q + eta * Dg)
        CE = int(os.environ.get('CANON_EVERY', '1'))
        if CE and it % CE == 0: Q, info = canon(Q)
        if it % 250 == 0 or it == iters - 1:
            f = float(w[:r].sum()); log.append((it, f))
            _, Qf_, _ = barthe(Q, its=300); mg_ = float(pl_margin(Qf_).min())
            print(json.dumps(dict(it=it, f=f, lam=lam, mlog=mg_, w_r=float(w[r - 1]), w_r1=float(w[r]), info=info,
                                  sec=round(time.time() - t0))), file=sys.stderr, flush=True)
    Q, info = canon(Q)
    K = buildK(Q); w = np.linalg.eigvalsh(K)
    xf, Qf, res = barthe(Q, its=500)
    mlog = pl_margin(Qf)
    levP = lev(Q[:nP]); levL = lev(Q[nP:])
    out = dict(graph=name, E=E, r=r, seed=seed, iters=iters, lam1=lam1, f=float(w[:r].sum()), w_r=float(w[r - 1]),
               gap=float(w[r]), nullity_1e8=int((w < 1e-8).sum()), nullity_1e5=int((w < 1e-5).sum()),
               min_log_flat_margin=float(mlog.min()), barthe_res=float(res), levL=[float(levL.min()), float(levL.max())],
               r_over_E=round(r / E, 5), Mg_E=round(E / 6, 2), sec=round(time.time() - t0))
    if int(os.environ.get('SAVE', '0')):
        np.save(f'sol_{name}_r{r}_s{seed}.npy', Q)
    return out

res_ = run(); print(json.dumps(res_), flush=True)
if os.environ.get('PASS_F'):   # limits-2 replay check: converged (f < PASS_F) with a healthy flattened MDS margin (>= e^-3)
    good = res_['f'] < float(os.environ['PASS_F']) and res_['min_log_flat_margin'] >= -3
    print(('PASS' if good else 'FAIL') + f"_EXPXP_CMDS_{name}_r{r}")
