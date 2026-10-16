"""Non-equivariant complex MDS Tanner-code search with a gauge-invariant MDS barrier and a canonical gauge.

Codes: A_v = col-span of Q_v (d x k, orthonormal), d-regular bipartite graph, k = d/2.  W = U_P cap U_L,
f_r = sum of the r smallest eigenvalues of K = 2I - Pi_P - Pi_L (= 0 iff dim W >= r).

Protocol (against Lemma D / equivariant degeneracy):
 * canonical gauge, re-imposed every iteration: P-codes Barthe-flat (all leverages k/d), then one scalar per point so that the
   mean L-leverage on each point star is 1/2 (doubly normalised).  Gauge drift cannot fake small eigenvalues of K.
 * MDS barrier: for every code v, its OWN Barthe-flattened gauge D_v (column scaling of that code alone, gauge invariant);
   barrier term -lam * sum_S log( N * |det (D_v Q_v)[S]|^2 / det(Q^H D_v^2 Q) ), annealed lam0 -> lam1 (floor kept > 0).
 * certificate: flattened Pluecker margin  mu_v = N * min_S |p_S|^2  (N = C(d,k)), spectral gap of K in canonical gauge.

expXV variant: barrier ONLY on the scheduled in/out minors of a fixed Eulerian orientation (orient_NAME.json).
usage: cmds_io.py GRAPHFILE NAME r seed [iters] [lam0] [lam1] [eta]   (env SAVE=1 saves the final codes)"""
import sys, json, itertools, os, time, numpy as np

gf, name, r, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
iters = int(sys.argv[5]) if len(sys.argv) > 5 else 3000
lam0 = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-2
lam1 = float(sys.argv[7]) if len(sys.argv) > 7 else 1e-6
eta = float(sys.argv[8]) if len(sys.argv) > 8 else 1.0
G = json.load(open(gf))[name]
EV = np.array(G['Pts'] + G['Lns']); nP = len(G['Pts']); V, d = EV.shape; k = d // 2; E = G['E']
SUBS = np.array(list(itertools.combinations(range(d), k)))   # all k-subsets: used only to REPORT the full MDS margin
OR = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'orient_{name}.json')))
TAIL, HEAD = np.array(OR['tail']), np.array(OR['head'])
SUBSV = np.zeros((V, 2, k), dtype=np.int64)       # per code: local positions of in-ports, out-ports
for v in range(V):
    loc_in = [i for i, e in enumerate(EV[v]) if HEAD[e] == v]; loc_out = [i for i, e in enumerate(EV[v]) if TAIL[e] == v]
    assert len(loc_in) == k and len(loc_out) == k
    SUBSV[v, 0], SUBSV[v, 1] = loc_in, loc_out
N = 2; VI = np.arange(V)[:, None, None]
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

def pl_margin(Qf):  # in/out margin: log of 2*min(|det in|^2, |det out|^2) per code (flattened gauge)
    dets = np.linalg.det(Qf[VI, SUBSV, :])
    return np.log(N * np.min(np.abs(dets) ** 2, axis=1))

def full_margin(Qf):  # full MDS margin (report only): log of C(d,k)*min_S |p_S|^2
    dets = np.linalg.det(Qf[:, SUBS, :])
    return np.log(len(SUBS) * np.min(np.abs(dets) ** 2, axis=1)), (np.abs(dets) ** 2 * len(SUBS) < 1e-8).sum(axis=1)

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
    QS = Q[VI, SUBSV, :]                                  # V x 2 x k x k
    inv = np.linalg.inv(QS)
    g = np.zeros_like(Q)
    gS = np.conj(np.swapaxes(inv, 2, 3))
    for t in range(k):
        np.add.at(g, (VI[:, :, 0], SUBSV[:, :, t]), gS[:, :, t, :])
    D2 = np.exp(2 * x)[:, :, None]
    Mq = np.conj(np.swapaxes(Q, 1, 2)) @ (D2 * Q)
    g -= N * (D2 * Q) @ np.linalg.inv(Mq)
    return g

def run():
    t0 = time.time()
    Q = orth(rng.standard_normal((V, d, k)) + 1j * rng.standard_normal((V, d, k)))
    if os.environ.get('INIT'): Q = np.load(os.environ['INIT'])
    Q, info = canon(Q)
    xall = None; log = []
    for it in range(iters):
        FR = float(os.environ.get('FRAC', '0.7')); lam = lam0 * (lam1 / lam0) ** min(1.0, it / (FR * iters))
        K = buildK(Q); w, U = np.linalg.eigh(K); X = U[:, :r]
        if it % 10 == 0 or xall is None:
            xall, _, _ = barthe(Q, its=100, x0=xall)
            if int(os.environ.get('RAWBAR', '0')): xall = np.zeros_like(xall)   # barrier in the common canonical gauge
        Xv = X[EV]                                        # V x d x r
        Z = Xv @ np.conj(np.swapaxes(Xv, 1, 2))
        Dg = Z @ Q + lam * barrier_grad(Q, xall)
        Dg -= Q @ (np.conj(np.swapaxes(Q, 1, 2)) @ Dg)
        Q = orth(Q + eta * Dg)
        CE = int(os.environ.get('CANON_EVERY', '1'))
        if CE and it % CE == 0: Q, info = canon(Q)
        if it % int(os.environ.get('PEVERY', '250')) == 0 or it == iters - 1:
            f = float(w[:r].sum()); log.append((it, f))
            _, Qf_, _ = barthe(Q, its=300); mg_ = float(pl_margin(Qf_).min())
            dets0 = np.linalg.det(Q[VI, SUBSV, :]); raw_io = float(np.log(2 * np.min(np.abs(dets0) ** 2)))
            print(json.dumps(dict(it=it, f=f, lam=lam, mlog=mg_, raw_io=raw_io, xspread=float(np.ptp(xall, axis=1).max()), n12=int((w < 1e-12).sum()), n9=int((w < 1e-9).sum()), n6=int((w < 1e-6).sum()), w_r=float(w[r - 1]), w_r1=float(w[r]), info=info,
                                  sec=round(time.time() - t0))), file=sys.stderr, flush=True)
    Q, info = canon(Q)
    K = buildK(Q); w = np.linalg.eigvalsh(K)
    xf, Qf, res = barthe(Q, its=500)
    mlog = pl_margin(Qf); fm, nsing = full_margin(Qf)
    raw = pl_margin(Q)
    levP = lev(Q[:nP]); levL = lev(Q[nP:])
    out = dict(graph=name, E=E, r=r, seed=seed, iters=iters, lam1=lam1, f=float(w[:r].sum()), w_r=float(w[r - 1]),
               gap=float(w[r]), nullity_1e8=int((w < 1e-8).sum()), nullity_1e5=int((w < 1e-5).sum()),
               min_log_io_margin=float(mlog.min()), min_log_raw_io_margin=float(raw.min()), rawbar=int(os.environ.get('RAWBAR', '0')),
               n_1e12=int((w < 1e-12).sum()), n_1e9=int((w < 1e-9).sum()), min_log_full_margin=float(fm.min()),
               singular_ksubsets_total=int(nsing.sum()), nonMDS_codes=int((nsing > 0).sum()), barthe_res=float(res), levL=[float(levL.min()), float(levL.max())],
               r_over_E=round(r / E, 5), Mg_E=round(E / 6, 2), sec=round(time.time() - t0))
    if int(os.environ.get('SAVE', '0')):
        tag = os.environ.get('TAG', ''); np.save(f'sol_{name}_r{r}_s{seed}{tag}.npy', Q); np.save(f'solw_{name}_r{r}_s{seed}{tag}.npy', w)
    return out

print(json.dumps(run()), flush=True)
