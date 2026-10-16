"""Riemannian search over Singer-equivariant complex [18,9] code pairs (A,B) on PG(2,17).
Loss = sum_chi sum_{i<r_chi} lambda_i(2I - P_A - Λ_chi P_B Λ_chi^H)  (zero  <=>  dim(A∩Λ_chi B) >= r_chi).
Gauge: common positive diagonal; A is Barthe-flattened (leverage 1/2) after every step, so B carries the
non-flat leverage.  MDS barrier  -mu * mean_S log|det A_S|^2/|...|  over all 48620 9-subsets (both codes).
usage: numopt.py PROFILE seed iters [mu0]   PROFILE: c<r> (constant r), m<r>:<m> (r on m chars), tri<W>"""
import sys, itertools, time, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from geom import singer
N, D = singer(17); n = 18; k = 9; E = N*n
D = np.array(D)
prof = sys.argv[1]; seed = int(sys.argv[2]); iters = int(sys.argv[3]); mu0 = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-2
rng = np.random.default_rng(seed)
chis = np.arange(N)
if prof.startswith('c'): r = np.full(N, int(prof[1:]))
elif prof.startswith('m'):
    rr, m = map(int, prof[1:].split(':')); r = np.zeros(N, int); r[rng.permutation(N)[:m]] = rr
elif prof.startswith('tri'):
    W = int(prof[3:]); dist = np.minimum(chis, N-chis); r = np.maximum(0, 9 - np.ceil(9*dist/W)).astype(int)
r = r.astype(int); print('profile', prof, 'target sum', r.sum(), flush=True)
Lam = np.exp(2j*np.pi*np.outer(chis, D)/N)            # (N,n)
SUBS = np.array(list(itertools.combinations(range(n), k)))
def orth(X):  # rows span -> orthonormal columns Q (n,k)
    Q, _ = np.linalg.qr(X.conj().T); return Q
def lev(Q): return (np.abs(Q)**2).sum(1)
def flatten(QA, QB, it=200):
    c = np.ones(n)
    for _ in range(it):
        l = lev(QA)
        if np.abs(l-0.5).max() < 1e-12: break
        s = np.sqrt(0.5/l); c *= s
        QA = orth((QA*s[:, None]).T); QB = orth((QB*s[:, None]).T)
    return QA, QB
def plucker_logs(Q):
    M = Q[SUBS].transpose(0, 1, 2)                     # (S,k,k) rows = coords
    sign, ld = np.linalg.slogdet(M)
    return ld
def barrier_grad(Q):
    M = Q[SUBS]                                          # (S,k,k)
    inv = np.linalg.inv(M)                               # d log|det M|^2 / d conj(M) = inv^H
    G = np.zeros_like(Q)
    np.add.at(G, SUBS, np.conj(inv.transpose(0, 2, 1)))
    return G / len(SUBS)
def loss_grad(QA, QB, want_grad=True):
    PA = QA @ QA.conj().T; PB = QB @ QB.conj().T
    LPB = Lam[:, :, None] * PB[None] * Lam.conj()[:, None, :]
    K = 2*np.eye(n)[None] - PA[None] - LPB
    w, V = np.linalg.eigh(K)
    mask = np.arange(n)[None, :] < r[:, None]
    f = (w*mask).sum()
    if not want_grad: return f, w
    wt = mask / np.sqrt(np.maximum(w, 0) + EPS**2) / 2    # gradient of sum sqrt(lambda+eps^2)
    Vm = V * np.sqrt(wt)[:, None, :]
    S = Vm @ Vm.conj().transpose(0, 2, 1)                # (N,n,n) sum v v^H over targeted
    GA = -2 * S.sum(0) @ QA
    SB = (Lam.conj()[:, :, None] * S * Lam[:, None, :]).sum(0)
    GB = -2 * SB @ QB
    return f, w, GA, GB
def hproj(Q, G): return G - Q @ (Q.conj().T @ G)
QA = orth(rng.normal(size=(k, n)) + 1j*rng.normal(size=(k, n)))
QB = orth(rng.normal(size=(k, n)) + 1j*rng.normal(size=(k, n)))
QA, QB = flatten(QA, QB)
EPS = 1e-2; eta = 0.01; mu = mu0; mA = 0; mB = 0; t0 = time.time(); fprev = None
for it in range(iters):
    f, w, GA, GB = loss_grad(QA, QB)
    if it % 10 == 0:
        bA = barrier_grad(QA); bB = barrier_grad(QB)
    GA = GA - mu*bA; GB = GB - mu*bB
    GA = hproj(QA, GA); GB = hproj(QB, GB)
    QA2 = orth((QA - eta*GA).T.conj().conj().T.T) if False else None
    mA = 0.9*mA + GA; mB = 0.9*mB + GB
    mA = hproj(QA, mA); mB = hproj(QB, mB)
    QA = np.linalg.qr(QA - eta*mA)[0]; QB = np.linalg.qr(QB - eta*mB)[0]
    QA, QB = flatten(QA, QB, 30)
    if it % 100 == 0 or it == iters-1:
        mu = max(mu*0.7, 1e-9); EPS = max(EPS*0.6, 1e-7)
        pa = plucker_logs(QA).min(); pb = plucker_logs(QB).min()
        f2, w2 = loss_grad(QA, QB, False)
        nz = (w2 < 1e-8).sum(1); gap = np.where(nz < n, w2[np.arange(N), np.minimum(nz, n-1)], 0)
        print(f"it {it} f={f2:.3e} dims(<1e-8)={nz.sum()} minGapNext={gap.min():.2e} levB=[{lev(QB).min():.3f},{lev(QB).max():.3f}] "
              f"minlogPl A={pa:.2f} B={pb:.2f} mu={mu:.1e} ({time.time()-t0:.0f}s)", flush=True)
np.save(f'sol_{prof}_s{seed}.npy', np.stack([QA, QB]))
