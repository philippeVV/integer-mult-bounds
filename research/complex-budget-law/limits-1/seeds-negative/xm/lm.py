"""Gauss-Newton/LM feasibility search for Singer-equivariant complex [18,9] pairs on PG(2,17) with prescribed
block dims r_chi:  residuals = the r_chi smallest singular values of M(chi) = H Λ_chi^{-1} G^T, G=[I|X], H=[-Y^T|I].
M(chi) is invariant under the common diagonal gauge; the gauge is fixed by X[0,:]=X[:,0]=1 (145 complex unknowns).
usage: lm.py PROFILE seed [maxfev]   PROFILE: m<r>:<m> | c<r> | tri<W> | orb<r>:<m> (r on m ×17-orbits)"""
import sys, itertools, time, numpy as np, json
from scipy.optimize import least_squares
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from geom import singer
N, D = singer(17); n = 18; k = 9; E = N*n; D = np.array(D)
prof = sys.argv[1]; seed = int(sys.argv[2]); maxfev = int(sys.argv[3]) if len(sys.argv) > 3 else 200
rng = np.random.default_rng(seed); chis = np.arange(N)
orbs = []; seen = set()
for c in range(N):
    if c in seen: continue
    o = sorted({c, 17*c % N, 289*c % N}); seen |= set(o); orbs.append(o)
if prof.startswith('m'):
    rr, m = map(int, prof[1:].split(':')); r = np.zeros(N, int); r[rng.permutation(N)[:m]] = rr
elif prof.startswith('orb'):
    rr, m = map(int, prof[3:].split(':')); r = np.zeros(N, int)
    for o in [orbs[i] for i in rng.permutation(len(orbs))[:m]]: r[o] = rr
elif prof.startswith('c'): r = np.full(N, int(prof[1:]))
elif prof.startswith('tri'):
    W = int(prof[3:]); dist = np.minimum(chis, N-chis); r = np.maximum(0, 9 - np.ceil(9*dist/W)).astype(int)
act = np.nonzero(r)[0]; ra = r[act]
Linv = np.exp(-2j*np.pi*np.outer(act, D)/N)            # (m,n) diag of Λ^{-1}
free = np.ones((k, k), bool); free[0, :] = False; free[:, 0] = False
nx = free.sum(); nv = nx + k*k
def unpack(p):
    z = p[:nv] + 1j*p[nv:]
    X = np.ones((k, k), complex); X[free] = z[:nx]; Y = z[nx:].reshape(k, k)
    return X, Y
def Ms(X, Y):
    G = np.hstack([np.eye(k), X]); H = np.hstack([-Y.T, np.eye(k)])
    return np.einsum('ij,mj,lj->mil', H, Linv, G), G, H        # (m,k,k)
def resid(p):
    X, Y = unpack(p); M, _, _ = Ms(X, Y)
    s = np.linalg.svd(M, compute_uv=False)                     # descending
    return np.concatenate([s[i, k-ra[i]:] for i in range(len(act))])
def jac(p):
    X, Y = unpack(p); M, G, H = Ms(X, Y)
    U, s, Vh = np.linalg.svd(M)
    rows = []
    for i in range(len(act)):
        for j in range(k-ra[i], k):
            u = U[i, :, j]; v = Vh[i, j, :].conj()
            w = (u.conj() @ H) * Linv[i]                            # row vector length n
            cX = np.outer(v, w[k:])                                 # d/dX_{ab} = w_{9+b} v_a
            y = (Linv[i] * (G.T @ v))[:k]
            cY = -np.outer(y, u.conj())                             # d/dY_{ji} = -conj(u_i) y_j
            c = np.concatenate([cX[free], cY.ravel()])
            rows.append(np.concatenate([c.real, -c.imag]))
    return np.array(rows)
def adj(M):
    U, s, Vh = np.linalg.svd(M); ph = np.linalg.det(U) * np.linalg.det(Vh)
    pr = np.array([[np.prod(np.delete(si, j)) for j in range(k)] for si in s])
    return np.einsum('mij,mj,mjl->mil', Vh.conj().transpose(0, 2, 1), pr, U.conj().transpose(0, 2, 1)) * ph[:, None, None]
DET = all(ra == 1)
if DET:
    def resid(p):
        X, Y = unpack(p); M, _, _ = Ms(X, Y)
        nrm = np.prod(np.linalg.svd(M, compute_uv=False)[:, :k-1], axis=1)
        z = np.linalg.det(M) / nrm
        return np.concatenate([z.real, z.imag])
    def jac(p):
        X, Y = unpack(p); M, G, H = Ms(X, Y)
        nrm = np.prod(np.linalg.svd(M, compute_uv=False)[:, :k-1], axis=1)
        Ad = adj(M) / nrm[:, None, None]                 # d z = tr(Ad dM)
        # dM = H Λ^{-1} [0|dX]^T + [-dY^T|0] Λ^{-1} G^T
        # tr(Ad H Λ^{-1} [0;dX^T]) = sum_ab (Ad H Λ^{-1})_{a, 9+b} dX_{ab}
        T1 = np.einsum('mai,ij,mj->maj', Ad, H, Linv)[:, :, k:]           # (m,k,k): coeff of dX_{ab}
        T2 = -np.einsum('mj,jl,mli->mji', Linv[:, :k], G.T[:k], Ad)           # coeff of dY_{ji}: tr(Ad[-dY^T|0]Λ^{-1}G^T)
        C = np.concatenate([T1[:, free], T2.reshape(len(act), -1)], axis=1)
        return np.vstack([np.hstack([C.real, -C.imag]), np.hstack([C.imag, C.real])])
p0 = rng.normal(size=2*nv)
t0 = time.time()
class _S: pass
sol = _S(); p = p0.copy(); lam = 1e-3; f = resid(p); c = (f**2).sum(); hist = []
for it in range(maxfev):
    J = jac(p)
    # damped minimal-norm Gauss-Newton step (works for under- and over-determined systems)
    Uj, sj, Vjh = np.linalg.svd(J, full_matrices=False)
    step = -(Vjh.T * (sj/(sj**2 + lam))) @ (Uj.T @ f)
    p2 = p + step; f2 = resid(p2); c2 = (f2**2).sum()
    if c2 < c: p, f, c = p2, f2, c2; lam = max(lam/3, 1e-15)
    else: lam *= 4
    if it % 20 == 0: hist.append((it, float(c)))
    if c < 1e-26 or lam > 1e8: break
sol.x = p; sol.nfev = it+1
res = resid(sol.x); X, Y = unpack(sol.x)
# full spectrum of block dims at the end (all chi), with gap
Lall = np.exp(-2j*np.pi*np.outer(chis, D)/N)
G = np.hstack([np.eye(k), X]); H = np.hstack([-Y.T, np.eye(k)])
Mall = np.einsum('ij,mj,lj->mil', H, Lall, G)
s = np.linalg.svd(Mall, compute_uv=False); sn = s / s[:, :1]
dims = (sn < 1e-9).sum(1)
nxt = np.array([sn[c, k-1-dims[c]] if dims[c] < k else 1 for c in range(N)])
# MDS check (all minors of the non-identity part, i.e. all k-subsets) in a flattened gauge
SUB = np.array(list(itertools.combinations(range(n), k)))
def minpl(Gm):
    Q, _ = np.linalg.qr(Gm.conj().T); return np.linalg.slogdet(Q[SUB])[1].min()
Bg = np.hstack([np.eye(k), Y])
c = np.ones(n)
for _ in range(500):   # Barthe-flatten A (common gauge applied to both codes)
    Q, _ = np.linalg.qr((G*c).conj().T); l = (abs(Q)**2).sum(1); c *= np.sqrt(0.5/l)
    if abs(l-0.5).max() < 1e-10: break
G = G*c; Bg = Bg*c
QB, _ = np.linalg.qr(Bg.conj().T); levB = (abs(QB)**2).sum(1)
out = dict(profile=prof, seed=seed, target=int(r.sum()), nres=len(res), unknowns=int(nv), cost=float((res**2).sum()), maxres=float(res.max()),
           achieved=int(dims.sum()), achieved_on_target=int(np.minimum(dims, r).sum()), min_gap_next=float(nxt.min()),
           minlogPl_A=float(minpl(G)), minlogPl_B=float(minpl(Bg)), levB=[round(float(levB.min()),3), round(float(levB.max()),3)], maxabs=float(max(abs(X).max(), abs(Y).max())),
           nfev=sol.nfev, hist=hist[-4:], t=round(time.time()-t0, 1))
print(json.dumps(out), flush=True)
np.save(f"lm_{prof.replace(':', '_')}_s{seed}.npy", np.stack([X, Y]))
# limits-1 replay marker: target reached exactly on target blocks with an MDS (non-degenerate) Pluecker margin
reached = out['achieved'] == out['achieved_on_target'] == out['target'] and min(out['minlogPl_A'], out['minlogPl_B']) > -20
print('PASS_EXPXM_LM' if reached else 'NOT_REACHED (target missed or degenerate)')
