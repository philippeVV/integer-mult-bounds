"""Continuation for rank-1 targets: add chars in steps, re-solve from previous MDS solution, monitor MDS margin.
usage: inc.py seed step [order: rand|orb|ap]"""
import sys, itertools, json, time, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from geom import singer
N, D = singer(17); n = 18; k = 9; D = np.array(D)
seed = int(sys.argv[1]); step = int(sys.argv[2]); order = sys.argv[3] if len(sys.argv) > 3 else 'rand'
rng = np.random.default_rng(seed)
if order == 'rand': seq = list(rng.permutation(N))
elif order == 'orb':
    seen = set(); seq = []
    for c in rng.permutation(N):
        if c in seen: continue
        o = [c, 17*c % N, 289*c % N]; seen |= set(o); seq += o
elif order == 'ap': seq = [(i * int(sys.argv[4])) % N for i in range(N)]
free = np.ones((k, k), bool); free[0, :] = False; free[:, 0] = False
nx = free.sum(); nv = nx + k*k
SUB = np.array(list(itertools.combinations(range(n), k)))
def unpack(p):
    z = p[:nv] + 1j*p[nv:]; X = np.ones((k, k), complex); X[free] = z[:nx]; return X, z[nx:].reshape(k, k)
def GH(X, Y): return np.hstack([np.eye(k), X]), np.hstack([-Y.T, np.eye(k)])
def Ms(X, Y, Linv):
    G, H = GH(X, Y); return np.einsum('ij,mj,lj->mil', H, Linv, G), G, H
def adj(M):
    U, s, Vh = np.linalg.svd(M); ph = np.linalg.det(U) * np.linalg.det(Vh)
    pr = np.stack([np.prod(np.delete(s, j, axis=1), axis=1) for j in range(k)], 1)
    return np.einsum('mji,mj,mlj->mil', Vh.conj(), pr, U.conj()) * ph[:, None, None]
def resid(p, Linv):
    X, Y = unpack(p); M, _, _ = Ms(X, Y, Linv)
    s = np.linalg.svd(M, compute_uv=False); z = np.linalg.det(M)/np.prod(s[:, :k-1], 1)
    return np.concatenate([z.real, z.imag])
def jac(p, Linv):
    X, Y = unpack(p); M, G, H = Ms(X, Y, Linv)
    s = np.linalg.svd(M, compute_uv=False); Ad = adj(M)/np.prod(s[:, :k-1], 1)[:, None, None]
    T1 = np.einsum('mai,ij,mj->maj', Ad, H, Linv)[:, :, k:]
    T2 = -np.einsum('mj,jl,mli->mji', Linv[:, :k], G.T[:k], Ad)
    C = np.concatenate([T1[:, free], T2.reshape(len(Linv), -1)], axis=1)
    return np.vstack([np.hstack([C.real, -C.imag]), np.hstack([C.imag, C.real])])
def solve(p, Linv, iters=60):
    lam = 1e-6; f = resid(p, Linv); c = (f**2).sum()
    for it in range(iters):
        J = jac(p, Linv); Uj, sj, Vjh = np.linalg.svd(J, full_matrices=False)
        st = -(Vjh.T * (sj/(sj**2 + lam))) @ (Uj.T @ f)
        p2 = p + st; f2 = resid(p2, Linv); c2 = (f2**2).sum()
        if c2 < c: p, f, c = p2, f2, c2; lam = max(lam/3, 1e-16)
        else: lam *= 4
        if c < 1e-28 or lam > 1e6: break
    return p, c
def flatpl(X, Y):
    G, _ = GH(X, Y); B = np.hstack([np.eye(k), Y]); c = np.ones(n)
    for _ in range(300):
        Q, _ = np.linalg.qr((G*c).conj().T); l = (abs(Q)**2).sum(1); c *= np.sqrt(0.5/l)
        if abs(l-.5).max() < 1e-10: break
    out = []
    for Mx in (G*c, B*c):
        Q, _ = np.linalg.qr(Mx.conj().T); out.append(np.linalg.slogdet(Q[SUB])[1].min())
    QB, _ = np.linalg.qr((B*c).conj().T); lb = (abs(QB)**2).sum(1)
    return out, (lb.min(), lb.max())
def alldims(X, Y):
    Lall = np.exp(-2j*np.pi*np.outer(np.arange(N), D)/N); M, _, _ = Ms(X, Y, Lall)
    s = np.linalg.svd(M, compute_uv=False); sn = s/s[:, :1]; return (sn < 1e-9).sum(1), sn
p = rng.normal(size=2*nv); m = 0; best = None; t0 = time.time()
while m < N:
    m2 = min(N, m + step); act = np.array(seq[:m2]); Linv = np.exp(-2j*np.pi*np.outer(act, D)/N)
    p2, c = solve(p.copy(), Linv)
    X, Y = unpack(p2); dims, sn = alldims(X, Y); (pa, pb), lb = flatpl(X, Y)
    ok = c < 1e-24 and min(pa, pb) > -25
    print(json.dumps(dict(m=int(m2), cost=float(c), S=int(dims.sum()), on=int((dims[act] >= 1).sum()), plA=round(float(pa), 2), plB=round(float(pb), 2),
                          levB=[round(float(x), 3) for x in lb], ok=bool(ok), t=round(time.time()-t0))), flush=True)
    if not ok: break
    p = p2; m = m2; best = (m, int(dims.sum()))
print('LAST_OK', best)
