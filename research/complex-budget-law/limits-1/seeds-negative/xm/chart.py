"""Rank-r targets via kernel charts: residual_chi = M(chi) R_chi [I_r; Z_chi] / ||M(chi)||_F  (bilinear, smooth).
usage: chart.py PROFILE seed iters    PROFILE: tri<W> | m<r>:<m> | orb<r>:<m> | list:chi=r,chi=r,..."""
import sys, itertools, json, time, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from geom import singer
N, D = singer(17); n = 18; k = 9; D = np.array(D)
prof = sys.argv[1]; seed = int(sys.argv[2]); iters = int(sys.argv[3]); rng = np.random.default_rng(seed); chis = np.arange(N)
extra = 0
if '+' in prof: prof, extra = prof.split('+'); extra = int(extra)
if prof.startswith('tri'):
    W = float(prof[3:]); dist = np.minimum(chis, N-chis); r = np.maximum(0, 9 - np.ceil(9*dist/W)).astype(int)
elif prof.startswith('m'):
    rr, m = map(int, prof[1:].split(':')); r = np.zeros(N, int); r[rng.permutation(N)[:m]] = rr
elif prof.startswith('orb'):
    rr, m = map(int, prof[3:].split(':')); r = np.zeros(N, int); seen = set(); cnt = 0
    for c in rng.permutation(N)[1:]:
        if c in seen or c == 0: continue
        o = [c, 17*c % N, 289*c % N]; seen |= set(o); r[o] = rr; cnt += 1
        if cnt == m: break
elif prof.startswith('list:'):
    r = np.zeros(N, int)
    for t in prof[5:].split(','): c, v = map(int, t.split('=')); r[c % N] = v
if extra:
    z0 = np.nonzero(r == 0)[0]; r[rng.permutation(z0)[:extra]] = 1
    prof = prof + '+' + str(extra)
act = np.nonzero(r)[0]; ra = r[act]; Lk = np.exp(-2j*np.pi*np.outer(act, D)/N)
print('target', int(r.sum()), 'blocks', len(act), 'sum r^2', int((ra**2).sum()), flush=True)
free = np.ones((k, k), bool); free[0, :] = False; free[:, 0] = False; nx = free.sum(); nv = nx + k*k
R = [np.linalg.qr(rng.normal(size=(k, k)) + 1j*rng.normal(size=(k, k)))[0] for _ in act]
zsz = [(k-x)*x for x in ra]; off = np.concatenate([[0], np.cumsum(zsz)]); nz = int(off[-1]); ntot = nv + nz
def unpack(p):
    z = p[:ntot] + 1j*p[ntot:]
    X = np.ones((k, k), complex); X[free] = z[:nx]; Y = z[nx:nv].reshape(k, k); return X, Y, z[nv:]
def Mof(X, Y, L):
    G = np.hstack([np.eye(k), X]); H = np.hstack([-Y.T, np.eye(k)]); return np.einsum('ij,mj,lj->mil', H, L, G), G, H
def resid_jac(p, want_j=True):
    X, Y, Z = unpack(p); M, G, H = Mof(X, Y, Lk)
    nr = sum(k*x for x in ra); res = np.zeros(nr, complex)
    if want_j: J = np.zeros((nr, ntot), complex)
    row = 0
    for i, c in enumerate(act):
        rr = ra[i]; Zi = Z[off[i]:off[i+1]].reshape(k-rr, rr)
        K = R[i] @ np.vstack([np.eye(rr), Zi])                  # (k, rr)
        s = np.linalg.norm(H)*np.linalg.norm(G)/k; Mi = M[i]/s
        F = Mi @ K                                              # (k, rr)
        res[row:row+k*rr] = F.ravel()
        if want_j:
            # dF = (dM/s) K + Mi R [0; dZ]  (normalisation frozen)
            # dM = H Λ [0|dX]^T + [-dY^T|0] Λ G^T ; F_{a,j} : d/dX_{bc}: (HΛ)_{a,9+c} K_{b,j}; d/dY_{ba}: -(ΛG^T)_{b,:}K_{:,j}
            HL = (H * Lk[i])[:, k:] / s                         # (k,k): [a, c]
            dX = np.einsum('ac,bj->ajbc', HL, K)                # F[a,j] wrt X[b,c]
            LG = ((Lk[i][:k, None] * G.T[:k]) @ K) / s          # (k_b, rr): row b
            dY = np.zeros((k, rr, k, k), complex)               # F[a,j] wrt Y[b,a']: -delta(a,a') LG[b,j]
            for a in range(k): dY[a, :, :, a] = -LG.T
            MR = Mi @ R[i][:, rr:]                              # (k, k-rr)
            dZ = np.einsum('ab,jl->ajbl', MR, np.eye(rr))       # F[a,j] wrt Z[b,l]
            blk = slice(row, row+k*rr)
            J[blk, :nx] = dX.reshape(k*rr, k, k)[:, free]
            J[blk, nx:nv] = dY.reshape(k*rr, k*k)
            J[blk, nv+off[i]:nv+off[i+1]] = dZ.reshape(k*rr, (k-rr)*rr)
        row += k*rr
    fr = np.concatenate([res.real, res.imag])
    if not want_j: return fr
    return fr, np.vstack([np.hstack([J.real, -J.imag]), np.hstack([J.imag, J.real])])
p = rng.normal(size=2*ntot); lam = 1e-4; f, J = resid_jac(p); c = (f**2).sum(); t0 = time.time()
for it in range(iters):
    Uj, sj, Vjh = np.linalg.svd(J, full_matrices=False)
    st = -(Vjh.T * (sj/(sj**2 + lam))) @ (Uj.T @ f)
    p2 = p + st; f2 = resid_jac(p2, False); c2 = (f2**2).sum()
    if c2 < c: p = p2; c = c2; lam = max(lam/3, 1e-16); f, J = resid_jac(p)
    else: lam *= 4
    if c < 1e-26 or lam > 1e8: break
X, Y, _ = unpack(p)
Lall = np.exp(-2j*np.pi*np.outer(chis, D)/N); M, G, H = Mof(X, Y, Lall)
s = np.linalg.svd(M, compute_uv=False); sn = s/(np.linalg.norm(H)*np.linalg.norm(G)/k); dims = (sn < 1e-10).sum(1)
nxt = min(sn[c, k-1-dims[c]] for c in range(N) if dims[c] < k)
SUB = np.array(list(itertools.combinations(range(n), k))); cg = np.ones(n); B = np.hstack([np.eye(k), Y])
for _ in range(500):
    Q, _ = np.linalg.qr((G*cg).conj().T); l = (abs(Q)**2).sum(1); cg *= np.sqrt(0.5/l)
    if abs(l-.5).max() < 1e-10: break
pl = [np.linalg.slogdet(np.linalg.qr((Mx*cg).conj().T)[0][SUB])[1].min() for Mx in (G, B)]
QB = np.linalg.qr((B*cg).conj().T)[0]; lb = (abs(QB)**2).sum(1)
print(json.dumps(dict(profile=prof, seed=seed, target=int(r.sum()), cost=float(c), it=it, S=int(dims.sum()), on_target=int(np.minimum(dims, r).sum()),
      gap=float(nxt), plA=round(float(pl[0]), 2), plB=round(float(pl[1]), 2), levB=[round(float(lb.min()), 3), round(float(lb.max()), 3)],
      maxabs=round(float(max(abs(X).max(), abs(Y).max())), 2), t=round(time.time()-t0))), flush=True)
np.save(f'chart_{prof.replace(":", "_")[:40]}_s{seed}.npy', np.stack([X, Y]))
