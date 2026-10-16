"""Ray test: for a frame direction N0 (source frames -t P_k, sink frames I - t P_k, ALL t), Z_A(tN0) = (I - tK_A) A~,
K_A = A~ N0 A~^{-1} - N0 (A~ = A (x) I_m, N0 = blockdiag(-P_k)). Z singular for every A in V at some t <=> the charpolys of
K_A share a nonzero root for all A. We compute f_A = charpoly(K_A) mod p (x-power stripped) for 3 independent random G and
their gcd. Trivial gcd => for every scaling t (any algebraic t, i.e. any denominator, via restriction of scalars) some A gives
full rank. usage: ray.py mode m fam [nG]"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, json, time, numpy as np, flint
from circ import *
from families import family, conflict_graph

mode, m, fam = sys.argv[1], int(sys.argv[2]), sys.argv[3]
nG = int(sys.argv[4]) if len(sys.argv) > 4 else 3
c = Circuit(mode=mode, seed=0); nb = conflict_graph(c); W = c.W
Min, Mout, info = family(fam, m, c, np.random.default_rng(7 * 0 + m), nb)
N0 = Min  # (W,m,m), M_in = -P = N0, M_out = I + N0 at t=1
polys = []
for g in range(nG):
    A = c.transfer(np.random.default_rng(5000 + g))
    t0 = time.time()
    Ai = np.array(flint.nmod_mat(A.tolist(), P).inv().entries(), dtype=np.int64).reshape(W, W)
    # K = A~ N0 A~^{-1} - N0 ; build A~ N0 : block (j,k) = A[j,k] N0[k];  then times A~^{-1}
    AN = (A[:, :, None, None] * N0[None, :, :, :]) % P            # (W,W,m,m)
    AN = AN.transpose(0, 2, 1, 3).reshape(W * m, W * m)
    Ait = np.kron(Ai, np.eye(m, dtype=np.int64))
    K = flint.nmod_mat(AN.tolist(), P) * flint.nmod_mat(Ait.tolist(), P)
    K = np.array(K.entries(), dtype=np.int64).reshape(W * m, W * m)
    for k in range(W):
        K[k*m:(k+1)*m, k*m:(k+1)*m] = (K[k*m:(k+1)*m, k*m:(k+1)*m] - N0[k]) % P
    f = flint.nmod_mat(K.tolist(), P).charpoly()
    co = [int(x) for x in f.coeffs()]; z = 0
    while co[z] == 0: z += 1
    f = flint.nmod_poly(co[z:], P)
    polys.append(f)
    print(json.dumps(dict(mode=mode, m=m, fam=fam, g=g, deg=f.degree(), xpow=z, sec=round(time.time() - t0, 1))), flush=True)
gg = polys[0]
for f in polys[1:]: gg = gg.gcd(f)
pair = [polys[i].gcd(polys[j]).degree() for i in range(nG) for j in range(i + 1, nG)]
rec = dict(mode=mode, m=m, fam=fam, W=W, nG=nG, gcd_deg_all=gg.degree(), gcd_deg_pairs=pair, **info)
print(json.dumps(rec), flush=True)
print('PASS_EXPQR_RAY' if gg.degree() == 0 else 'FAIL_EXPQR_RAY')   # trivial gcd: no rescaling survives GRS
