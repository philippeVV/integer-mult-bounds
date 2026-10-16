"""Fallback test: linear Schur gauge rule x = L delta, L = -Gamma^+ C. Maximise g(delta) = f(L d, (L+I) d) over delta (gradient ascent, many scales).
usage: lin_adv.py GRAPH base seed nstarts"""
import sys, time
from xr_tools import *
name, base, seed, ns = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
S, Q0, W, G, C, D, lmin, rng = setup(name, base, seed); E = S.E
print(f'{name} base {base} seed {seed}: lambda_min(Gamma) {lmin:.3e}', flush=True)
L = -np.linalg.pinv(G, rcond=1e-10) @ C
def g(d):
    th, gu, gv = fuv(S, Q0, L @ d, L @ d + d); return th, L.T @ (gu + gv) + gv
best = 0
for k in range(ns):
    sc = [0.1, 0.5, 1, 2, 4][k % 5]
    d0 = rng.standard_normal(E); d0 *= sc * np.sqrt(E) / np.linalg.norm(d0)
    r = minimize(lambda d: tuple(-np.asarray(z) for z in g(d)), d0, jac=True, method='L-BFGS-B', options=dict(maxiter=400))
    val = -r.fun; best = max(best, val)
    d = r.x; inf = inner(S, Q0, d, [L @ d, -d / 2, np.zeros(E)])[0]
    print(f'start {k} scale {sc}: max g = {val:.8f}  |d|_rms {np.linalg.norm(d)/np.sqrt(E):.3f}  true inf {inf:.8f} {"EXCEEDS 1/8" if val > 0.125 + 1e-9 else ""}', flush=True)
    np.save(f'lin_{name}_{base}_{seed}_{k}.npy', d)
print(f'SUMMARY lin {name} base {base}: max g = {best:.9f}')
