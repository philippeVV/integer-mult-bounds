"""Linear Schur gauge rule on spheres: max over |d|_rms = r of g(d) = f(L d, (L+I) d) (projected gradient ascent, multistart).
usage: lin_sphere.py GRAPH base seed nstarts r1 r2 ..."""
import sys
from xr_tools import *
name, base, seed, ns = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]); rs = [float(z) for z in sys.argv[5:]]
S, Q0, W, G, C, D, lmin, rng = setup(name, base, seed); E = S.E
L = -np.linalg.pinv(G, rcond=1e-10) @ C
Bv = np.zeros((E, len(S.V)))
for i, es in enumerate(S.V): Bv[es, i] = 1
Ub, sb, _ = np.linalg.svd(Bv, full_matrices=False); Ub = Ub[:, sb > 1e-9]
proj = lambda d: d - Ub @ (Ub.T @ d)
print(f'{name} base {base}: lambda_min(Gamma) {lmin:.3e}', flush=True)
def g(d):
    th, gu, gv = fuv(S, Q0, L @ d, L @ d + d); return th, L.T @ (gu + gv) + gv
for r in rs:
    nrm = lambda d: proj(d) * r * np.sqrt(E) / np.linalg.norm(proj(d))
    best = -1
    for k in range(ns):
        d = nrm(rng.standard_normal(E)); v, gr = g(d); eta = 0.3
        for it in range(150):
            gq = proj(gr); gt = gq - (gq @ d) / (d @ d) * d
            if np.linalg.norm(gt) < 1e-13: break
            dn = nrm(d + eta * r * np.sqrt(E) * gt / np.linalg.norm(gt)); vn, grn = g(dn)
            if vn > v: d, v, gr = dn, vn, grn; eta = min(1.0, eta * 1.3)
            else: eta *= 0.5
            if eta < 1e-7: break
        if v > best: best = v; np.save(f'linsph_{name}_{base}_{seed}_r{r:g}.npy', d)
    dbest = np.load(f'linsph_{name}_{base}_{seed}_r{r:g}.npy'); inf = inner(S, Q0, dbest, [L @ dbest, -dbest / 2, np.zeros(E)])[0]
    import xg_core as Dc
    Ys = [orth(np.exp((L @ dbest + (dbest if i >= S.nP else 0))[es])[:, None] * q) for i, (es, q) in enumerate(zip(S.V, Q0))]
    PP, PL = Dc.projectors(S.V, S.nP, E, Ys); m = Dc.moments(PP, PL, 3); th_dense = Dc.markov_ref(m, 3)
    print(f'   dense Theta3 (independent, xg_core) = {th_dense}, m = {np.round(m,6)}, true inf_x at that code = {inf:.9f}', flush=True)
    print(f'r = {r}: max over sphere of linear-rule Theta3 = {best:.9f} {"> 1/8" if best > 0.125 + 1e-10 else "<= 1/8"}', flush=True)
