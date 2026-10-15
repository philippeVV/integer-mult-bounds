"""Generic masked LM on a port set L (subset of the 120 lines)."""
import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import coo_matrix
from e8 import *
def setup(L):
    idx = list(bits(L)); n = len(idx)
    REQ = np.array([[a == b or IP[a][b] == 0 for b in idx] for a in idx])
    return idx, n, REQ
def lm(L, Vm, Wm, seed=0, tries=4, maxfev=300, sym=False):
    idx, n, REQ = setup(L); K = Vm.shape[1]
    vid = -np.ones((n, K), int); c = 0
    for a in range(n):
        for k in range(K):
            if Vm[a, k]: vid[a, k] = c; c += 1
    if sym: wid = vid
    else:
        wid = -np.ones((n, K), int)
        for a in range(n):
            for k in range(K):
                if Wm[a, k]: wid[a, k] = c; c += 1
    nv = c; I = []; J = []; Lk = []; tg = []; e = 0
    for a in range(n):
        for b in range(n):
            if not REQ[a, b]: continue
            ks = [k for k in range(K) if vid[a, k] >= 0 and wid[b, k] >= 0]
            if not ks:
                if a == b: return None
                continue
            for k in ks: I.append(e); J.append(vid[a, k]); Lk.append(wid[b, k])
            tg.append(1.0 if a == b else 0.0); e += 1
    I = np.array(I); J = np.array(J); Lk = np.array(Lk); tg = np.array(tg); ne = e
    def f(x): return np.bincount(I, x[J] * x[Lk], minlength=ne) - tg
    rows = np.concatenate([I, I]); cols = np.concatenate([J, Lk])
    def jac(x): return coo_matrix((np.concatenate([x[Lk], x[J]]), (rows, cols)), shape=(ne, nv)).tocsr()
    rng = np.random.default_rng(seed); best = None
    for t in range(tries):
        r = least_squares(f, rng.normal(size=nv), jac=jac, method='trf', tr_solver='lsmr', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=maxfev)
        res = float(np.linalg.norm(r.fun))
        if best is None or res < best[0]: best = (res, r.x)
        if res < 1e-9: break
    x = best[1]
    V = np.zeros((n, K)); W = np.zeros((n, K))
    V[vid >= 0] = x[vid[vid >= 0]]; W[wid >= 0] = x[wid[wid >= 0]]
    return best[0], V, W
