"""Inertia of the second-order code/gauge form K at two-sided flat codes of any degree (girth >= 8)."""
from model import *
import xg_core as Dc
def flatten(q, it=4000):
    d = np.zeros(q.shape[0])
    for _ in range(it):
        Yq = orth(np.exp(d)[:, None] * q); a = np.sum(abs(Yq)**2, 1); d += 0.5*(np.log(.5) - np.log(a))
    return orth(np.exp(d)[:, None] * q)
def XY_codes(V, nP, E, Qs):
    X = np.zeros((E, E)); Y = np.zeros((E, E))
    for i, (es, q) in enumerate(zip(V, Qs)):
        P = q @ q.conj().T; O = 4 * abs(P)**2; np.fill_diagonal(O, 0)
        (X if i < nP else Y)[np.ix_(es, es)] = O
    return X, Y
def inert(M, tol=1e-9):
    ev = np.linalg.eigvalsh((M+M.T)/2); s = abs(ev).max(); return (int((ev > tol*s).sum()), int((ev < -tol*s).sum()), int((abs(ev) <= tol*s).sum()))
name = sys.argv[1]; V, nP, E, g = load_any(name); rng = np.random.default_rng(2)
print(name, 'E', E, 'V', len(V), 'girth', g, 'degrees', sorted(set(len(es) for es in V)))
for kind in ['fourier', 'rand', 'rand', 'realrand']:
    if kind == 'fourier': Qs = Dc.fourier_codes(V, rng)
    else: Qs = [flatten(q) for q in Dc.rand_codes(V, rng, cplx=(kind != 'realrand'))]
    X, Y = XY_codes(V, nP, E, Qs); K = Kmat(X, Y); G, C, D = xdelta(K)
    ev, es_, _ = schur(K)
    print(f'{kind}: K inertia {inert(K)}  (expect ({E-1}, {E-len(V)+1}, {len(V)}))  Gamma inertia {inert(G)}  Schur max {es_[-1]:.2e}', flush=True)
