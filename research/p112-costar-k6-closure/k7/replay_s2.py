"""Independent re-replay of all stored k=7 certificates (s2/) with verify_cert.replay.
Checks: every (rep, composition) of every covering rep appears exactly once; every non-survivor has a
certificate that replays (filter kills: ("L1", clique) rebuilt from the stored clique index; engine
kills: stored (final, steps)).  usage: replay_s2.py NWORKERS"""
import glob, lzma, pickle, sys, itertools, numpy as np, os
from multiprocessing import Pool
COMPS = [c for c in itertools.product(range(1, 10), repeat=7) if sum(c) == 9]
def loadcerts(f):
    out = []
    with lzma.open(f, "rb") as fh:
        while True:
            try: x = pickle.load(fh)
            except EOFError: break
            if isinstance(x, list): out.extend(x)
            else: out.append(x)
    return out
def job(ci):
    import verify_cert as V
    from model import bits
    cl = pickle.load(open("cliques.pkl", "rb"))
    reps = np.load("reps7.npz")["cov_A"]
    rows = np.load("s2/d%04d.npz" % ci)["rows"]
    certs = loadcerts("s2/e%04d.xz" % ci)
    eng = {(tuple(A), tuple(c)): (fin, st) for A, c, fin, st in certs}
    ok = bad = surv = 0
    seen = set()
    for ridx, ii, kind, ns, rok, qi in rows.tolist():
        A = tuple(int(a) for a in reps[ridx]); c = COMPS[ii]
        assert (ridx, ii) not in seen; seen.add((ridx, ii))
        if kind == 7: surv += 1; continue
        if kind == 0: fin, st = ("L1", tuple(bits(cl[qi]))), []
        else: fin, st = eng.pop((A, c))
        try: V.replay(list(A), list(c), fin, st); ok += 1
        except V.Fail: bad += 1
    assert not eng
    return ci, ok, bad, surv, len(seen), len(set(r for r, _ in seen))
if __name__ == "__main__":
    n = len(glob.glob("s2/d*.npz")); tot = np.zeros(5, dtype=np.int64)
    with Pool(int(sys.argv[1]), maxtasksperchild=1) as p:
        for k, r in enumerate(p.imap_unordered(job, range(n))):
            tot += np.array(r[1:]); 
            if k % 20 == 0 or k == n - 1: print(k + 1, n, "ok/bad/surv/pairs/reps", tot.tolist(), flush=True)
    print("REPLAY", "ALL OK" if tot[1] == 0 else "FAILURES", tot.tolist())
