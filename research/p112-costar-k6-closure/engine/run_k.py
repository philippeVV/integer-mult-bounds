"""k-normal extension: for every covering k-set orbit rep and every composition of 9
into k positive parts, apply the initial clique filter (L0+L1 at S_alpha supports),
then the full certificate engine on survivors."""
import os
JOBS = int(os.environ.get("JOBS", "4"))   # worker processes (original runs used 15)
import pickle, sys, time
import numpy as np
from multiprocessing import Pool
from collections import Counter
from model import COSTAR, FULL, bits, NL

def compositions(n, k):
    if k == 1:
        yield (n,)
        return
    for first in range(1, n - k + 2):
        for rest in compositions(n - first, k - 1):
            yield (first,) + rest

_MQ = None
_SZ = None
def init():
    global _MQ, _SZ
    cl = pickle.load(open("cliques.pkl", "rb"))
    M = np.zeros((len(cl), NL), dtype=np.float32)
    for i, q in enumerate(cl):
        for b in bits(q):
            M[i, b] = 1
    V = np.zeros((NL, 120), dtype=np.float32)
    for a in range(120):
        for b in bits(COSTAR[a]):
            V[b, a] = 1
    _MQ = (M @ V).astype(np.int8)
    _SZ = M.sum(axis=1).astype(np.int8)

_CLL = None
def filt(A):
    global _CLL
    if _CLL is None:
        _CLL = pickle.load(open("cliques.pkl", "rb"))
    k = len(A)
    MA = _MQ[:, list(A)]
    out = []
    dead = {}
    for c in compositions(9, k):
        tot = np.minimum(MA, np.array(c, dtype=np.int8)).sum(axis=1, dtype=np.int16)
        bad = np.nonzero(tot < _SZ)[0]
        if len(bad) == 0:
            out.append(c)
        else:
            dead[c] = tuple(bits(_CLL[int(bad[0])]))
    return A, out, dead

def work(args):
    normals, caps = args
    from engine4 import propagate
    t = time.time()
    status, fin, steps = propagate(list(normals), list(caps))
    return (normals, caps, status, fin, steps, time.time() - t)

if __name__ == "__main__":
    k = int(sys.argv[1])
    reps = pickle.load(open("orbk.pkl", "rb"))[k]
    cov = []
    for A in reps:
        U = 0
        for a in A:
            U |= COSTAR[a]
        if U == FULL:
            cov.append(A)
    print("k", k, "orbit reps", len(reps), "covering reps", len(cov), flush=True)
    t = time.time()
    with Pool(JOBS, initializer=init) as p:
        fr = p.map(filt, cov, chunksize=4)
    surv = [(A, c) for A, cs, dd in fr for c in cs]
    print("initial L1 filter: (rep,composition) pairs total %d, surviving %d, time %.1f"
          % (len(cov) * len(list(compositions(9, k))), len(surv), time.time() - t), flush=True)
    print("surviving capacity patterns:", Counter(tuple(sorted(c, reverse=True)) for A, c in surv))
    t = time.time()
    with Pool(JOBS) as p:
        res = p.map(work, surv, chunksize=1)
    print("engine time %.1f" % (time.time() - t))
    print(Counter(x[2] if x[2] == "SURVIVE" else x[3][0] for x in res))
    print("rule usage:", Counter(s[0] for x in res for s in x[4]))
    for x in res:
        if x[2] == "SURVIVE":
            print("SURV", x[0], x[1])
    pickle.dump({"cov": cov, "filtered": fr, "res": res}, open(f"certk{k}.pkl", "wb"))
