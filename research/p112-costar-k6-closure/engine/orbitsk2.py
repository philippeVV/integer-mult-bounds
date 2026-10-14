"""Orbit reps of k-subsets of lines under H via extension + fast canonical form (parallel)."""
import os
JOBS = int(os.environ.get("JOBS", "4"))   # worker processes (original runs used 15)
import pickle, sys, time
from multiprocessing import Pool
from model import COSTAR, FULL

def ext_canon(B):
    import canonfast as cf
    out = set()
    for x in range(120):
        if x not in B:
            out.add(cf.canon(tuple(sorted(B + (x,)))))
    return out

def stab(A):
    import canonfast as cf
    return cf.stab_size(A)

if __name__ == "__main__":
    kmax = int(sys.argv[1])
    try:
        reps = pickle.load(open("orbk.pkl", "rb"))
    except Exception:
        reps = {}
    if 1 not in reps:
        reps[1] = [(0,), (2,)]
    with Pool(JOBS) as p:
        for k in range(2, kmax + 1):
            if k in reps:
                print(k, "cached", len(reps[k]))
                continue
            t = time.time()
            new = set()
            for s in p.imap_unordered(ext_canon, reps[k - 1], chunksize=8):
                new |= s
            reps[k] = sorted(new)
            st = p.map(stab, reps[k], chunksize=16)
            tot = sum(172032 // s for s in st)
            cov = []
            totcov = 0
            for A, s in zip(reps[k], st):
                U = 0
                for a in A:
                    U |= COSTAR[a]
                if U == FULL:
                    cov.append(A)
                    totcov += 172032 // s
            print(k, "orbits", len(reps[k]), "sets", tot, "covering orbits", len(cov), "covering sets", totcov,
                  "time %.1f" % (time.time() - t), flush=True)
            pickle.dump(reps, open("orbk.pkl", "wb"))
