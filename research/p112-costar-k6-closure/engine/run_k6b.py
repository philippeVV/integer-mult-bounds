"""k=6 engine pass over initial-filter survivors with incremental logging (imap).
Jobs ordered by capacity pattern (4,1^5), (3,2,1^4), (2^3,1^3)."""
import os
JOBS = int(os.environ.get("JOBS", "4"))   # worker processes (original runs used 15)
import pickle, sys, time
from multiprocessing import Pool
from collections import Counter
import run_k

def work(args):
    normals, caps = args
    from engine4 import propagate
    t = time.time()
    try:
        status, fin, steps = propagate(list(normals), list(caps), ztime=20.0)
    except Exception as e:
        return (normals, caps, "ERROR", repr(e), [], time.time() - t)
    return (normals, caps, status, fin, steps, time.time() - t)

if __name__ == "__main__":
    k = 6
    reps = pickle.load(open("orbk.pkl", "rb"))[k]
    from model import COSTAR, FULL
    cov = []
    for A in reps:
        U = 0
        for a in A:
            U |= COSTAR[a]
        if U == FULL:
            cov.append(A)
    t = time.time()
    with Pool(JOBS, initializer=run_k.init) as p:
        fr = p.map(run_k.filt, cov, chunksize=8)
    surv = [(A, c) for A, cs, dd in fr for c in cs]
    order = {(4, 1, 1, 1, 1, 1): 0, (3, 2, 1, 1, 1, 1): 1, (2, 2, 2, 1, 1, 1): 2}
    surv.sort(key=lambda x: order[tuple(sorted(x[1], reverse=True))])
    print("filter done %.0fs, survivors %d" % (time.time() - t, len(surv)), flush=True)
    pickle.dump({"cov": cov, "filtered": fr}, open("certk6_filter.pkl", "wb"))
    res = []
    t = time.time()
    with Pool(JOBS) as p:
        for i, x in enumerate(p.imap(work, surv, chunksize=4)):
            res.append(x)
            if (i + 1) % 2000 == 0 or i + 1 == len(surv):
                cnt = Counter(y[2] if y[2] != "INFEASIBLE" else y[3][0] for y in res)
                pat = Counter((tuple(sorted(y[1], reverse=True)), y[2]) for y in res)
                print("%d/%d %.0fs %s" % (i + 1, len(surv), time.time() - t, dict(cnt)), flush=True)
                print("    by pattern:", dict(pat), flush=True)
                pickle.dump({"cov": cov, "filtered": fr, "res": res}, open("certk6.pkl", "wb"))
    print("done all")
