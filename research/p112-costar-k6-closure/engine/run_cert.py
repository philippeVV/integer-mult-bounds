"""Run engine2 over orbit reps x capacity compositions; store certificates."""
import os
JOBS = int(os.environ.get("JOBS", "4"))   # worker processes (original runs used 15)
import pickle, sys, time
from multiprocessing import Pool
from collections import Counter

def compositions(n, k):
    if k == 1:
        yield (n,)
        return
    for first in range(1, n - k + 2):
        for rest in compositions(n - first, k - 1):
            yield (first,) + rest

def work(args):
    normals, caps = args
    from engine2 import propagate
    t = time.time()
    status, fin, steps = propagate(list(normals), list(caps))
    if status == "SURVIVE":
        R, C = fin
        return (normals, caps, status, (R, C), steps, time.time() - t)
    return (normals, caps, status, fin, steps, time.time() - t)

if __name__ == "__main__":
    which = sys.argv[1]
    reps = pickle.load(open(f"orb{which}.pkl", "rb"))
    k = int(which)
    jobs = [(r, caps) for r, n in reps for caps in compositions(9, k)]
    print("jobs", len(jobs), flush=True)
    t = time.time()
    with Pool(JOBS) as p:
        res = p.map(work, jobs, chunksize=2)
    print("time %.1f" % (time.time() - t))
    print(Counter(x[2] if x[2] == "SURVIVE" else x[3][0] for x in res))
    print("kill-steps used: max %d, total %d" % (max(len(x[4]) for x in res), sum(len(x[4]) for x in res)))
    print("rule usage:", Counter(s[0] for x in res for s in x[4]))
    print("cases needing kills:", sum(1 for x in res if x[4]))
    for x in res:
        if x[2] == "SURVIVE":
            print("SURV", x[0], x[1])
    pickle.dump(res, open(f"cert{which}.pkl", "wb"))
