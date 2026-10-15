"""k=7 stage 1 (resumable).  Every 7-set is H-equivalent to B+{x}, B one of the 32,292 6-orbit reps
(orbk.pkl[6]); stage 1 enumerates all (B, x), x notin B:
  non-covering A -> canonical code recorded (orbit bookkeeping only; L0 at full supports kills all comps)
  covering A     -> initial L1 filter (ext7.filt_fast) over all 28 compositions; canonical code,
                    A itself and the 28-bit mask of surviving compositions recorded (first A per code).
Completeness of the orbit list is checked later by summing orbit sizes to C(120,7).
Output: s1/cNNNN.npz (+ .done marker).  usage: stage1.py NWORKERS   (env ORBK: path of orbk.pkl,
default ../engine/orbk.pkl; group_elems.npy and cliques.pkl are read from the working directory)"""
import os, sys, time, pickle
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool
import ext7
COMPS = ext7.COMPS
CI = {c: i for i, c in enumerate(COMPS)}
W = [128 ** i for i in range(7)][::-1]
def code(t): return sum(a * w for a, w in zip(t, W))
def init():
    ext7.init()
def job(ci):
    import canonfast as cf
    reps = pickle.load(open(ORBK, "rb"))[6][ci * CH:(ci + 1) * CH]
    cov = {}; non = set(); cnt = np.zeros(4, dtype=np.int64)   # noncov ext, cov ext, L1 kills, surv pairs
    for B in reps:
        for x, A, surv, dead in ext7.filt_fast(B):
            if A is None:
                cnt[0] += 1
                non.add(code(cf.canon(tuple(sorted(B + (x,))))))
                continue
            cnt[1] += 1; cnt[2] += len(dead); cnt[3] += len(surv)
            assert len(surv) + len(dead) == 28
            k = code(cf.canon(A))
            m = 0
            for c in surv: m |= 1 << CI[c]
            if k in cov:
                assert bin(cov[k][1]).count("1") == len(surv), ("orbit-invariance of filter violated", A)
            else:
                cov[k] = (A, m)
    ks = sorted(cov)
    np.savez_compressed("s1/c%04d.npz" % ci, cov_code=np.array(ks, dtype=np.int64),
                        cov_A=np.array([cov[k][0] for k in ks], dtype=np.uint8).reshape(-1, 7),
                        cov_mask=np.array([cov[k][1] for k in ks], dtype=np.int32),
                        non_code=np.array(sorted(non), dtype=np.int64), cnt=cnt)
    open("s1/c%04d.done" % ci, "w").write(" ".join(map(str, cnt)))
    return ci, cnt
CH = 100
ORBK = os.environ.get("ORBK", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "engine", "orbk.pkl"))
if __name__ == "__main__":
    os.makedirs("s1", exist_ok=True)
    n = (32292 + CH - 1) // CH
    todo = [i for i in range(n) if not os.path.exists("s1/c%04d.done" % i)]
    t0 = time.time(); agg = np.zeros(4, dtype=np.int64)
    with Pool(int(sys.argv[1]), initializer=init) as p:
        for k, (ci, cnt) in enumerate(p.imap_unordered(job, todo)):
            agg += cnt
            if k % 10 == 0 or k == len(todo) - 1:
                print(k + 1, len(todo), "noncov/cov/L1kills/surv", agg.tolist(), "%.0fs" % (time.time() - t0), flush=True)
    print("stage1 done", "%.0fs" % (time.time() - t0))
