"""k=7 stage 2 (resumable): every covering 7-orbit rep x all 28 compositions.
  composition killed by the initial L1 filter -> certificate ("L1", clique), no steps
  otherwise engineC (= engine4 rules, C clique scan) -> certificate (final, steps) or SURVIVE
Every certificate is replayed in-process by verify_cert.replay (independent of model/engine code).
Also records |Stab_H(A)| for the orbit-size sum check.
Output per chunk: s2/dNNNN.npz (per-pair table), s2/eNNNN.xz (engine certificates), s2/dNNNN.done
usage: stage2.py NWORKERS [CH]"""
import os, sys, time, pickle, lzma
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool
KIND = {"FILTER_L1": 0, "L0": 1, "L1": 2, "L1Z": 3, "L3odd": 4, "L2rank": 5, "L3GZodd": 6, "SURVIVE": 7, "ERROR": 8}
REPS = os.environ.get("REPS", "reps7.npz"); OUT = os.environ.get("OUT", "s2")
CH = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
def init():
    global ext7, engineC, V, EL, cl, bits, COMPS
    import ext7, engineC, verify_cert as V
    from model import bits
    ext7.init(); COMPS = ext7.COMPS
    EL = np.load("group_elems.npy").astype(np.int16)
    cl = pickle.load(open("cliques.pkl", "rb"))
def job(ci):
    z = np.load(REPS); As = z["cov_A"][ci * CH:(ci + 1) * CH].tolist(); Ms = z["cov_mask"][ci * CH:(ci + 1) * CH].tolist()
    rows = []; surv = []; stabs = []; t0 = time.time()
    cf = lzma.open(OUT + "/e%04d.xz.tmp" % ci, "wb", preset=6)   # stream of pickles (A, c, final, steps)
    for j, (A, m) in enumerate(zip(As, Ms)):
        A = tuple(A); ridx = ci * CH + j
        inA = np.zeros(120, dtype=bool); inA[list(A)] = True
        stabs.append(int(inA[EL[:, list(A)]].all(axis=1).sum()))
        sv, dead = ext7.filt(A)
        assert sum(1 << ext7.COMPS.index(c) for c in sv) == m, ("filter mismatch", A)
        for ii, c in enumerate(COMPS):
            if c in dead:
                fin = ("L1", tuple(bits(cl[dead[c]])))
                try: ok = V.replay(list(A), list(c), fin, [])
                except V.Fail: ok = False
                rows.append((ridx, ii, 0, 0, int(ok), dead[c]))
                continue
            try:
                st, fin, steps = engineC.propagate(list(A), list(c), ztime=20.0)
            except Exception as e:
                rows.append((ridx, ii, KIND["ERROR"], 0, 0, -1)); surv.append((A, c, "ERROR " + repr(e))); continue
            if st == "INFEASIBLE":
                try: ok = V.replay(list(A), list(c), fin, steps)
                except V.Fail as e: ok = False
                rows.append((ridx, ii, KIND[fin[0]], len(steps), int(ok), -1))
                pickle.dump((A, c, fin, steps), cf)
            else:
                rows.append((ridx, ii, KIND["SURVIVE"], len(steps), 0, -1))
                surv.append((A, c, fin))
    cf.close()
    os.replace(OUT + "/e%04d.xz.tmp" % ci, OUT + "/e%04d.xz" % ci)
    np.savez_compressed(OUT + "/d%04d.npz" % ci, rows=np.array(rows, dtype=np.int32).reshape(-1, 6),
                        stab=np.array(stabs, dtype=np.int32))
    pickle.dump(surv, open(OUT + "/v%04d.pkl" % ci, "wb"))
    open(OUT + "/d%04d.done" % ci, "w").write("%.0f" % (time.time() - t0))
    r = np.array(rows).reshape(-1, 6)
    return ci, np.bincount(r[:, 2], minlength=9), int((r[:, 4] == 0).sum() - (r[:, 2] == 7).sum() - (r[:, 2] == 8).sum()), len(surv)
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    n = (len(np.load(REPS)["cov_code"]) + CH - 1) // CH
    todo = [i for i in range(n) if not os.path.exists(OUT + "/d%04d.done" % i)]
    t0 = time.time(); agg = np.zeros(9, dtype=np.int64); bad = 0; ns = 0
    names = sorted(KIND, key=KIND.get)
    with Pool(int(sys.argv[1]), initializer=init, maxtasksperchild=1) as p:
        for k, (ci, bc, nbad, nsv) in enumerate(p.imap_unordered(job, todo)):
            agg += bc; bad += nbad; ns += nsv
            if k % 5 == 0 or k == len(todo) - 1:
                print(k + 1, len(todo), dict(zip(names, agg.tolist())), "replayFAIL", bad, "surv", ns, "%.0fs" % (time.time() - t0), flush=True)
    print("stage2 done this session", "%.0fs" % (time.time() - t0))
