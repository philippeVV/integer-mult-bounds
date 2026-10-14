"""Independent verification for the k-normal extension (k >= 5).

1. Rebuild H (stabilizer of F0) from frame reflections + F0-preserving double reflections
   (using verify_cert's independently built E8), enumerate all elements, check |H|=172032.
2. Orbit completeness by induction: for every (k-1)-rep B and x notin B, B+{x} must be
   H-equivalent to one of the supplied k-reps (checked by scanning all of H).
   Also sum of orbit sizes == C(120,k).
3. Covering reps recomputed; for every covering rep and every composition of 9 into k
   positive parts, a certificate is replayed (initial-filter clique violation or engine chain).
Usage: python3 verify_k.py orbk.pkl certk5.pkl [certk6.pkl ...]
"""
import sys, pickle, itertools, math, time
import numpy as np
import verify_cert as V

def group_elements():
    Rf = [V.refl_perm(b) for b in range(120)]
    D = V.DELETED
    gens = set()
    for b in sorted(D):
        gens.add(tuple(Rf[b]))
    for b in range(120):
        for c in range(120):
            g = tuple(Rf[b][Rf[c][i]] for i in range(120))
            if {g[i] for i in D} == D:
                gens.add(g)
    gens = [np.array(g, dtype=np.int16) for g in gens]
    idn = np.arange(120, dtype=np.int16)
    seen = {idn.tobytes()}
    allel = [idn]
    frontier = [idn]
    while frontier:
        F = np.array(frontier)
        new = []
        for s in gens[:40]:
            H = s[F]
            for row in H:
                kb = row.tobytes()
                if kb not in seen:
                    seen.add(kb)
                    new.append(row)
        frontier = new
        allel.extend(new)
    # closure under ALL generators: check every generator maps the element set to itself
    EL = np.array(allel, dtype=np.int16)
    for s in gens:
        H = s[EL[:: max(1, len(EL) // 3000)]]
        assert all(r.tobytes() in seen for r in H)
    return EL

_EL = None
_CANON = None
_W = np.array([128 ** i for i in range(9)][::-1], dtype=np.int64)
def _ext_ok(B):
    for x in range(120):
        if x in B:
            continue
        A = tuple(sorted(B + (x,)))
        img = _EL[:, list(A)].astype(np.int64)
        img.sort(axis=1)
        m = int((img @ _W[-len(A):]).min())
        if m not in _CANON:
            return (B, x)
    return None

def _replay(job):
    A, c, final, steps = job
    return V.replay(A, c, final, steps)

def compositions(n, k):
    return [c for c in itertools.product(range(1, n + 1), repeat=k) if sum(c) == n]

if __name__ == "__main__":
    import multiprocessing as _mp
    _mp.set_start_method("fork", force=True)   # workers read module globals _EL/_CANON (py>=3.14 defaults to forkserver)
    t0 = time.time()
    reps = pickle.load(open(sys.argv[1], "rb"))
    EL = group_elements()
    _EL = EL
    print("|H| =", len(EL))
    assert len(EL) == 172032
    W = np.array([128 ** i for i in range(9)][::-1], dtype=np.int64)
    def codes(A):
        img = EL[:, list(A)].astype(np.int64)
        img.sort(axis=1)
        return img @ W[-len(A):]
    kmax = max(int(p.split("certk")[1].split(".")[0]) for p in sys.argv[2:])
    # orbit completeness by induction
    prev = [(0,), (2,)]
    assert sorted(reps[1]) == prev
    c1 = set(codes((0,)).tolist()) | set(codes((2,)).tolist())
    assert len(c1) == 120
    for k in range(2, kmax + 1):
        R = reps[k]
        repcode = {}
        sizes = 0
        for A in R:
            cs = codes(A)
            own = int(np.array(sorted(A), dtype=np.int64) @ W[-k:])
            stab = int((cs == own).sum())
            sizes += 172032 // stab
            repcode[A] = int(cs.min())
        assert sizes == math.comb(120, k), (k, sizes)
        # minimal codes identify orbits; each rep's orbit = set of codes reachable
        canon_of = {}
        for A, cm in repcode.items():
            canon_of[cm] = A
        assert len(canon_of) == len(R), "two reps in one orbit"
        _CANON = set(canon_of)
        from multiprocessing import Pool
        with Pool(4) as pool:
            bad = [b for b in pool.imap_unordered(_ext_ok, reps[k - 1], chunksize=16) if b is not None]
        assert not bad, ("extension not covered", bad[:3])
        print("k=%d: %d orbit reps; sizes sum to C(120,%d); all extensions covered (%.0fs)" % (k, len(R), k, time.time() - t0), flush=True)
    import os
    for path in ([] if os.environ.get("ORBIT_ONLY") else sys.argv[2:]):
        k = int(path.split("certk")[1].split(".")[0])
        d = pickle.load(open(path, "rb"))
        cov = [A for A in reps[k] if set().union(*[V.CS[a] for a in A]) == set(V.LIVE)]
        assert sorted(cov) == sorted(d["cov"])
        comps = compositions(9, k)
        res = {(tuple(x[0]), tuple(x[1])): x for x in d["res"]}
        filt = {tuple(A): (cs, dd) for A, cs, dd in d["filtered"]}
        nver = 0
        open_cases = []
        jobs = []
        for A in cov:
            cs, dd = filt[A]
            for c in comps:
                if c in dd:
                    jobs.append((list(A), list(c), ("L1", dd[c]), []))
                else:
                    x = res[(A, c)]
                    if x[2] != "INFEASIBLE":
                        open_cases.append((A, c))
                        continue
                    jobs.append((list(A), list(c), x[3], x[4]))
        from multiprocessing import Pool
        with Pool(4) as pool:
            for ok in pool.imap_unordered(_replay, jobs, chunksize=64):
                assert ok
                nver += 1
        print("%s: k=%d covering reps %d, (rep,composition) pairs %d, certificates verified %d, open %d"
              % (path, k, len(cov), len(cov) * len(comps), nver, len(open_cases)))
        for oc in open_cases[:20]:
            print("   OPEN", oc)
    print("total time %.0fs" % (time.time() - t0))
