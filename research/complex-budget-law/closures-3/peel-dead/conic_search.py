"""Abstract (A_Z) search for many-deficient configurations with k extra climbs.

Climbs: E_1..E_k (integer m x m, I not in their Q-span, Q-independent) and main climbs
D = I + sum_i P_i E_i, P in Z^k (forced by A_Z + telescoping, README sec. 2).
A_Z  <=>  every (k+1)-multiset of climbs Z-spans I.
Lemma R (README sec. 3) caps the number of non-identity main climbs at nmax (4 for k=2).
We maximise  delta = sum_D (m - rank D) - R  over admissible sets of deficient points P in a box.
usage: conic_search.py k m nsamples seed [box] [mode]
"""
import sys, random, itertools
from fractions import Fraction

def rank(M):
    M = [[Fraction(x) for x in r] for r in M]; n = len(M); c = len(M[0]); r = 0
    for j in range(c):
        p = next((i for i in range(r, n) if M[i][j] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(n):
            if i != r and M[i][j] != 0:
                f = M[i][j] / M[r][j]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r

def hnf_contains(vecs, target):
    """Is target in the Z-span of vecs (integer vectors)?  Integer row reduction."""
    rows = [list(v) for v in vecs if any(v)]
    n = len(target); basis = []
    col = 0
    while rows and col < n:
        nz = [r for r in rows if r[col] != 0]
        z = [r for r in rows if r[col] == 0]
        while len(nz) > 1:
            nz.sort(key=lambda r: abs(r[col]))
            p = nz[0]
            new = [p]
            for r in nz[1:]:
                q = r[col] // p[col]
                r2 = [a - q * b for a, b in zip(r, p)]
                (new if r2[col] != 0 else z).append(r2)
            nz = new
        if nz:
            basis.append(nz[0])
        rows = [r for r in z if any(r)]
        col += 1
    t = list(target)
    for b in basis:
        j = next(i for i, x in enumerate(b) if x != 0)
        if t[j] % b[j]: return False
        q = t[j] // b[j]
        t = [a - q * c for a, c in zip(t, b)]
    return not any(t)

def det_int(M):
    M = [list(r) for r in M]; n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k] != 0), None)
            if sw is None: return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def matmix(E, P, m):
    return [[(1 if i == j else 0) + sum(p * e[i][j] for p, e in zip(P, E)) for j in range(m)] for i in range(m)]

def run(k, m, nsamp, seed, box=4, mode='rank1', pbox=4, nmax=None):
    rng = random.Random(seed)
    nmax = nmax or {1: 2, 2: 4, 3: 7}[k]
    target = (1,) + (0,) * k
    best = (-10**9, None)
    hist = {}
    for it in range(nsamp):
        if mode == 'rank1':
            E = []
            for i in range(k):
                u = [rng.randint(-box, box) for _ in range(m)]; v = [rng.randint(-box, box) for _ in range(m)]
                E.append([[u[a] * v[b] for b in range(m)] for a in range(m)])
        else:
            E = [[[rng.randint(-box, box) for _ in range(m)] for _ in range(m)] for _ in range(k)]
        vecsQ = [[1 if i == j else 0 for i in range(m) for j in range(m)]] + [[x for r in e for x in r] for e in E]
        if rank(vecsQ) < k + 1: continue
        R = sum(rank(e) for e in E)
        # deficient points
        pts = []
        for P in itertools.product(range(-pbox, pbox + 1), repeat=k):
            if not any(P): continue
            Mx = matmix(E, P, m)
            if det_int(Mx) != 0: continue
            d = m - rank(Mx)
            if d > 0: pts.append((P, d))
        evec = [tuple(1 if j == i + 1 else 0 for j in range(k + 1)) for i in range(k)]
        def ok(chosen):  # chosen: list of P; check all (k+1)-multisets drawn from E's + chosen containing newest
            climbs = evec + [(1,) + tuple(P) for P in chosen]
            new = len(climbs) - 1
            for sub in itertools.combinations(range(len(climbs) - 1), k):
                if not hnf_contains([climbs[i] for i in sub] + [climbs[new]], target):
                    return False
            return True
        # also require singletons with E's ok
        cand = [(P, d) for P, d in pts if ok([P])]
        bestset = []; bestw = 0
        def dfs(start, chosen, w):
            nonlocal bestset, bestw
            if w > bestw: bestw, bestset = w, list(chosen)
            if len(chosen) == nmax: return
            for idx in range(start, len(cand)):
                P, d = cand[idx]
                if ok([c for c, _ in chosen] + [P]) if False else ok_seq(chosen, P):
                    chosen.append((P, d)); dfs(idx, chosen, w + d); chosen.pop()
        def ok_seq(chosen, P):
            return ok([c for c, _ in chosen] + [P])
        dfs(0, [], 0)
        delta = bestw - R
        key = (R, bestw)
        hist[key] = hist.get(key, 0) + 1
        if delta - 2 * R > best[0] - 2 * (best[1][1] if best[1] else 0) or best[1] is None or delta / R > best[0] / best[1][1]:
            if best[1] is None or delta / R > best[0] / best[1][1]:
                best = (delta, (E, R, bestset, len(pts), len(cand)))
    print(f"k={k} m={m} mode={mode} samples={nsamp}: hist (R, max sum d) -> count:")
    for key in sorted(hist): print("   ", key, hist[key])
    d, (E, R, S, npts, ncand) = best
    print(f"  best delta/R = {d}/{R} = {d/R:.3f}; deficient pts {npts}, admissible singles {ncand}; set {S}; E={E}")
    return d, R

if __name__ == '__main__':
    k, m, n, seed = map(int, sys.argv[1:5])
    box = int(sys.argv[5]) if len(sys.argv) > 5 else 4
    mode = sys.argv[6] if len(sys.argv) > 6 else 'rank1'
    pbox = int(sys.argv[7]) if len(sys.argv) > 7 else 4
    nmax = int(sys.argv[8]) if len(sys.argv) > 8 else None
    d, R = run(k, m, n, seed, box, mode, pbox, nmax)
    print('PASS_EXPQX_CONIC' if d <= 0 else 'FAIL_EXPQX_CONIC')   # expQX: delta_max = 0 in every sample (stop test: delta <= 2R)
