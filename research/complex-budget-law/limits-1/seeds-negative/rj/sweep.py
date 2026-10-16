#!/usr/bin/env python3
"""expRJ sweep: python sweep.py <group> <out.jsonl>.  groups: rr, rrbig, struct, tanner"""
import sys, json, time, math
sys.path.insert(0, '.')
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from cc2 import *
from cc2 import _finish

def pg2(q):
    """incidence graph of PG(2,q), q prime: points 0..N-1, lines N..2N-1."""
    pts = []
    for x in range(q):
        for y in range(q):
            pts.append((x, y, 1))
    for x in range(q): pts.append((x, 1, 0))
    pts.append((1, 0, 0))
    N = len(pts); edges = []
    for j, l in enumerate(pts):
        for i, p in enumerate(pts):
            if (l[0] * p[0] + l[1] * p[1] + l[2] * p[2]) % q == 0: edges.append((i, N + j))
    return _finish(2 * N, edges)

def tanner_sessions(n, E, cap, f, b0, seed):
    """choose session edges R (fraction f of edges), keep only those whose endpoints are >= b0 apart in G - R."""
    rng = np.random.default_rng(seed)
    m = len(E); isR = np.zeros(m, bool)
    deg = np.zeros(n); np.add.at(deg, E[:, 0], 1); np.add.at(deg, E[:, 1], 1)
    quota = np.ceil(f * deg).astype(int); used = np.zeros(n, int)
    for j in rng.permutation(m):   # balanced: every gate hosts <= ceil(f deg) session ports
        u, v = E[j]
        if used[u] < quota[u] and used[v] < quota[v]:
            isR[j] = True; used[u] += 1; used[v] += 1
    for _ in range(20):
        Ei = E[~isR]; ci = cap[~isR]
        ip, ind, eid = csr(n, Ei, ci)
        dist = np.empty(n, np.int64); order = np.empty(n, np.int64)
        bad = []
        for j in np.flatnonzero(isR):
            bfs(ip, ind, E[j, 0], dist, order)
            if dist[E[j, 1]] < b0: bad.append(j)
        if not bad: break
        isR[rng.choice(bad, size=max(1, len(bad) // 2), replace=False)] = False
    R = np.flatnonzero(isR)
    Ei = E[~isR]; ci = cap[~isR]
    ip, ind, eid = csr(n, Ei, ci)
    dist = np.empty(n, np.int64); order = np.empty(n, np.int64)
    P = []
    for j in R:
        bfs(ip, ind, E[j, 0], dist, order)
        if dist[E[j, 1]] >= b0: P.append((E[j, 0], E[j, 1], dist[E[j, 1]]))
    # a Tanner session needs one port at each end; a vertex may host several sessions
    return Ei, ci, np.array(P, dtype=np.int64).reshape(-1, 3)

def far(name, n, E, cap, b0, d, seed, kmax=512, lp=True, rate_cap=None):
    ip, ind, eid = csr(n, E, cap)
    P = far_pairs(ip, ind, n, b0, seed)
    if len(P) < 4: return None
    if len(P) > kmax:
        P = P[np.random.default_rng(seed).choice(len(P), kmax, replace=False)]
    rc = (16 - d) if rate_cap is None else rate_cap
    return evaluate(name, n, E, cap, P, b0, rate_cap=rc, arity_base=lambda r: d + r, lp=lp, seed=seed)

def emit(out, r, t0):
    if r is None: return
    r = {k: (float(v) if isinstance(v, (np.floating,)) else v) for k, v in r.items()}
    r['sec'] = round(time.time() - t0, 1)
    print(json.dumps(r), flush=True)
    with open(out, 'a') as f: f.write(json.dumps(r) + '\n')

def main():
    grp, out = sys.argv[1], sys.argv[2]
    if grp == 'rr':
        for d in (3, 4, 6, 8, 10, 12):
            for n in (32, 64, 128, 256, 512):
                if n * d % 2: continue
                nn, E, cap = rand_regular(n, d, 7 + n)
                for b0 in (2, 3, 4, 5):
                    t0 = time.time()
                    try: emit(out, far(f"rr{d}", nn, E, cap, b0, d, 1), t0)
                    except Exception as e: print("ERR", d, n, b0, e, flush=True)
    elif grp == 'rrbig':
        for d, n in ((6, 1024), (4, 1024), (8, 1024), (3, 1024), (10, 1024)):
            nn, E, cap = rand_regular(n, d, 7 + n)
            for b0 in (2, 3):
                t0 = time.time()
                try: emit(out, far(f"rr{d}", nn, E, cap, b0, d, 1), t0)
                except Exception as e: print("ERR", d, n, b0, e, flush=True)
    elif grp == 'struct':
        G = lps(13, 5)
        for b0 in (2, 3): t0 = time.time(); emit(out, far("LPS(13,5)", *G, b0, 14, 1), t0)
        for m_ in (8, 12, 16, 22):
            G = margulis(m_)
            for b0 in (2, 3, 4): t0 = time.time(); emit(out, far(f"Margulis{m_}", *G, b0, 8, 1), t0)
        for (n1, t, D) in ((128, 3, 6), (256, 3, 6), (128, 4, 8), (256, 4, 8), (192, 4, 6), (384, 4, 6), (128, 6, 12), (256, 8, 8), (512, 3, 3)):
            G = rand_biregular(n1, t, D, 3)
            for b0 in (2, 3, 4): t0 = time.time(); emit(out, far(f"bireg({t},{D})n1={n1}", *G, b0, max(t, D), 1), t0)
        for q in (5, 7, 11, 13):
            G = pg2(q)
            for b0 in (2, 3): t0 = time.time(); emit(out, far(f"PG(2,{q})-inc", *G, b0, q + 1, 1), t0)
        G = lps(5, 13)
        for b0 in (3, 4, 5, 6): t0 = time.time(); emit(out, far("LPS(5,13)", *G, b0, 6, 1, kmax=256), t0)
    elif grp == 'lp2':
        for d, n, bl in ((4, 256, (2, 3)), (6, 128, (2,)), (8, 256, (2, 3)), (10, 256, (2,)), (12, 256, (2,)), (4, 512, (2,)), (6, 512, (2,)), (3, 1024, (2,))):
            nn, E, cap = rand_regular(n, d, 7 + n)
            for b0 in bl:
                t0 = time.time()
                try: emit(out, far(f"rr{d}", nn, E, cap, b0, d, 1), t0)
                except Exception as e: print("ERR", d, n, b0, e, flush=True)
    elif grp == 'tanner':
        base = [("PG(2,5)", pg2(5), 6), ("PG(2,7)", pg2(7), 8), ("PG(2,11)", pg2(11), 12), ("PG(2,13)", pg2(13), 14),
                ("LPS(13,5)", lps(13, 5), 14), ("rr6n128", rand_regular(128, 6, 5), 6), ("rr8n96", rand_regular(96, 8, 5), 8),
                ("rr16n64", rand_regular(64, 16, 5), 16), ("Margulis12", margulis(12), 8)]
        for nm, G, d in base:
            n, E, cap = G
            ip, ind, eid = csr(n, E, cap)
            dist = np.empty(n, np.int64); order = np.empty(n, np.int64)
            for f in (0.18, 0.25, 0.33):
                for b0 in (3, 5):
                    Ei, ci, P = tanner_sessions(n, E, cap, f, b0, 1)
                    if len(P) < 4 or len(P) > 512: continue
                    t0 = time.time()
                    try:
                        r = evaluate(f"Tanner:{nm}", n, Ei, ci, P, b0, rate_cap=1.0, arity_base=lambda r: d, lp=True)
                        r['f_target'] = f; r['S_over_E'] = round(len(P) / (len(P) + ci.sum()), 4)
                        # repo uniform form a = 1 - ln((E-S)/S)/ln b (expQK), assumes routing distance-limited
                        r['a_expQK'] = round(1 - math.log(ci.sum() / len(P)) / math.log(b0), 4)
                        emit(out, r, t0)
                    except Exception as e: print("ERR", nm, f, b0, e, flush=True)

if __name__ == '__main__':
    main()
