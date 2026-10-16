#!/usr/bin/env python3
"""expRJ large-n pass: consistent uniform-hop ceiling a_H = 1 - ln(E/(r k))/ln b_min, r = min(sigma_hat, 16 - d).
No LP needed. python big.py <group> <out.jsonl>"""
import sys, json, time, math
sys.path.insert(0, '.')
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from cc2 import *
from cc2 import _finish

def diam_est(ip, ind, n, rng):
    dist = np.empty(n, np.int64); order = np.empty(n, np.int64); D = 0
    s = int(rng.integers(n))
    for _ in range(4):
        t = bfs(ip, ind, s, dist, order); far = order[t - 1]; D = max(D, dist[far]); s = far
    return D

def run(name, G, d, out, rng, bspan=(1, 2, 3, 4)):
    n, E, cap = G
    ip, ind, eid = csr(n, E, cap)
    D = diam_est(ip, ind, n, rng)
    C = float(cap.sum())
    for db in bspan:
        b0 = D - db
        if b0 < 2: continue
        t0 = time.time()
        P = far_pairs_fast(ip, ind, n, b0, 1, 400)
        if len(P) < max(8, n // 20): continue
        sig = sigma_hat(n, E, cap, P, seed=1, nvec=3, nballs=8, passes=15)
        r = min(sig, 16 - d)
        aH = a_H(r, len(P), b0, C)
        rec = dict(name=name, n=n, d=d, E=C, diam=int(D), b0=int(b0), k=len(P), sigma=round(sig, 4), rate=round(r, 4),
                   roles=round(r * len(P)), arity=round(d + r, 2), E_over_rk=round(C / (r * len(P)), 3), a_H=round(aH, 4),
                   sec=round(time.time() - t0, 1))
        print(json.dumps(rec), flush=True)
        with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')

def main():
    grp, out = sys.argv[1], sys.argv[2]
    rng = np.random.default_rng(5)
    if grp == 'rr':
        # roles = r k ~ (sigma/2) n ; keep roles <= ~1e5
        for d, emax in ((10, 15), (8, 16), (6, 16), (4, 17), (12, 14), (5, 16)):
            for e in range(8, emax + 1):
                n = 1 << e
                run(f"rr{d}", rand_regular(n, d, e), d, out, rng, bspan=(1, 2, 3))
    elif grp == 'struct':
        for q in (13, 17, 29, 37, 41):
            G = lps(5, q); run(f"LPS(5,{q})", G, 6, out, rng)
        for q in (5, 17, 29):
            G = lps(13, q); run(f"LPS(13,{q})", G, 14, out, rng)
        for m_ in (16, 32, 64, 128, 256):
            run(f"Margulis{m_}", margulis(m_), 8, out, rng)
        for e in range(8, 16):
            n1 = 1 << e
            run(f"bireg(3,6)n1={n1}", rand_biregular(n1, 3, 6, e), 6, out, rng)
            run(f"bireg(4,8)n1={n1}", rand_biregular(n1, 4, 8, e), 8, out, rng)
main()
