"""E_n roots in I_{n,1} (form diag(-1,1^n)), k = (3;1^n): x.x = 2, x.k = 0  <=>  sum_i x_i = 3 x0, sum x_i^2 = x0^2+2.
One representative per line (x0 > 0, or x0 = 0 and first nonzero > 0). Heights 0..Hmax."""
import numpy as np


def en_roots(n, Hmax, vals=(-2, -1, 0, 1, 2, 3)):
    out = []
    for x0 in range(Hmax + 1):
        S, Q = 3 * x0, x0 * x0 + 2
        cur = [0] * n
        vmax = max(vals)

        def rec(i, s, q):
            rem = n - i
            if rem == 0:
                if s == S and q == Q:
                    out.append([x0] + cur.copy())
                return
            # prune: remaining sum must be reachable with remaining square budget: |S-s| <= sqrt(rem*(Q-q))
            d = S - s
            if d * d > rem * (Q - q):
                return
            for v in vals:
                if q + v * v <= Q:
                    cur[i] = v
                    rec(i + 1, s + v, q + v * v)
            cur[i] = 0
        rec(0, 0, 0)
    X = np.array(out, dtype=np.int64)
    keep = []
    for r in X:
        if r[0] > 0:
            keep.append(r)
        else:
            nz = r[np.nonzero(r)[0][0]]
            if nz > 0:
                keep.append(r)
    return np.array(keep, dtype=np.int64)


if __name__ == '__main__':
    import sys
    import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
    from lib import evaluate
    X = en_roots(8, 3)
    eps = np.array([-1] + [1] * 8)
    print('E8 lines', len(X), evaluate(X, eps))
    for n in (9, 10, 12, 16):
        for H in (2, 3):
            print(n, H, len(en_roots(n, H)))
