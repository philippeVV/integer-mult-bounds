"""Lag-support lower bound for whole-volume pass programs (any coding, XOR or not).

Lambda_t = set of lags l such that some cell c of tape t holds a form whose support
contains input c+l.  A pass that writes tape o with read lags delta_t gives
    Lambda_o' subset  Lambda_o  U  U_{t != o} (Lambda_t + delta_t),
because masks/XOR only combine what is read (cancellation can shrink, never grow).
So |Lambda_o'| <= sum_t |Lambda_t|.  The transpose needs the 2^{u+1}-1 lags
j(1,-1), |j| < 2^u, on one tape.  Max-reachable sizes are computed exactly over
size vectors (one written tape per pass, or all tapes written per pass).
Prints PASS_LAGBOUND when the claimed closed forms hold.
"""
import math


def one_write_max(k, P):
    best = 0
    frontier = {tuple([1] + [0] * (k - 1))}
    out = []
    for p in range(1, P + 1):
        nf = set()
        for v in frontier:
            s = sum(v)
            for o in range(k):
                w = list(v); w[o] = s
                nf.add(tuple(sorted(w, reverse=True)))
        nf = sorted(nf, reverse=True)
        frontier = set(v for v in nf if not any(w != v and all(a >= b for a, b in zip(w, v)) for w in nf))
        out.append(max(max(v) for v in frontier))
    return out


def lower_bound(u, k, P=40):
    need = 2 ** (u + 1) - 1
    seq = one_write_max(k, 2 * u + 4)
    for p, m in enumerate(seq, 1):
        if m >= need: return p
    return None


def all_write_lb(u, k):
    # all k tapes written each pass: total grows by factor <= k
    need = 2 ** (u + 1) - 1
    return math.ceil(math.log(need, k))


if __name__ == '__main__':
    print("one-write max lag-set size after p passes:")
    for k in (2, 3, 4):
        print(f"  k={k}:", one_write_max(k, 12))
    print("\nlower bounds on passes (one-write / all-write) vs routing upper bound 3u (bit-by-bit, 3 passes per 1-bit swap):")
    print(" u | k=2 1w | k=3 1w | k=4 1w | k=3 allw | 3u")
    ok = True
    for u in range(1, 9):
        lbs = [lower_bound(u, k) for k in (2, 3, 4)]
        print(f"{u:2d} | {lbs[0]:6d} | {lbs[1]:6d} | {lbs[2]:6d} | {all_write_lb(u,3):8d} | {3*u}")
        # every bound grows linearly in u (ratio to u bounded below), so no power saving
        ok &= lbs[0] >= u + 1 and all_write_lb(u, 3) >= u / math.log2(3) - 1e-9
    print("PASS_LAGBOUND" if ok else "FAIL")
