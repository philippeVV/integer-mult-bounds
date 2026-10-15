"""Exact exclusion of 10-term label decoders with savings >= 19 on an effective puncture P' (T1 cost <= 61).

From M = I + v w^T (see sav10.py): savings >= 19 forces nine cost-6 pair terms (lambda_x, lambda_y), x _|_ y, and one
more term t with saving >= 1; a pair term j has B0[x_j, y_j] = 1, so j is not in V n W; |V n W| is odd, hence
V n W = {t}.  Pair terms get types 00, 10 (in V only), 01 (in W only).  Conditions (j != k):
    B0[x_j, y_k] = [k in V][j in W],
distinct readers / sources among pair terms (a repeated reader forces v = e_j + e_j', impossible with t in V unless
the repeat involves t), and for t: B0[x_t, y_t] = 0 if both sides are reduced.  Exhaustive DFS; standard library."""
import core
from itertools import product


def exists19(Pp):
    Pp = list(Pp)
    B = core.B0
    pairs = [(x, y) for x in Pp for y in Pp if x != y and core.G[x][y] == 0]
    TYPES = ((0, 0), (1, 0), (0, 1))          # (in V, in W)
    cand = [(tp, x, y) for tp in TYPES for (x, y) in pairs]
    cand.sort()
    def compat(p, q):
        (vp, wp), xp, yp = p; (vq, wq), xq, yq = q
        if xp == xq or yp == yq: return False
        if B[xp][yq] != (vq & wp): return False
        if B[xq][yp] != (vp & wq): return False
        return True
    N = len(cand)
    comp = [[compat(cand[i], cand[j]) for j in range(N)] for i in range(N)]
    opts = Pp + [None]
    def t_ok(chosen):
        # term t has type (1,1); W = {t} u {01 pairs}, V = {t} u {10 pairs}
        Vp = [c for c in chosen if c[0][0] == 1]; Wp = [c for c in chosen if c[0][1] == 1]
        for xt, yt in product(opts, opts):
            if xt is None and yt is None: continue
            if xt is not None and yt is not None and B[xt][yt] != 0: continue
            ok = True
            for (vj, wj), xj, yj in chosen:
                if xt is not None:
                    # reader of t (t in W) vs source of j: B0 = [j in V]
                    if B[xt][yj] != vj: ok = False; break
                    if xt == xj and not (len(Vp) == 1 and vj == 1): ok = False; break
                if yt is not None:
                    # reader of j vs source of t (t in V): B0 = [j in W]
                    if B[xj][yt] != wj: ok = False; break
                    if yt == yj and not (len(Wp) == 1 and wj == 1): ok = False; break
            if ok: return True
        return False
    found = []
    def rec(start, chosen, idxs):
        if len(chosen) == 9:
            if t_ok(chosen): found.append(list(chosen)); return True
            return False
        for i in range(start, N):
            if all(comp[i][j] for j in idxs):
                chosen.append(cand[i]); idxs.append(i)
                if rec(i + 1, chosen, idxs): return True
                chosen.pop(); idxs.pop()
        return False
    rec(0, [], [])
    return bool(found), found


if __name__ == '__main__':
    import json, sys, time
    P = list(range(111, 120))
    sets = [s for s in json.load(open('sets_all.json')) if s[0] >= 111 and len(s[1]) == 5]
    step = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    t0 = time.time(); bad = 0; cnt = 0
    for r, S, rk in sets[::step]:
        ex, f = exists19(P + S); cnt += 1
        if ex: bad += 1; print('FOUND', r, S, f[0], flush=True)
    print('checked', cnt, 'found', bad, round(time.time() - t0), 's', flush=True)
