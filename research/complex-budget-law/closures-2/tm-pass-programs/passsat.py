"""Exact minimum PASS count for interchanging two u-bit address fields (h,d) on a
k-tape machine, in the lag-pass model (strategist idea 5, 'masked-lag pass').

Layout (P, h, g, d, z): records z are opaque (one F2/linear 'symbol' per cell of
the 2^u x 2^u grid), the spectator fields P, g are long.  A pass is one sweep of
all heads over the whole volume (cost ~V, independent of the lags).  In a pass:
  * exactly one tape o is written (in place, at its own head);
  * every other tape t is read at a fixed 2D lag delta_t = (dh, dd) relative to
    o (head offset); o reads its own old cell (lag 0);
  * mask: per output cell c, the written value is
        route: one of { T_t[c+delta_t] } or 0          (select)
        xor  : any F2 combination of those             (masked XOR)
    masks may depend on (h,d) arbitrarily (finite control tracks the u+u bits).
  * reads that leave the grid in h or d (carry into P or g) are not usable.
Start: tape 0 = identity, others 0.  Goal: some tape holds the transpose.
"""
import sys, time, json, itertools
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType
from tmsat import Enc


def build(u, k, P, mode, lags=None):
    n = 1 << u
    N = n * n
    cells = [(h, d) for h in range(n) for d in range(n)]
    idx = {cd: i for i, cd in enumerate(cells)}
    if lags is None:
        lags = [(a, b) for a in range(-(n - 1), n) for b in range(-(n - 1), n)]
    e = Enc()
    T = [[[e.const(j == 0 and v == i) for v in range(N)] for i in range(N)] for j in range(k)]
    for p in range(P):
        wsel = [e.new() for _ in range(k)]
        e.cl.append(wsel)
        for a, b in itertools.combinations(wsel, 2): e.cl.append([-a, -b])
        # lag one-hot per tape (lag of tape t relative to the written tape); written tape must have lag 0
        lsel = []
        for t in range(k):
            ls = [e.new() for _ in lags]
            enc = CardEnc.equals(ls, 1, top_id=e.n, encoding=EncType.seqcounter)
            e.cl.extend(enc.clauses)
            e.n = max(e.n, enc.nv)
            z = lags.index((0, 0))
            e.cl.append([-wsel[t], ls[z]])
            lsel.append(ls)
        # sources at output cell c: R_t[c] (form) and valid_t[c]
        newval = [None] * N  # value written at cell c (computed for 'the' written tape)
        for c, (h, d) in enumerate(cells):
            srcs = []
            for t in range(k):
                terms = [[] for _ in range(N)]
                valid = []
                for li, (dh, dd) in enumerate(lags):
                    hh, ddd = h + dh, d + dd
                    if 0 <= hh < n and 0 <= ddd < n:
                        j = idx[(hh, ddd)]
                        valid.append(lsel[t][li])
                        for v in range(N):
                            terms[v].append(e.AND(lsel[t][li], T[t][j][v]))
                form = [e.OR(terms[v]) for v in range(N)]
                ok = e.OR(valid)
                coef = e.new()
                e.cl.append([-coef, ok])
                srcs.append((coef, form))
            coefs = [s[0] for s in srcs]
            if mode == 'route':
                for a, b in itertools.combinations(coefs, 2): e.cl.append([-a, -b])
                newval[c] = [e.OR([e.AND(cf, f[v]) for cf, f in srcs]) for v in range(N)]
            else:
                newval[c] = [e.XOR([e.AND(cf, f[v]) for cf, f in srcs]) for v in range(N)]
        T = [[[e.MUX(wsel[t], newval[c][v], T[t][c][v]) for v in range(N)] for c in range(N)]
             for t in range(k)]
    goal = []
    for t in range(k):
        g = e.new(); goal.append(g)
        for c, (h, d) in enumerate(cells):
            want = idx[(d, h)]
            for v in range(N):
                e.cl.append([-g, T[t][c][v] if v == want else -T[t][c][v]])
    e.cl.append(goal)
    return e


def minP(u, k, mode, Pmax=12, lags=None, log=True):
    for P in range(1, Pmax + 1):
        t0 = time.time()
        e = build(u, k, P, mode, lags)
        with Solver(name='cadical153', bootstrap_with=e.cl) as s:
            r = s.solve()
        if log: print(f"  u={u} k={k} {mode} P={P} {'SAT' if r else 'unsat'} vars={e.n} {time.time()-t0:.1f}s", flush=True)
        if r: return P
    return None


if __name__ == '__main__':
    u, k = int(sys.argv[1]), int(sys.argv[2]); mode = sys.argv[3]
    P = minP(u, k, mode)
    rec = dict(model='lagpass', u=u, k=k, mode=mode, minP=P)
    print("RESULT", json.dumps(rec), flush=True)
    open('results_pass.jsonl', 'a').write(json.dumps(rec) + "\n")
