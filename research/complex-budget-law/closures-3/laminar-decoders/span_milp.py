"""Joint decoder-support + laminar-correction MILP in the span of the known compressing pieces.
Decoder N = sum_j lam_j * (G/2 restricted to L_j) + diag(delta)   (exact on every pair: G/2 vanishes on required pairs)
  N_aa = sum_{j ni a} lam_j + delta_a = 1
  N_ab (free) = (ip_ab/2) * sum_{j ni a,b} lam_j  must be 0 unless {a,b} in the allowed laminar pair set A_lam
Cost = sum_j C_j [lam_j != 0] + sum_a [delta_a != 0]  (C_j = gauge cost of piece j, an upper bound on its own cost;
identity terms cost 1).  Pool on P120: P120 (72), P112_F = P120 - F for frames F (64), E7_alpha for all 120 alpha (49).
A_lam variants: 'none' (N = I off-diagonal: trivially laminar) and 'greedy' (the bidirected pairs of the greedy
laminar schedule of maxlam.py, all simultaneously orderable).  scipy/HiGHS MILP, big-M |lam| <= 4."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, random, json
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from sets import *
from maxlam import greedy

def piece_pool(nF, seed):
    rnd = random.Random(seed); fr = frames(); pool = []
    pool.append(('P120', FULL, gauge_cost([ROOT[a] for a in fr[0]], FULL)))
    for fi in rnd.sample(range(len(fr)), nF):
        F = sum(1 << a for a in fr[fi]); L = FULL & ~F
        pool.append((f'P112_{fi}', L, gauge_cost([ROOT[a] for a in fr[fi]], L)))
    byline = {}
    for f in fr:
        for a in f: byline.setdefault(a, []).append(f)
    for al in range(NL):
        L = sum(1 << b for b in range(NL) if IP[b][al] == 0 and b != al)
        c = min(gauge_cost([ROOT[a] for a in f], L) for f in byline[al][:4])
        pool.append((f'E7_{al}', L, c))
    return pool

def greedy_pairs(seed=0):
    # rerun the greedy schedule and return the unordered pairs with both directions scheduled
    idx = list(range(NL)); rnd = random.Random(seed)
    E = [(b, a) for a in idx for b in idx if a != b and not (BROW[a] >> b) & 1]
    xs = {b: 1 << b for b in idx}; ts = {a: 0 for a in idx}; done = set(); rest = E[:]
    for p in range(50):
        rnd.shuffle(rest); nxt = []
        for b, a in rest:
            U = xs[b] | ts[a]
            if U & BROW[a]: nxt.append((b, a)); continue
            xs[b] = U; ts[a] = U; done.add((b, a))
        if len(nxt) == len(rest): break
        rest = nxt
    return {frozenset((a, b)) for (b, a) in done if (a, b) in done}, done

def solve(pool, allowed, tl):
    J = len(pool); nv = 2 * J + 2 * NL      # lam_j, y_j, delta_a, id_a
    il = lambda j: j; iy = lambda j: J + j; idl = lambda a: 2 * J + a; iid = lambda a: 2 * J + NL + a
    rows = []; lo = []; hi = []
    cons = lil_matrix((NL + NL * NL + 2 * J + 2 * NL, nv)); r = 0
    mem = [[(L >> a) & 1 for a in range(NL)] for _, L, _ in pool]
    for a in range(NL):
        for j in range(J):
            if mem[j][a]: cons[r, il(j)] = 1
        cons[r, idl(a)] = 1; lo.append(1); hi.append(1); r += 1
    nfree = 0
    for a in range(NL):
        for b in range(a + 1, NL):
            if (BROW[a] >> b) & 1 or frozenset((a, b)) in allowed: continue
            js = [j for j in range(J) if mem[j][a] and mem[j][b]]
            if not js: continue
            for j in js: cons[r, il(j)] = 1
            lo.append(0); hi.append(0); r += 1; nfree += 1
    Mb = 4.0
    for j in range(J):
        cons[r, il(j)] = 1; cons[r, iy(j)] = -Mb; lo.append(-np.inf); hi.append(0); r += 1
        cons[r, il(j)] = 1; cons[r, iy(j)] = Mb; lo.append(0); hi.append(np.inf); r += 1
    for a in range(NL):
        cons[r, idl(a)] = 1; cons[r, iid(a)] = -Mb; lo.append(-np.inf); hi.append(0); r += 1
        cons[r, idl(a)] = 1; cons[r, iid(a)] = Mb; lo.append(0); hi.append(np.inf); r += 1
    cons = cons[:r].tocsr()
    c = np.zeros(nv)
    for j in range(J): c[iy(j)] = pool[j][2]
    for a in range(NL): c[iid(a)] = 1
    integ = np.zeros(nv); integ[[iy(j) for j in range(J)]] = 1; integ[[iid(a) for a in range(NL)]] = 1
    lb = np.r_[-Mb * np.ones(J), np.zeros(J), -Mb * np.ones(NL), np.zeros(NL)]
    ub = np.r_[Mb * np.ones(J), np.ones(J), Mb * np.ones(NL), np.ones(NL)]
    res = milp(c, constraints=LinearConstraint(cons, lo, hi), integrality=integ, bounds=Bounds(lb, ub),
               options=dict(time_limit=tl, disp=False))
    used = [(pool[j][0], round(res.x[il(j)], 4)) for j in range(J) if res.x is not None and abs(res.x[il(j)]) > 1e-7]
    return dict(status=res.message, C=res.fun, bound=getattr(res, 'mip_dual_bound', None), pieces=used[:20],
                n_id=int(round(sum(res.x[iid(a)] for a in range(NL)))) if res.x is not None else None, free_rows=nfree)

if __name__ == '__main__':
    nF = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    pool = piece_pool(nF, 0)
    print('pool', len(pool), 'costs', sorted({p[2] for p in pool}), flush=True)
    allowed_g, done = greedy_pairs()
    print('greedy laminar schedule: directed', len(done), 'bidirected pairs', len(allowed_g), flush=True)
    out = {}
    # disjoint free triangles (bidirected, components separate => laminar): greedy packing
    rnd = random.Random(1); free = lambda a, b: not (BROW[a] >> b) & 1
    for t in range(200):
        left = list(range(NL)); rnd.shuffle(left); tris = []
        while left:
            a = left.pop(); nb = [b for b in left if free(a, b)]
            pr = next(((b, c) for i, b in enumerate(nb) for c in nb[i + 1:] if free(b, c)), None)
            if pr is None: continue
            left.remove(pr[0]); left.remove(pr[1]); tris.append((a,) + pr)
        if t == 0 or len(tris) > len(best_t): best_t = tris
    allowed_t = {frozenset(p) for tr in best_t for p in [(tr[0], tr[1]), (tr[0], tr[2]), (tr[1], tr[2])]}
    print('disjoint free triangles', len(best_t), 'pairs', len(allowed_t), flush=True)
    for name, allowed in (('none', set()), ('greedy', allowed_g), ('triangles', allowed_t)):
        r = solve(pool, allowed, tl=float(sys.argv[2]) if len(sys.argv) > 2 else 300)
        out[name] = r; print(name, json.dumps(r), flush=True)
    json.dump(out, open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'span_milp.json'), 'w'), indent=1)
    print('PASS_SPAN_MILP')
