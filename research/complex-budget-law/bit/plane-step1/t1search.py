"""Best 9-term label decoder on L' = L \\ S (effective puncture P u S), by reduced-side assignment + explicit completion.
Cost verified exactly by rational pairing rank. Then class-P decoder = T1 + five private plane terms {a}x{a} (cost 2 each)."""
import core, json, sys, itertools, random
L = list(range(111)); P = list(range(111, 120))
QF = core.QF
lab = [core.U[t] for t in range(120)]   # u_t in F2^9 (w.r.t. IB coords)
def f2_inv(M): return core._f2_inv(M)
def matmul(A, B): return [[sum(A[i][k] & B[k][j] for k in range(len(B))) % 2 for j in range(len(B[0]))] for i in range(len(A))]
def T(A): return [list(x) for x in zip(*A)]
def rank(vs): return core.rank2([sum(b << i for i, b in enumerate(v)) for v in vs])
# functional for line x: reader set C(y) = {t : u_t . y = 1}; lambda_x has C = {t: B0[t][x]=1} -> y = QF u_x  (u_t^T QF u_x = B0[t,x])
def lam(x): return [sum(QF[i][j] & lab[x][j] for j in range(9)) % 2 for i in range(9)]
def build(Lp, assign):
    """assign: list of 9 (xr or None, ys or None). readers y_j, sources z_j with Ymat^T Zmat = QF  <=> Zmat = Ymat^{-T} QF.
    Fill unassigned readers randomly s.t. constraints hold; return terms or None."""
    rnd = random.Random(0)
    for attempt in range(200):
        Y = []
        for xr, ys in assign:
            Y.append(lam(xr) if xr is not None else [rnd.randrange(2) for _ in range(9)])
        if rank(Y) < 9: continue
        # Z^T = Y^{-1}... sum_j y_j z_j^T = Y^T Z = QF  (rows of Y are y_j) -> Z = (Y^T)^{-1} QF
        Yi = f2_inv(T(Y)); Z = matmul(Yi, QF)
        ok = all(ys is None or Z[j] == lam(ys) for j, (xr, ys) in enumerate(assign))
        if ok: return Y, Z
    return None
def cost_terms(Lp, Y, Z):
    E = core.line_E(); tot = 0; cs = []
    for y, z in zip(Y, Z):
        R = sorted(core.codeword(y) & Lp); A = sorted(core.codeword(z) & Lp)
        c = core.pair_rank(R, A, E); cs.append(c); tot += c
    return tot, cs
if __name__ == '__main__':
    sets = [s for s in json.load(open('sets_all.json')) if len(s[1]) == 5 and s[2] == 4]
    idx = int(sys.argv[1]); r, S, rk = sets[idx]
    Lp = set(L) - set(S); Pp = P + list(S)
    # incumbent-like: readers lambda_x for x in P, sources dual
    res = build(Lp, [(x, None) for x in P]); print('incumbent T1 on Lprime', cost_terms(Lp, *res))
    json.dump(dict(r=r, S=S), open('/dev/null', 'w'))

QFi = f2_inv(QF)
def pi(z, y): return sum(z[i] & QFi[i][j] & y[j] for i in range(9) for j in range(9)) % 2
def complete(assign, tries=300):
    rnd = random.Random(1)
    JR = {j: lam(x) for j, (x, _) in enumerate(assign) if x is not None}
    JS = {k: lam(y) for k, (_, y) in enumerate(assign) if y is not None}
    # linear constraints on free y_j: pi(z_k, y_j) = delta_jk for k in JS
    for j, y in JR.items():
        for k, z in JS.items():
            if pi(z, y) != (j == k): return None
    free = [j for j in range(9) if j not in JR]
    # basis of solutions for y: rows c_k = z_k^T QFi ; constraint c_k . y = delta
    C = [[sum(z[i] & QFi[i][j] for i in range(9)) % 2 for j in range(9)] for k, z in JS.items()]
    ks = list(JS)
    sols = [y for y in itertools.product([0, 1], repeat=9)]
    for attempt in range(tries):
        Y = [None] * 9
        for j in JR: Y[j] = JR[j]
        okall = True
        for j in free:
            cand = [list(y) for y in sols if all(sum(a & b for a, b in zip(C[i], y)) % 2 == (j == ks[i]) for i in range(len(ks)))]
            if not cand: okall = False; break
            Y[j] = rnd.choice(cand)
        if not okall: return None
        if rank(Y) < 9: continue
        Z = matmul(f2_inv(T(Y)), QF)
        if all(Z[k] == JS[k] for k in JS): return Y, Z
    return None
def best_assign(Pp, S_lines_all, budget_nodes=200000):
    """backtracking over term assignments (reader line or None, source line or None); returns list of (savings, assign)."""
    Pp = list(Pp); out = []; nodes = [0]
    B0 = core.B0
    def save(x, y):
        if x is None and y is None: return 0
        if x is None or y is None: return 1
        return 2 if (x != y and core.G[x][y] == 0) else 1
    best = [0]
    def rec(j, assign, usedR, usedS, sc):
        nodes[0] += 1
        if nodes[0] > budget_nodes: return
        if j == 9:
            if sc >= best[0] - 1: out.append((sc, list(assign)))
            best[0] = max(best[0], sc); return
        if sc + 2 * (9 - j) < best[0] - 1: return
        opts = []
        for x in [None] + [p for p in Pp if p not in usedR and p not in usedS]:
            for y in [None] + [p for p in Pp if p not in usedS and (p not in usedR or p == x)]:
                if x is not None and any(B0[x][yy] != 0 for k, (xx, yy) in enumerate(assign) if yy is not None): continue
                if y is not None and any(B0[xx][y] != 0 for k, (xx, yy) in enumerate(assign) if xx is not None): continue
                if x is not None and y is not None and B0[x][y] != 1: continue
                opts.append((save(x, y), x, y))
        opts.sort(key=lambda t: -t[0])
        for s_, x, y in opts:
            # symmetry: terms ordered by (x,y) to avoid permutations: require nonincreasing savings
            if assign and s_ > save(*assign[-1]): continue
            assign.append((x, y))
            rec(j + 1, assign, usedR | ({x} - {None}), usedS | ({y} - {None}), sc + s_)
            assign.pop()
    rec(0, [], set(), set(), 0)
    out.sort(key=lambda t: -t[0])
    return out, nodes[0]
