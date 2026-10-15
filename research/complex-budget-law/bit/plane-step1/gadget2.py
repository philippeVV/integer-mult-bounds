"""Class U, r punctured, |S| = 2 window sets (maxsav9 = 10): a cost-64 decoder (DeltaC = 1) would need a 9-term T1 of
cost 62 plus exactly one cost-1 row-anchored term (R = {a, b}, delta empty) and one cost-1 column-anchored term.
For n1 = 9 the plane block is N[S,S]_{x,x'} = pi(Z_x, Y_x') = Z_x^T QF^{-1} Y_x' (independent of T1), where Z_x / Y_x are
the functionals of row / column x on the line ports.  This script enumerates all (Z_a, Z_b) admitting a shared
correction set A of pairing cost 1 (A is forced on Req(a) u Req(b); the free region is left out, which can only lower
the cost), the same for columns, and asks whether pi(Z_x, Y_x') = delta can hold.  If no set admits it, the n1 = 9
branch of the window is closed for that set."""
import core, json, itertools, sys, time
L = list(range(111)); QFi = core._f2_inv(core.QF)
funcs = [tuple(y) for y in itertools.product([0, 1], repeat=9)]
CW = {y: core.codeword(list(y)) for y in funcs}
def pi(z, y): return sum(z[i] & QFi[i][j] & y[j] for i in range(9) for j in range(9)) % 2
def gadgets(r, a, b, Lp, side):
    """side 'row': pairs (Z_a, Z_b) with a cost-1 shared row correction; 'col' is identical by symmetry of the
    requirement and of the pairing rank, so the same list serves columns."""
    E = core.line_E(); E[a] = [core.LINES[a], core.LINES[r]]; E[b] = [core.LINES[b], core.LINES[r]]
    req = {x: {c for c in Lp if any(core.dot(u, core.LINES[c]) for u in E[x])} for x in (a, b)}
    out = []
    for Za in funcs:
        ca = CW[Za] & req[a]
        for Zb in funcs:
            if Za == Zb == tuple([0] * 9): continue
            cb = CW[Zb] & req[b]
            inter = req[a] & req[b]
            if (ca & inter) != (cb & inter): continue
            A = sorted(ca | cb)
            if not A: continue
            c = core.pair_rank([a, b], A, E)
            if c == 1: out.append((Za, Zb))
    return out
if __name__ == '__main__':
    sets = json.load(open('sets_all.json'))
    P = list(range(111, 120))
    import sav9fast
    pairs2 = {}
    for r, S, rk in sets:
        if r < 111: continue
        for c in itertools.combinations(S, 2): pairs2.setdefault(c, r)
    win = [(c, r) for c, r in pairs2.items() if sav9fast.max_savings(P + list(c)) == 10]
    print('window 2-sets', len(win), flush=True)
    step = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    t = time.time(); found = 0
    for (a, b), r in win[::step]:
        Lp = set(L) - {a, b}
        G = gadgets(r, a, b, Lp, 'row')
        ok = False
        for Za, Zb in G:
            for Ya, Yb in G:
                if pi(Za, Ya) == 1 and pi(Zb, Yb) == 1 and pi(Za, Yb) == 0 and pi(Zb, Ya) == 0:
                    ok = True; break
            if ok: break
        found += ok
        print((a, b), r, 'row gadgets', len(G), 'identity possible' if ok else 'closed', round(time.time() - t), flush=True)
    print('sets with a possible identity:', found)
