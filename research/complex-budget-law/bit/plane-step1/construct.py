"""Explicit class-P plane decoders (upper bounds), verified entry by entry and priced exactly:
 T1 = 9-term label decoder on L' = L \\ S found by reduced-side assignment (t1search), S excluded from all T1 supports,
 plus |S| private plane terms {a} x {a} (cost 2 each)."""
import core, gsat, t1search as T, json, sys
sets = [s for s in json.load(open('sets_all.json')) if len(s[1]) == 5 and s[2] == 4]
out = {}
for idx in map(int, sys.argv[1:]):
    r, S, rk = sets[idx]; E = gsat.plane_E(r, S); Lp = set(gsat.L) - set(S)
    cands, nodes = T.best_assign(T.P + list(S), S, budget_nodes=200000)
    best = None
    for sc, assign in cands[:200]:
        res = T.complete(assign)
        if res is None: continue
        tot, cs = T.cost_terms(Lp, *res)
        if best is None or tot < best[0]: best = (tot, res, assign)
        if tot <= 72 - 2 * len(S) + 0: break
    tot, (Y, Z), assign = best
    terms = [(sorted(core.codeword(y) & Lp), sorted(core.codeword(z) & Lp)) for y, z in zip(Y, Z)]
    terms += [([a], [a]) for a in S]
    ok, costs = gsat.check(E, terms)
    # control: same decoder fails the plane target if S is put back into the label supports (incumbent-style)
    print('set', idx, 'r', r, 'S', S, 'T1 cost', tot, 'verified', ok, 'costs', costs, 'total', sum(costs), flush=True)
    out[idx] = dict(r=r, S=S, terms=terms, costs=costs, total=sum(costs))
json.dump(out, open('construct_classP.json', 'w'))
