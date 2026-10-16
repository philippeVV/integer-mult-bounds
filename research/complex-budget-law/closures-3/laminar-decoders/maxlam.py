"""How many of a G/2-type decoder's corrections can be direct, helper-free gates?
Greedy append-only laminar schedule on the full bidirected free graph of P120 / P112 / P104 / E7: repeatedly scan the
remaining corrections in random order and append one whenever the joined frame stays inside p_a^perp.
Lower bound on the max laminar sub-edge-set (the remainder must be relayed, i.e. needs helpers).
Upper-bound side: a greedy packing of edge-disjoint bowtie obstructions (each forces >= 1 relayed correction)."""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import random, sys, itertools, json
from sets import *

def greedy(L, seed, passes=50, strategy='rand'):
    idx = list(bits(L)); rnd = random.Random(seed)
    E = [(b, a) for a in idx for b in idx if a != b and not (BROW[a] >> b) & 1]
    xs = {b: 1 << b for b in idx}; ts = {a: 0 for a in idx}
    # forbidden label masks for target a: required/diag labels = BROW[a]
    done = []; rest = E[:]
    for p in range(passes):
        rnd.shuffle(rest)
        if strategy == 'small':
            rest.sort(key=lambda e: bin(xs[e[0]] | ts[e[1]]).count('1'))
        nxt = []
        for b, a in rest:
            U = xs[b] | ts[a]
            if U & BROW[a] & L:
                nxt.append((b, a)); continue
            xs[b] = U; ts[a] = U; done.append((b, a))
        if len(nxt) == len(rest): break
        rest = nxt
    return len(E), len(done)

def bowties(L, limit=10**6):
    """greedy edge-disjoint packing of (undirected) bowties in the free graph of L; each needs >= 1 relayed corr."""
    idx = list(bits(L)); free = lambda a, b: a != b and not (BROW[a] >> b) & 1
    used = set(); cnt = 0
    for c in idx:
        nb = [a for a in idx if free(a, c)]
        tris = [(a, b) for a, b in itertools.combinations(nb, 2) if free(a, b)]
        for (a, b), (d, e) in itertools.combinations(tris, 2):
            if len({a, b, d, e}) < 4: continue
            if any(free(x, y) for x in (a, b) for y in (d, e)): continue
            ed = {frozenset(p) for p in [(c, a), (c, b), (a, b), (c, d), (c, e), (d, e)]}
            if ed & used: continue
            used |= ed; cnt += 1
    return cnt

if __name__ == '__main__':
    P104, _ = p104(); fr = frames(); F0 = sum(1 << a for a in fr[0])
    systems = {'P120': FULL, 'P112': FULL & ~F0, 'P104': P104, 'E7': e7()}
    out = {}
    for name, L in systems.items():
        best = 0
        for s in range(int(sys.argv[1]) if len(sys.argv) > 1 else 4):
            for strat in ('rand', 'small'):
                tot, d = greedy(L, s, strategy=strat); best = max(best, d)
        bt = bowties(L)
        V = bin(L).count('1')
        out[name] = dict(V=V, corrections=tot, direct_greedy=best, frac=best / tot, per_port=best / V, bowties_disjoint=bt)
        print(name, out[name], flush=True)
    json.dump(out, open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'maxlam.json'), 'w'), indent=1)
    ok = [out[k]['bowties_disjoint'] for k in ('P120', 'P112', 'P104', 'E7')] == [465, 407, 344, 136]   # expQS table
    print('PASS_MAXLAM' if ok else 'FAIL_MAXLAM')
