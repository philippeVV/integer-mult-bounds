"""Exact census of reader spans on L' = L \\ S for all 511 functionals (premise of #1084 Lemma 1 / Theorem (b) on puncture P u S).
dim >= 8 certified by rank mod a large prime (lower bound for rational rank); dim = 7 certified by an integer normal vector.
Report: min dim, and the set of reduced (dim 7) functionals, which must be exactly lambda_x for x in P u S."""
import core, json, sys, itertools, time
PR = 1000003
L = list(range(111)); P = list(range(111, 120))
LINES = core.LINES
def rank_mod(vs):
    M = [[x % PR for x in v] for v in vs]; r = 0
    for c in range(8):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]; inv = pow(M[r][c], PR - 2, PR)
        M[r] = [x * inv % PR for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(a - f * b) % PR for a, b in zip(M[i], M[r])]
        r += 1
    return r
roots_lines = list(range(120))
funcs = [y for y in itertools.product([0, 1], repeat=9) if any(y)]
CW = [core.codeword(y) for y in funcs]
lamx = {}
for x in range(120):
    c = frozenset(t for t in range(120) if core.B0[t][x]); lamx[c] = x
def census(S):
    Lp = set(L) - set(S); Pp = set(P) | set(S); red = []; mind = 9
    for cw in CW:
        R = sorted(cw & Lp)
        d = rank_mod([LINES[t] for t in R])
        if d < 8:
            # certify upper bound: find a line x orthogonal to all of R (normal); reduced sides are root hyperplanes here
            normals = [x for x in range(120) if all(core.G[x][t] == 0 for t in R)]
            assert d == 7 and normals, (d, len(R))
            red.append(cw)
        mind = min(mind, d)
    ok = all(cw in lamx and lamx[cw] in Pp for cw in red) and len(red) == len(Pp)
    return mind, len(red), ok
if __name__ == '__main__':
    sets = [s for s in json.load(open('sets_all.json')) if len(s[1]) == 5 and s[2] == 4]
    step = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    t = time.time(); bad = 0; cnt = 0
    for r, S, rk in sets[::step]:
        mind, nred, ok = census(S); cnt += 1
        if not (mind >= 7 and ok): bad += 1; print('BAD', r, S, mind, nred)
    print('checked', cnt, 'bad', bad, round(time.time() - t), 's')
