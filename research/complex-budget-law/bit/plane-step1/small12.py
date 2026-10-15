"""Close |S| = 1, 2 in class P (n1 = 10 residual): need a 10-term label decoder on L \\ S with savings >= 17 (s=1) / >= 18 (s=2).
Terms: 'pairs' (lambda_x, lambda_y), x _|_ y in P' (cost 6); extras: R_x (reduced reader, cost 7: u_x w^T), S_y (reduced
source: w u_y^T), G_u (generic, cost 8: u w^T with u enumerated).  Sum must equal QF^{-1} (label coordinates u).
Case lists by savings; exact F2 enumeration of each affine solution space (nullity cap -> 'undecided')."""
import core, json, itertools, sys, time
from small_s import outer, TARGET, solve_affine
U = core.U
E9 = [[int(i == k) for i in range(9)] for k in range(9)]
def blocks_for(kind, v):
    if kind == 'R': return [outer(U[v], e) for e in E9]
    if kind == 'S': return [outer(e, U[v]) for e in E9]
    if kind == 'G': return [outer(v, e) for e in E9]
def hit(part, null, npairs, nextra, maxdim=24):
    d = len(null)
    if d > maxdim: return None
    mask = (1 << npairs) - 1; cur = part
    wm = [((1 << 9) - 1) << (npairs + 9 * i) for i in range(nextra)]
    for g in range(1 << d):
        if g:
            b = (g & -g).bit_length() - 1; cur ^= null[b]
        nz = sum(1 for m in wm if cur & m)
        if bin(cur & mask).count('1') + nz == 10 and nz == nextra: return cur
    return False
def check(S, s):
    Pp = list(range(111, 120)) + list(S)
    pairs = [(x, y) for x in Pp for y in Pp if x != y and core.G[x][y] == 0]
    cols = [outer(U[x], U[y]) for x, y in pairs]; npairs = len(pairs)
    RS = [('R', x) for x in Pp] + [('S', x) for x in Pp]
    Gs = [('G', list(u)) for u in itertools.product([0, 1], repeat=9) if any(u)]
    cases = [[]] + [[g] for g in Gs] + [list(c) for c in itertools.combinations(RS, 2)]
    if s == 1:
        cases += [[a, g] for a in RS for g in Gs] + [list(c) for c in itertools.combinations(RS, 3)]
    und = 0
    for ex in cases:
        allc = cols[:]
        for kind, v in ex: allc += blocks_for(kind, v)
        part, null = solve_affine(allc, TARGET)
        if part is None: continue
        res = hit(part, null, npairs, len(ex))
        if res is None: und += 1
        elif res: return 'FOUND', ex
    return ('undecided' if und else 'none'), und
if __name__ == '__main__':
    s = int(sys.argv[1]); step = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    sets = json.load(open('sets_all.json'))
    allS = sorted({c for r, S, rk in sets for c in itertools.combinations(S, s)})[::step]
    t = time.time(); from collections import Counter; cnt = Counter()
    for S in allS:
        res, info = check(S, s); cnt[res] += 1
        if res != 'none': print(res, S, info, flush=True)
    print('s', s, 'sets', len(allS), dict(cnt), round(time.time() - t), 's', flush=True)
