"""Close |S| = 3, 4 in class P: the only surviving case is n1 = 10 label terms on L' = L \\ S with
  T1 cost <= 60 (|S| = 4)  i.e. ten cost-6 terms (lambda_x, lambda_y), x _|_ y, x != y, x, y in P' = P u S
  T1 cost <= 61 (|S| = 3)  i.e. nine such terms plus one term with a reduced reader lambda_x (x in P') or reduced source.
Over label coordinates: sum_j y_j z_j^T = QF with lambda_x = QF u_x, so sum_j u_{x_j} u_{y_j}^T = QF^{-1}.
Exact F2 linear algebra: enumerate the affine solution space of the pair-indicator vector c (and the free factor w)."""
import core, json, itertools, sys, time
QF = core.QF; QFi = core._f2_inv(QF)
U = core.U
def outer(u, v): return sum((u[i] & v[j]) << (9 * i + j) for i in range(9) for j in range(9))
TARGET = sum(QFi[i][j] << (9 * i + j) for i in range(9) for j in range(9))
def solve_affine(cols, target):
    """cols: list of 81-bit ints. returns (particular solution as bitmask over cols or None, nullspace basis bitmasks)."""
    piv = {}  # pivot bit -> (vec, combo)
    null = []
    for idx, v in enumerate(cols):
        combo = 1 << idx
        while v:
            h = v.bit_length() - 1
            if h in piv: pv, pc = piv[h]; v ^= pv; combo ^= pc
            else: piv[h] = (v, combo); break
        if not v: null.append(combo)
    v = target; combo = 0
    while v:
        h = v.bit_length() - 1
        if h not in piv: return None, null
        pv, pc = piv[h]; v ^= pv; combo ^= pc
    return combo, null
def min_weight_hit(part, null, npairs, want, maxdim=26):
    """does part + span(null) contain a vector whose weight on the first npairs coordinates equals want (and, for
    extra coordinates, nonzero w)? brute force over the nullspace (dimension must be <= maxdim)."""
    d = len(null)
    if d > maxdim: return None
    mask = (1 << npairs) - 1
    # gray code enumeration
    cur = part; hits = []
    for g in range(1 << d):
        if g:
            b = (g & -g).bit_length() - 1; cur ^= null[b]
        if bin(cur & mask).count('1') == want(cur): return cur
    return False
def check_set(S, s):
    Pp = list(range(111, 120)) + list(S)
    pairs = [(x, y) for x in Pp for y in Pp if x != y and core.G[x][y] == 0]
    cols = [outer(U[x], U[y]) for x, y in pairs]
    npairs = len(pairs)
    if s == 4:
        part, null = solve_affine(cols, TARGET)
        if part is None: return 'nosol', 0
        res = min_weight_hit(part, null, npairs, lambda c: 10)
        return ('FOUND' if res else ('undecided' if res is None else 'none')), len(null)
    # s == 3: nine cost-6 pairs plus one cost-7 term: reader-reduced (u_x w^T, w free) or source-reduced (w u_y^T)
    worst = 0
    for x in Pp:
        for side in (0, 1):
            extra = [outer(U[x], [int(i == k) for i in range(9)]) if side == 0 else outer([int(i == k) for i in range(9)], U[x]) for k in range(9)]
            part, null = solve_affine(cols + extra, TARGET)
            worst = max(worst, len(null))
            if part is None: continue
            wmask = ((1 << 9) - 1) << npairs
            def want(c):  # nine pairs and w != 0, or ten pairs and w == 0 (pure cost-6 case)
                return 9 if c & wmask else 10
            res = min_weight_hit(part, null, npairs, want)
            if res is None: return 'undecided', len(null)
            if res: return 'FOUND', len(null)
    return 'none', worst
if __name__ == '__main__':
    s = int(sys.argv[1])
    sets = json.load(open('sets_all.json'))
    allS = set()
    for r, S, rk in sets:
        for c in itertools.combinations(S, s): allS.add(c)
    allS = sorted(allS)
    if len(sys.argv) > 2: allS = allS[::int(sys.argv[2])]
    t = time.time(); from collections import Counter; cnt = Counter(); dims = Counter()
    for S in allS:
        res, d = check_set(S, s); cnt[res] += 1; dims[d] += 1
        if res in ('FOUND', 'undecided'): print(res, S, d, flush=True)
    print('s', s, 'sets', len(allS), dict(cnt), 'nullity dist', dict(dims), round(time.time() - t), 's')
