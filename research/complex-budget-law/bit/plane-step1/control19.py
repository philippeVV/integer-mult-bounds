"""Positive control for ten19.py: build an explicit 10-term label decoder with savings 19 (eight orthogonal-pair terms
from a nine-pair dual system, the split term (x1, y') + (x1, lambda_y1 + lambda_y')), verify sum y z^T = QF exactly and
its costs on L' = complement of P' by rational pairing rank, then check that ten19 finds a configuration on P'."""
import core, random, sys, itertools
import ten19, t1search as T
B = core.B0; G = core.G
rnd = random.Random(3)
pairs = [(x, y) for x in range(120) for y in range(120) if x != y and G[x][y] == 0]
def search():
    for attempt in range(2000):
        sysm = []
        cand = pairs[:]
        rnd.shuffle(cand)
        for x, y in cand:
            if any(x in p or y in p for p in sysm): continue
            if all(B[x][yy] == 0 and B[xx][y] == 0 for xx, yy in sysm):
                sysm.append((x, y))
                if len(sysm) == 9: return sysm
    return None
sysm = search(); print('dual system', sysm)
lam = T.lam
Y = [lam(x) for x, y in sysm]; Z = [lam(y) for x, y in sysm]
X = [[sum(Y[k][i] & Z[k][j] for k in range(9)) % 2 for j in range(9)] for i in range(9)]
assert X == core.QF
x1, y1 = sysm[0]
used = {l for p in sysm for l in p}
yp = next(y for y in range(120) if y not in used and G[x1][y] == 0 and y != x1)
zsum = [a ^ b for a, b in zip(lam(y1), lam(yp))]
terms = list(zip(Y[1:], Z[1:])) + [(Y[0], lam(yp)), (Y[0], zsum)]
X2 = [[sum(t[0][i] & t[1][j] for t in terms) % 2 for j in range(9)] for i in range(9)]
assert X2 == core.QF
Pp = sorted(used | {yp})
Lp = set(range(120)) - set(Pp)
E = core.line_E()
costs = [core.pair_rank(sorted(core.codeword(y) & Lp), sorted(core.codeword(z) & Lp), E) for y, z in terms]
print('10-term costs', costs, sum(costs), 'savings', 80 - sum(costs), '|P\'|', len(Pp))
ex, f = ten19.exists19(Pp)
print('ten19 finds a configuration:', ex)
assert ex
