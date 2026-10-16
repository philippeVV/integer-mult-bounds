"""Catalogue of test port sets L (bitmasks over the 120 lines)."""
from subsys import *
from fractions import Fraction as Fr
R240 = ROOT + [tuple(-x for x in r) for r in ROOT]
g = (1, 2, 3, 5, 8, 13, 21, 34)   # generic-ish (check no root orthogonal)
assert all(sum(x * y for x, y in zip(r, g)) != 0 for r in R240)
POS = [r for r in R240 if sum(x * y for x, y in zip(r, g)) > 0]
def add(a, b): return tuple(x + y for x, y in zip(a, b))
PS = set(POS)
SIMPLE = [r for r in POS if not any(add(a, b) == r for a in POS for b in POS)]
assert len(SIMPLE) == 8
def solve_coeffs(v):
    # coefficients of v in simple-root basis (rational)
    M = [[Fr(SIMPLE[j][i]) for j in range(8)] + [Fr(v[i])] for i in range(8)]
    for c in range(8):
        p = next(r for r in range(c, 8) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
        M[c] = [x / M[c][c] for x in M[c]]
        for r in range(8):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][8] for i in range(8)]
HIGH = max(POS, key=lambda r: sum(x * y for x, y in zip(r, g)))
MARK = [int(c) for c in solve_coeffs(HIGH)]
# fundamental coweights (= weights, simply laced, unimodular): omega_i . alpha_j = delta (true ip) -> doubled coords dot = 4 delta
def fundamental():
    W = []
    for i in range(8):
        M = [[Fr(SIMPLE[j][c]) for c in range(8)] + [Fr(4 * (i == j))] for j in range(8)]
        for c in range(8):
            p = next(r for r in range(c, 8) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
            M[c] = [x / M[c][c] for x in M[c]]
            for r in range(8):
                if r != c and M[r][c] != 0:
                    f = M[r][c]; M[r] = [x - f * y for x, y in zip(M[r], M[c])]
        W.append(tuple(M[c][8] for c in range(8)))
    return W
OMEGA = fundamental()
def coeff_line(a):
    return solve_coeffs(ROOT[a])
COEF = [coeff_line(a) for a in range(NL)]
def bds(i):
    """Borel-de Siebenthal: lines whose alpha_i-coefficient is 0 mod mark_i (extended diagram minus node i)"""
    return sum(1 << a for a in range(NL) if COEF[a][i] % MARK[i] == 0)
def levi(drop):
    """flat / Levi subsystem: lines with zero coefficient on dropped simple roots"""
    return sum(1 << a for a in range(NL) if all(COEF[a][i] == 0 for i in drop))
def catalog():
    C = {}
    for i in range(8):
        if MARK[i] > 1: C['bds:' + name(bds(i))] = bds(i)
    for i in range(8):
        L = levi([i]); C['levi:' + name(L)] = L
    for i, j in itertools.combinations(range(8), 2):
        L = levi([i, j]); C.setdefault('levi:' + name(L), L)
    F0 = frames()[0]
    C['P120'] = FULL; C['P112'] = FULL & ~sum(1 << a for a in F0)
    return C
if __name__ == '__main__':
    print('marks', MARK)
    for k, L in catalog().items(): print(k, bin(L).count('1'))
