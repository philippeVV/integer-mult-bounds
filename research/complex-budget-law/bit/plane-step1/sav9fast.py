"""Exact max of the 9-term savings relaxation (same model as savings9.py) by DFS with bitmasks; standard library."""
import core
def max_savings(Pp):
    Pp = list(Pp); n = len(Pp)
    comp = [0] * n   # comp[i] bitmask of j with B0[Pp[i]][Pp[j]] == 1 (includes i itself)
    for i in range(n):
        for j in range(n):
            if core.B0[Pp[i]][Pp[j]] == 1: comp[i] |= 1 << j
    best = [0]
    # state: Xm, Ym masks; degree checks: for each x in X, |comp[x] & Y| <= 1; for each y in Y, |comp[y] & X| <= 1
    def ok(Xm, Ym):
        m = Xm
        while m:
            i = (m & -m).bit_length() - 1; m &= m - 1
            if bin(comp[i] & Ym).count('1') > 1: return False
        m = Ym
        while m:
            i = (m & -m).bit_length() - 1; m &= m - 1
            if bin(comp[i] & Xm).count('1') > 1: return False
        return True
    def edges(Xm, Ym):
        e = 0; m = Xm
        while m:
            i = (m & -m).bit_length() - 1; m &= m - 1
            e += bin(comp[i] & Ym).count('1')
        return e
    def rec(i, Xm, Ym, sav):
        if sav + (n - i) <= best[0]: return
        if i == n:
            best[0] = max(best[0], sav); return
        b = 1 << i
        for addX, addY in ((1, 1), (1, 0), (0, 1), (0, 0)):
            X2 = Xm | (b if addX else 0); Y2 = Ym | (b if addY else 0)
            if (addX or addY) and not ok(X2, Y2): continue
            nX = bin(X2).count('1'); nY = bin(Y2).count('1')
            if nX + nY - edges(X2, Y2) > 9: continue
            rec(i + 1, X2, Y2, sav + addX + addY - (addX & addY))
    rec(0, 0, 0, 0)
    return best[0]
