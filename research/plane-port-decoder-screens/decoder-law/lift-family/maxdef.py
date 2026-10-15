"""Upper bound on total deficiency sum_k (delta_R(g_k) + delta_A(g'_k)) over dual bases (g, g') of K^9, given
the finite sets of deficient lines on each side (necessary conditions only => valid cost lower bound 81 - maxdef)."""
import numpy as np, itertools, math
def maxdef(XL, YL, tol=1e-7):
    """XL, YL: lists of (unit vector, delta). Returns max total delta and the choice."""
    nx, ny = len(XL), len(YL)
    P = np.array([[abs(x @ y) > tol for y, _ in YL] for x, _ in XL], bool).reshape(nx, ny)
    Xm = np.array([x for x, _ in XL]).reshape(nx, -1); Ym = np.array([y for y, _ in YL]).reshape(ny, -1)
    best = (0, None)
    for Ix in range(1 << nx):
        I = [i for i in range(nx) if Ix >> i & 1]
        if I and np.linalg.matrix_rank(Xm[I], 1e-8) < len(I): continue
        dI = sum(XL[i][1] for i in I)
        for Jy in range(1 << ny):
            J = [j for j in range(ny) if Jy >> j & 1]
            tot = dI + sum(YL[j][1] for j in J)
            if tot <= best[0]: continue
            sub = P[np.ix_(I, J)] if I and J else np.zeros((len(I), len(J)), bool)
            # nonzero-pairing graph must be a partial matching
            if (sub.sum(0) > 1).any() or (sub.sum(1) > 1).any(): continue
            m = int(sub.sum())
            if len(I) + len(J) - m > 9: continue
            if J and np.linalg.matrix_rank(Ym[J], 1e-8) < len(J): continue
            best = (tot, (I, J, m))
    return best
def lines_from(defd):
    """defd: h -> unit kernel vector. group into lines with delta = log2(#classes+1)."""
    L = []
    for h, v in defd.items():
        for l in L:
            if abs(abs(l[0] @ v) - 1) < 1e-6: l[1].append(h); break
        else: L.append([v, [h]])
    return [(v, round(math.log2(len(hs) + 1))) for v, hs in L]
def span2(hs):
    S = {0}
    for h in hs: S |= {x ^ h for x in S}
    return S
def rho_def(DR, DA):
    """exact pairing rank from Def subspaces (sets of F2^9 ints incl 0): rho = (9 - dim DR) - dim(DR^perp cap DA)."""
    dR = len(DR).bit_length() - 1
    inter = [a for a in DA if all(bin(a & r).count('1') % 2 == 0 for r in DR)]
    return (9 - dR) - (len(inter).bit_length() - 1)
def mincost(XL, YL, tol=1e-7):
    """XL/YL: lists of (unit vector, list of deficient class h's). Lower bound (necessary feasibility conditions) on the
    cost of a square decoder N = sum (X g_k)(Y g'_k)^T over dual bases, using exact rho from Def subspaces."""
    DX = [span2(hs) for _, hs in XL]; DY = [span2(hs) for _, hs in YL]
    nx, ny = len(XL), len(YL)
    Xm = np.array([x for x, _ in XL]).reshape(nx, -1); Ym = np.array([y for y, _ in YL]).reshape(ny, -1)
    P = (np.abs(Xm @ Ym.T) > tol) if nx and ny else np.zeros((nx, ny), bool)
    best = (81, None)
    for Ix in range(1 << nx):
        I = [i for i in range(nx) if Ix >> i & 1]
        if I and np.linalg.matrix_rank(Xm[I], 1e-8) < len(I): continue
        for Jy in range(1 << ny):
            J = [j for j in range(ny) if Jy >> j & 1]
            sub = P[np.ix_(I, J)] if I and J else np.zeros((len(I), len(J)), bool)
            if (sub.sum(0) > 1).any() or (sub.sum(1) > 1).any(): continue
            pos = len(I) + len(J) - int(sub.sum())
            if pos > 9: continue
            if J and np.linalg.matrix_rank(Ym[J], 1e-8) < len(J): continue
            c = 9 * (9 - pos); mi = set(); mj = set()
            for a, i in enumerate(I):
                for b, j in enumerate(J):
                    if sub[a, b]: c += rho_def(DX[i], DY[j]); mi.add(i); mj.add(j)
            c += sum(rho_def(DX[i], {0}) for i in I if i not in mi) + sum(rho_def({0}, DY[j]) for j in J if j not in mj)
            if c < best[0]: best = (c, (I, J))
    return best
def lines_h(defd):
    L = []
    for h, v in defd.items():
        for l in L:
            if abs(abs(l[0] @ v) - 1) < 1e-6: l[1].append(h); break
        else: L.append([v, [h]])
    return [(v, hs) for v, hs in L]
