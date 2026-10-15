"""Relaxed exact-credit bound for square decoders of kappa-family members with |A|, |B| <= 8.
Lines (exact classification, kfam_ctrl.py): X: e9 [Def span(q_A)^perp], g_f (f in A, f odd to all of A or al_f = 0) [Def <q_f>];
Y: e9', h'_f (f in B, odd to all of B or alp_f = 0) [Def <q_f>].
Pairings: e9.e9' = 1; e9.h'_f = -2 alp_f; g_f.e9' = -2 al_f; g_f.h'_f = 2; g_f.h'_f' = ka kap (r_f.r_f') + 4 al_f alp_f'
 -> forced nonzero / forced zero / free (free = may be zero or not: relaxation in favour of the decoder).
Decoder (dual bases): a position carries <= 1 X-line and <= 1 Y-line; X-line and Y-line at different positions must pair to 0;
at the same position must pair nonzero; <= 9 positions.  credit(position) = 9 - rho_exact (rho from Def subspaces).
cost >= 81 - max credit.  Returns max credit (CP-SAT, exact integer optimum)."""
from ortools.sat.python import cp_model
from e8 import *
from maxdef import rho_def, span2
def lines(A, B, ZX, ZY):
    A = list(A); B = list(B)
    def side(K, Z):
        L = []
        if K:
            D = {h for h in range(512) if all(dotq(LAB[a], h) == 0 for a in K)} - {511}
            # 511 never lies in span(q_K)^perp unless... keep the exact F2 subspace: recompute without filtering
            D = {h for h in range(512) if all(dotq(LAB[a], h) == 0 for a in K)}
            L.append(('e9', D))
        for f in K:
            if all(IP[f][a] != 0 for a in K if a != f) or f in Z:
                L.append((f, {0, LAB[f]}))
        return L
    return side(A, ZX), side(B, ZY)
def status(x, y, ZX, ZY):
    """'nz' forced nonzero, 'z' forced zero, 'free'."""
    if x == 'e9' and y == 'e9': return 'nz'
    if x == 'e9': return 'z' if y in ZY else 'nz'
    if y == 'e9': return 'z' if x in ZX else 'nz'
    if x == y: return 'nz'
    if x in ZX or y in ZY: return 'nz'
    return 'free'
def maxcredit(A, B, ZX=(), ZY=(), want=False):
    LX, LY = lines(A, B, ZX, ZY)
    cx = [9 - rho_def(D, {0}) for _, D in LX]; cy = [9 - rho_def({0}, D) for _, D in LY]
    m = cp_model.CpModel()
    ux = [m.NewBoolVar(f'x{i}') for i in range(len(LX))]; uy = [m.NewBoolVar(f'y{j}') for j in range(len(LY))]
    M = {}; obj = []
    for i, (x, Dx) in enumerate(LX):
        for j, (y, Dy) in enumerate(LY):
            st = status(x, y, set(ZX), set(ZY))
            if st == 'z': continue
            v = m.NewBoolVar(f'm{i}_{j}'); M[i, j] = v
            m.AddImplication(v, ux[i]); m.AddImplication(v, uy[j])
            if st == 'nz': m.AddBoolOr([ux[i].Not(), uy[j].Not(), v])
            loss = cx[i] + cy[j] - (9 - rho_def(Dx, Dy))
            obj.append(-loss * v)
    for i in range(len(LX)): m.Add(sum(M[i, j] for j in range(len(LY)) if (i, j) in M) <= 1)
    for j in range(len(LY)): m.Add(sum(M[i, j] for i in range(len(LX)) if (i, j) in M) <= 1)
    m.Add(sum(ux) + sum(uy) - sum(M.values()) <= 9)
    m.Maximize(sum(c * u for c, u in zip(cx, ux)) + sum(c * u for c, u in zip(cy, uy)) + sum(obj))
    s = cp_model.CpSolver(); s.parameters.num_workers = 1
    r = s.Solve(m); assert r == cp_model.OPTIMAL
    val = int(round(s.ObjectiveValue()))
    if want:
        conf = ([LX[i][0] for i in range(len(LX)) if s.Value(ux[i])], [LY[j][0] for j in range(len(LY)) if s.Value(uy[j])],
                [(LX[i][0], LY[j][0]) for (i, j), v in M.items() if s.Value(v)])
        return val, conf
    return val
