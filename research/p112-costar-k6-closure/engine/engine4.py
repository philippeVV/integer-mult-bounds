"""engine3 + forced-zero ("Z") cliques.

Forced zero: M[r,s] = 0 is forced if r is orthogonal to s (masked) or if no group h has
r in R[h] and s in C[h] (every term vanishes there).  A Z-clique X is a set of live
lines with M[x,y] forced zero for all ordered pairs x != y in X; then M[X,X] = 4 I_X.
  L1Z : sum_g min(c_g, |X & R_g|, |X & C_g|) >= |X|   (else infeasible)
  L2GZ: Gamma(X) = groups with X&R_g and X&C_g nonempty; if sum_{Gamma} c = |X| then
        D_Gamma[X], T_Gamma[X] invertible; r notin X with M[r,x] forced zero and (r,x) not
        touched by any group outside Gamma (all x in X) -> D_g[r]=0 for g in Gamma; symm.
The Z rules are applied only when engine3's rules stall (they are slower).
Certificate records: final ("L1Z", X); step ("L2GZ", gamma, side, r, X).
"""
import time
from model import NL, ONB, COSTAR, popc, bits, LIVE
from engine2 import owners_tables
import engine3

def zgraph(R, C):
    rowg, colg = owners_tables(R, C)
    ZN = [0] * NL
    for r in range(NL):
        m = ONB[r]
        for s in range(NL):
            if s != r and not (m >> s & 1):
                if rowg[r] & colg[s] == 0 and rowg[s] & colg[r] == 0:
                    m |= 1 << s
        ZN[r] = m
    return ZN, rowg, colg

def zscan(R, C, caps, time_limit=120.0):
    """enumerate Z-cliques; return ('L1Z', X) on violation, else list of tight cliques"""
    G = len(R)
    ZN, rowg, colg = zgraph(R, C)
    tight = []
    t0 = time.time()
    res = [None]
    def bound(X):
        tot = 0
        for g in range(G):
            a = (X & R[g]).bit_count()
            if a == 0:
                continue
            b = (X & C[g]).bit_count()
            tot += min(caps[g], a, b)
        return tot
    def ext(X, cand, size):
        if res[0] is not None:
            return
        bd = bound(X)
        if bd < size:
            res[0] = X
            return
        # bound(X) is monotone in X, so a superset X' can only be violating/tight if
        # |X'| >= bound(X); |X'| <= size + |cand|.
        if size + cand.bit_count() < bd:
            return
        if bd == size and size >= 2:
            # tight under the per-group rank bound; record if capacity-tight (L2GZ premise)
            gam = [g for g in range(G) if X & R[g] and X & C[g]]
            if sum(caps[g] for g in gam) == size:
                tight.append((X, tuple(gam)))
        c = cand
        while c:
            b = c & -c
            i = b.bit_length() - 1
            c ^= b
            ext(X | b, c & ZN[i], size + 1)
            if time.time() - t0 > time_limit:
                raise TimeoutError
    for a in range(NL):
        ext(1 << a, ZN[a] & ~((1 << (a + 1)) - 1), 1)
        if res[0] is not None:
            return ("L1Z", tuple(bits(res[0]))), None, (ZN, rowg, colg)
    return None, tight, (ZN, rowg, colg)

def propagate(normals, caps, maxiter=200, ztime=120.0, R0=None, C0=None):
    G = len(normals)
    steps_all = []
    R = list(R0) if R0 is not None else [COSTAR[a] for a in normals]
    C = list(C0) if C0 is not None else [COSTAR[a] for a in normals]
    for outer in range(maxiter):
        st, fin, steps = engine3_from(normals, caps, R, C)
        steps_all += steps
        if st == "INFEASIBLE":
            return "INFEASIBLE", fin, steps_all
        R, C = [list(x) for x in fin]
        try:
            viol, tight, (ZN, rowg, colg) = zscan(R, C, caps, ztime)
        except TimeoutError:
            return "SURVIVE", (R, C), steps_all
        if viol is not None:
            return "INFEASIBLE", viol, steps_all
        newR, newC = list(R), list(C)
        new_steps = []
        for X, gam in tight:
            xs = bits(X)
            gm = sum(1 << g for g in gam)
            commonZ = (1 << NL) - 1
            for x in xs:
                commonZ &= ZN[x] | 0
            # ZN is symmetric forced zero; for a kill we need M[r,x] forced zero (ordered) -
            # use the stronger symmetric relation for simplicity.
            colX = 0
            rowX = 0
            for x in xs:
                colX |= colg[x]
                rowX |= rowg[x]
            out_col = colX & ~gm
            out_row = rowX & ~gm
            for r in bits(commonZ & ~X):
                if rowg[r] & out_col == 0:
                    hit = [g for g in gam if newR[g] >> r & 1]
                    if hit:
                        for g in hit:
                            newR[g] &= ~(1 << r)
                        new_steps.append(("L2GZ", gam, "R", r, tuple(xs)))
                if colg[r] & out_row == 0:
                    hit = [g for g in gam if newC[g] >> r & 1]
                    if hit:
                        for g in hit:
                            newC[g] &= ~(1 << r)
                        new_steps.append(("L2GZ", gam, "C", r, tuple(xs)))
        if not new_steps:
            return "SURVIVE", (R, C), steps_all
        steps_all += new_steps
        R, C = newR, newC
    return "SURVIVE", (R, C), steps_all

def engine3_from(normals, caps, R0, C0):
    return engine3.propagate(list(normals), list(caps), R0=R0, C0=C0)
