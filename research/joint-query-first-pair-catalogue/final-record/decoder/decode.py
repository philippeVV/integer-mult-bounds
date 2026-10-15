"""Explicit ordered decoders for the QQ16 selected inputs from the FINAL contract (+ retained records).

Notation (chart s): B_s z = R_s iota_U z, K_s = K_(s,U), P_sA, P_sT the selected sections,
ell_sA, ell_sT the selected queries (CONSERVED-QQ16.md).  Every product is applied in the written
order; no inverse of an ordinary query, of a chart route, or of any Schur complement is used.
The only scalar inverses are 1/4 (inside the sections) and, for banks 1-3, 1/11.

Inputs (observations):
  Xf[k], Yf[k]  final selected X/Y of bank k            C[k]  final Y_kU (completed-X descendant)
  f[s]          forward histories, s in {0,2,4,6,8}      INIT[(k,s)] = B_s X_kU^initial (INIT64)
  dA[k]         Y_kA right after bank k's first backward call (retained record)
  XB[k]         optional cross-bank initial reads B2 X_1U, B4 X_2U, B6 X_3U (alternative contract)

`mut` is a set of mutation labels used only by the controls.
"""
from qq16mat import A, T, U, vadd, vsub, dadd, dsub

TAIL = {0: 2, 1: 4, 2: 6, 3: 8}          # C_k = w_k + K_t f_t  (PR2174 tail identity, all banks)


def sums(o, k):
    pA = {p: o['Xf'][k][p] + o['Yf'][k][p] for p in A}
    pT = {p: o['Xf'][k][p] + o['Yf'][k][p] for p in T}
    return pA, pT


def w_init(M, o, k, mut=frozenset()):
    """Initial Y_kU from the tail identity C_k = w_k + K_t f_t."""
    t = TAIL[k]
    if 'tail_sign' in mut:
        return dadd(o['C'][k], M.Kop(t, o['f'][t]))
    if 'tail_chart' in mut:
        t = t - 1 if t > 2 else t + 1
        return dsub(o['C'][k], M.Kop(t, o['f'][TAIL[k]]))
    return dsub(o['C'][k], M.Kop(t, o['f'][t]))


def decode_bank0(M, o, dA0, mut=frozenset()):
    """PR2174 (5) with p_A, p_T read from the final selected fields (conservation) and the
    initial Y_0U from the tail identity."""
    pA, pT = sums(o, 0)
    f0, f2 = o['f'][0], o['f'][2]
    q0, q2 = o['INIT'][(0, 0)], o['INIT'][(0, 2)]
    w0 = w_init(M, o, 0, mut)
    inner = vadd(M.B(2, w0), q2)
    if 'drop_K0f' not in mut:
        inner = vsub(inner, M.B(2, M.Kop(0, f0)))
    inner = vadd(inner, M.lT(2, dsub(pT, M.PT(0, f0))))
    sec = M.PA(0 if 'section_chart' in mut else 2, inner)
    t = dadd(dsub(pA, dA0), sec)
    xT = M.PT(0, M.lA(0, t))
    xA = M.PA(0, vsub(vsub(f0, q0), M.lT(0, xT)))
    return {'XA': xA, 'XT': xT, 'YA': dsub(pA, xA), 'YT': dsub(pT, xT)}


def V(M, o, k, mut=frozenset()):
    """V_k = w_k - C_(k-1) + C_(k-2) - ... ; then Y_kU at bank k's first B input is V_k + X_kU^initial."""
    v = w_init(M, o, k, mut)
    sign = -1
    for j in range(k - 1, -1, -1):
        if 'drop_alternating' in mut and j < k - 1:
            break
        v = dsub(v, o['C'][j]) if sign < 0 else dadd(v, o['C'][j])
        sign = -sign
    return v


def beta(M, o, k, mut=frozenset()):
    """b + 12a of bank k's first backward call, read from the final X_A (ordered: ell_sT P_tT ell_tA)."""
    s, t = 2 * k + 1, 2 * k + 2
    if 'beta_swap_order' in mut:          # commuting-route-only variant: D_t^-1 D_s instead of D_s D_t^-1
        x = M.lA(t, {p: o['Xf'][k][p] for p in A})
        res = []
        K = x[0].ncols()
        qp = {}
        for r in T:
            acc = M.zero(K)
            for i in range(8):
                c = M.R[i][r]
                if c:
                    acc = acc + x[i] * (c % M.P)
            qp[r] = M.Di[t][r] * (M.D[s][r] * (acc * M.i4))
        for i in range(8):
            acc = M.zero(K)
            for r in T:
                c = M.R[i][r]
                if c:
                    acc = acc + qp[r] * (c % M.P)
            res.append(acc)
        return res
    return M.lT(s, M.PT(t, M.lA(t, {p: o['Xf'][k][p] for p in A})))


def decode_bank(M, o, k, dA=None, mut=frozenset()):
    """Banks k = 1, 2, 3 (first call B_s, then F_t, s = 2k+1, t = s+1).  dA is the retained
    post-B_s Y_A record; for k = 3 it is derived from f8 when dA is None.  For the alternative
    contract pass dA='XB' to derive it from f_t and the cross-bank read XB[k]."""
    s, t = 2 * k + 1, 2 * k + 2
    pA, pT = sums(o, k)
    v = V(M, o, k, mut)
    Bs_w = vadd(M.B(s, v), o['INIT'][(k, s)])
    Bt_w = vadd(M.B(t, v), o['INIT'][(k, t)])
    if 'drop_init' in mut:
        Bs_w = M.B(s, v)
    be = beta(M, o, k, mut)
    if dA is None or dA == 'XB':
        # f_t = B_t(X_kU at F_t) + ell_tA(p_A - u_A) + ell_tT P_sT beta, and ell_tT P_sT beta = ell_tA Xf_kA.
        # X_kU at F_t = X_kU^init - Y_kU^(B_s in) [+ X_(k+1)U^init for k = 1, 2].
        Bt_x = vsub(o['INIT'][(k, t)], Bt_w)
        if dA == 'XB':
            Bt_x = vadd(Bt_x, o['XB'][k])
        known = vadd(Bt_x, M.lA(t, {p: o['Xf'][k][p] for p in A}))
        dA = dadd(pA, M.PA(t, vsub(known, o['f'][t])))
    c = 12 if 'eleven_to_twelve' in mut else 11
    num = vsub(vsub(vsub(be, Bs_w), M.lA(s, dA)), M.lT(s, pT))
    inv = pow(c, M.P - 2, M.P)
    xT = M.PT(s, [n * inv for n in num])
    yT = dsub(pT, xT)
    yA = dsub(dA, M.PA(t, vadd(Bt_w, M.lT(t, yT))))
    xA = dsub(pA, yA)
    return {'XA': xA, 'XT': xT, 'YA': yA, 'YT': yT}


def dA0_from_XB(M, o, mut=frozenset()):
    """Alternative contract for bank 0: d_A0 = p_A + Xf_0A - P_2A(f2 + B2 w0 - B2 X_1U)."""
    pA, _ = sums(o, 0)
    w0 = w_init(M, o, 0, mut)
    inner = vsub(vadd(o['f'][2], M.B(2, w0)), o['XB'][0])
    return dsub(dadd(pA, {p: o['Xf'][0][p] for p in A}), M.PA(2, inner))


def truth(inputs, k):
    x, y = inputs[k]
    return {'XA': {p: x[p] for p in A}, 'XT': {p: x[p] for p in T},
            'YA': {p: y[p] for p in A}, 'YT': {p: y[p] for p in T}}


def equal(dec, tru):
    return all(dec[b][p] == tru[b][p] for b in tru for p in tru[b])
