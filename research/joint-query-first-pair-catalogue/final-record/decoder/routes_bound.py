"""Actual-cover rank bound behind the bank-k losses (route data from the pinned certificate).

Claim checked here.  Let D_v(p) = M(g_v(p), delta_v(p)) be the actual chart-v route of port p, read as
the first 112 cover maps (X alignment, ports L0..L111) of the certificate's F_v program, and
B_v = sum_p R[:,p] D_v(p) iota_U.  For a triple (r; s, t) put the 24 x 24 Gram
G = [B_r;B_s;B_t][B_r;B_s;B_t]^*.  Diagonal blocks are 48 I (Gram identity).  The off-diagonal block
B_a B_b^* = sum_p R_p R_p^T D_a(p) D_b(p)^-1, and its identity-coefficient l2 mass is
sum over equal relative elements.  When the 96 relative elements D_a(p) D_b(p)^-1 (p in U) are pairwise
distinct for each ordered pair, tr(G^2) = 3*8*48^2 + 6*96*16 and tr(G) = 3*8*48, so on every free orbit
of the (finite) cover action Cauchy-Schwarz gives
     rank [B_r;B_s;B_t] >= tr(G)^2 / tr(G^2) = 1152^2/64512 = 20.571... per owner,
while rank [B_s;B_t] <= 16.  Hence dim B_r(ker B_s cap ker B_t) >= 4.571 arrays per owner (> 0).

Relevant triples: (4;5,6) bank-1 loss, (6;7,8) bank-2 loss, (2;3,4) bank-0 loss without d_A0.
Control: (2;1,...) relative distinctness for stages 1,2 reproduces the 96 route distinctions
of QQ16-SELECTED-WITNESS; a substituted D_5 := D_4 drives the bound to <= 16 (rejected).

The owner parts are 73 x 73 binary matrices (columns as integers); composition follows the certificate:
(g,a)(h,b) = (gh, a XOR g b), inverse (g^-1, g^-1 a).  All actual routes are involutions (checked).
Distinctness is checked on owner parts AND on full (g, delta) elements.
"""
import gzip
import json
import sys
from fractions import Fraction
from itertools import permutations

from qq16mat import PACKET, ROOT, U, load, need

ID = tuple(1 << i for i in range(73))


def app(g, c):
    v, i = 0, 0
    while c:
        if c & 1:
            v ^= g[i]
        c >>= 1
        i += 1
    return v


def mul(g, h):
    return tuple(app(g, c) for c in h)


def compose(x, y):
    (g, a), (h, b) = x, y
    return (mul(g, h), a ^ app(g, b))


def inverse_involution(x):
    g, a = x
    need(mul(g, g) == ID, 'owner involution')
    return (g, app(g, a))


def routes(pkt):
    mr = pkt['map_registry']
    out = {}
    for v in range(1, 9):
        cm = [e for e in pkt['programs'][f'F_{v}'] if e[0] == 'cover_map'][:112]
        need([e[1] for e in cm] == [['L', p] for p in range(112)], f'F_{v} X alignment order')
        out[v] = [(tuple(mr[e[2]]['owner_right_matrix_columns']), mr[e[2]]['payload_translation_at_identity_owner'])
                  for e in cm]
    return out


def y_alignment_owner_check(pkt, D):
    """Owner parts of the Y alignments to the future chart in each B_s_t program equal rho_t(p)."""
    mr = pkt['map_registry']
    checked = 0
    for (s, t) in ((1, 2), (3, 4), (5, 6), (7, 8)):
        cm = [e for e in pkt['programs'][f'B_{s}_{t}'] if e[0] == 'cover_map' and 112 <= e[1][1] < 224]
        owners_t = {D[t][p][0] for p in range(112)}
        owners_s = {D[s][p][0] for p in range(112)}
        for e in cm:
            p = e[1][1] - 112
            g = tuple(mr[e[2]]['owner_right_matrix_columns'])
            need(g in (D[t][p][0], D[s][p][0]), f'B_{s}_{t} Y alignment owner of port {p}')
            checked += 1
        need(owners_t and owners_s, 'owners')
    return checked


def gram_trace_bound(R, D, triple):
    """Return (tr G, tr G^2, bound) with exact accounting of coinciding relative elements."""
    trG = Fraction(3 * 8 * 48)
    # diagonal blocks: (48 I)^2 -> 8 * 48^2 each (Gram identity, checked in load())
    trG2 = Fraction(3 * 8 * 48 * 48)
    detail = {}
    for a, b in permutations(triple, 2):
        groups = {}
        for p in U:
            rel = compose(D[a][p], inverse_involution(D[b][p]))
            groups.setdefault(rel, []).append(p)
        mass = 0
        for ports in groups.values():
            Mx = [[sum(R[i][p] * R[j][p] for p in ports) for j in range(8)] for i in range(8)]
            mass += sum(v * v for row in Mx for v in row)
        owner_distinct = len({compose(D[a][p], inverse_involution(D[b][p]))[0] for p in U})
        detail[f'{a},{b}'] = {'distinct_full': len(groups), 'distinct_owner': owner_distinct, 'l2_mass': mass}
        trG2 += mass
    return trG, trG2, trG * trG / trG2, detail


def main(argv):
    R, _, pkt = load()
    need(all(sum(R[i][p] ** 2 for i in range(8)) == 4 for p in range(112)), 'port column norms 4')
    D = routes(pkt)
    for v in D:
        for x in D[v]:
            inverse_involution(x)
            need(app(x[0], x[1]) == x[1], 'translation fixed by owner (involutive element)')
    ny = y_alignment_owner_check(pkt, D)
    print(f'actual routes: 8 charts x 112 ports, all involutions; {ny} Y-alignment owners match rho_s/rho_t')
    res = {}
    # control: stage-1/2 relative routes (the QQ16-SELECTED-WITNESS check)
    rel12 = {compose(D[1][p], inverse_involution(D[2][p]))[0] for p in U}
    need(len(rel12) == 96, 'control: 96 distinct stage-1/2 relative owner routes')
    print('control: 96 distinct stage-1/2 relative owner routes (QQ16-SELECTED-WITNESS premise) reproduced')
    for name, triple in (('bank1 (4;5,6)', (4, 5, 6)), ('bank2 (6;7,8)', (6, 7, 8)), ('bank0 (2;3,4)', (2, 3, 4))):
        trG, trG2, bound, detail = gram_trace_bound(R, D, triple)
        need(all(d['distinct_owner'] == 96 for d in detail.values()), name + ' relative owner distinctness')
        need(bound > 16, name + ' bound')
        res[name] = {'trG': str(trG), 'trG2': str(trG2), 'rank_lower_bound': str(bound),
                     'loss_lower_bound_per_owner': str(bound - 16), 'pairs': detail}
        print(f'{name}: tr G={trG} tr G^2={trG2} rank[B_r;B_s;B_t] >= {float(bound):.4f} > 16; '
              f'dim B_r(ker B_s cap ker B_t) >= {float(bound - 16):.4f} per owner')
    # bank 0 under FINAL+INIT64 (PR2346's theorem, reproduced): needs tt in P_1A B1(ker B0 cap ker B2)
    # intersected with P_2A B2(ker B3 cap ker B4); D_0 = I.
    D0 = dict(D)
    D0[0] = [(ID, 0)] * 112
    trG, trG2, b102, detail = gram_trace_bound(R, D0, (1, 0, 2))
    need(all(d['distinct_owner'] == 96 for d in detail.values()), 'bank0 (1;0,2) relative owner distinctness')
    b0 = res['bank0 (2;3,4)']
    inter = (b102 - 16) + Fraction(b0['loss_lower_bound_per_owner']) - 8
    need(inter > 0, 'bank-0 intersection bound')
    res['bank0 under FINAL+INIT64 (PR2346 reproduced)'] = {'(1;0,2) bound': str(b102 - 16),
                                                          'intersection_lower_bound_per_owner': str(inter)}
    print(f'bank0 (1;0,2): image >= {float(b102 - 16):.4f}; intersection with (2;3,4) image >= {inter} = '
          f'{float(inter):.4f} per owner (PR2346 ceil(8|G|/7) reproduced)')
    # mutation control: identify chart 5 with chart 4 -> bound must collapse
    Dm = dict(D)
    Dm[5] = D[4]
    _, _, bm, _ = gram_trace_bound(R, Dm, (4, 5, 6))
    need(bm <= 16, 'control: D5 := D4 must not certify')
    print(f'control D5:=D4: bound {float(bm):.4f} <= 16 (no certificate, as required)')
    res['control_D5_eq_D4_bound'] = str(bm)
    if '--json' in argv:
        print(json.dumps(res, indent=1))
    return res


if __name__ == '__main__':
    main(sys.argv[1:])
    print('PASS_QQ16_ACTUAL_ROUTE_RANK_BOUND')
