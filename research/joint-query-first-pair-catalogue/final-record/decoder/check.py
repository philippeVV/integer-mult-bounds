#!/usr/bin/env python3
"""QQ16 FINAL-contract decoder and bank-k loss witness: checker with controls.

usage: check.py [--quick]      (needs python-flint; REPO=<main checkout>)

A. Ordered decoders (decode.py) recover every selected input of the banks they claim, by randomized
   exact substitution over two primes, with random NONCOMMUTING matrix routes (m = 3, 8), with D_0 = I
   (actual) and D_0 random, plus permutation, diagonal and identity route families.
B. Decoder mutations are rejected; an order-swapped (commuting-only) variant passes with diagonal
   routes and is rejected with matrix routes (the test sees operator order).
C. Simulator mutations (bridge sign, F2/F8 pair reassignment, backward 11 -> 12) break the decoders.
D. Rank census of contracts in a generic matrix model (m = 2) and a permutation model (m = 3).
E. Bank-1 and bank-2 loss witnesses assembled from any z in ker B_s' cap ker B_t' with B_t z != 0,
   verified by full simulation: every FINAL+INIT64+dA0 observation is zero, bank-k selected != 0.
F. Actual-cover rank bound (routes_bound.py).

Prepared by Philippe Van Velzen with substantial Claude Code (Claude Opus 5.5) assistance. Apache-2.0.
"""
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

import flint

import routes_bound
from decode import decode_bank, decode_bank0, dA0_from_XB, equal, truth
from qq16mat import A, INIT_STAGES, PRIMES, S, T, U, Model, load, need, random_inputs, run

MUT0 = ['tail_sign', 'tail_chart', 'drop_K0f', 'section_chart']
MUTK = ['tail_sign', 'tail_chart', 'beta_swap_order', 'drop_init', 'eleven_to_twelve']


def all_decoders(M, o, inp, mut=frozenset()):
    r = {
        'bank0 (dA0)': equal(decode_bank0(M, o, o['dA'][0], mut), truth(inp, 0)),
        'bank1 (dA1)': equal(decode_bank(M, o, 1, o['dA'][1], mut), truth(inp, 1)),
        'bank2 (dA2)': equal(decode_bank(M, o, 2, o['dA'][2], mut), truth(inp, 2)),
        'bank3 (no record, f8)': equal(decode_bank(M, o, 3, None, mut), truth(inp, 3)),
        'bank0 alt (B2 X1U)': equal(decode_bank0(M, o, dA0_from_XB(M, o, mut), mut), truth(inp, 0)),
        'bank1 alt (B4 X2U)': equal(decode_bank(M, o, 1, 'XB', mut), truth(inp, 1)),
        'bank2 alt (B6 X3U)': equal(decode_bank(M, o, 2, 'XB', mut), truth(inp, 2)),
    }
    return r


def noncommuting(M):
    for s in range(1, 9):
        for p in (U[0], A[0], T[0]):
            X, Y = M.D[s][p], M.D[(s % 8) + 1][p + 1]
            if X * Y != Y * X:
                return True
    return False


def part_A(R, recipe, quick):
    fams = [('mat', 8, True), ('mat', 8, False), ('mat', 3, True), ('perm', 4, True), ('diag', 3, True),
            ('ident', 2, True)]
    seeds = (1, 2) if not quick else (1,)
    out = []
    for P in PRIMES:
        for kind, m, d0 in fams:
            for seed in seeds:
                M = Model(R, m, seed, P, kind=kind, d0=not d0)
                if kind == 'mat':
                    need(noncommuting(M), 'matrix routes noncommuting')
                inp = random_inputs(M, 3, random.Random(f'in{seed}{P}{kind}'))
                o = run(M, recipe, inp)
                res = all_decoders(M, o, inp)
                need(all(res.values()), f'decoder {kind} m={m} D0={"I" if d0 else "rand"} p={P}: {res}')
                out.append({'p': P, 'kind': kind, 'm': m, 'D0': 'I' if d0 else 'random', 'seed': seed})
    print(f'A. decoders: {len(out)} runs x 7 decoders, all selected inputs exact (incl. D0 random, '
          'noncommuting m=8 routes, both primes)')
    return out


def part_B(R, recipe):
    res = {}
    for kind in ('mat', 'diag'):
        M = Model(R, 4, 7, PRIMES[0], kind=kind, d0=True)
        inp = random_inputs(M, 2, random.Random(77))
        o = run(M, recipe, inp)
        for mut in MUT0:
            ok = equal(decode_bank0(M, o, o['dA'][0], frozenset([mut])), truth(inp, 0))
            res[f'{kind}/bank0/{mut}'] = ok
            need(not ok, f'mutation {mut} on bank 0 ({kind}) not rejected')
        for mut in MUTK:
            for k in (1, 2):
                ok = equal(decode_bank(M, o, k, o['dA'][k], frozenset([mut])), truth(inp, k))
                res[f'{kind}/bank{k}/{mut}'] = ok
                if mut == 'beta_swap_order' and kind == 'diag':
                    need(ok, 'commuting-only variant should pass with diagonal routes')
                else:
                    need(not ok, f'mutation {mut} on bank {k} ({kind}) not rejected')
        ok = equal(decode_bank(M, o, 2, o['dA'][2], frozenset(['drop_alternating'])), truth(inp, 2))
        need(not ok, 'drop_alternating not rejected')
        res[f'{kind}/bank2/drop_alternating'] = ok
        # wrong retained record: bank 1 decoded with bank 0's d_A
        ok = equal(decode_bank(M, o, 1, o['dA'][0]), truth(inp, 1))
        need(not ok, 'wrong record not rejected')
        res[f'{kind}/bank1/uses dA0'] = ok
    print(f'B. decoder mutations: {sum(1 for v in res.values() if not v)} rejected; '
          'order-swapped beta passes only with commuting (diagonal) routes')
    return res


def part_C(R, recipe):
    res = {}
    rng = random.Random(9)
    cases = {
        'bridge sign flipped': dict(model=dict(), run=dict(bridge_sign=-1)),
        'pair=stage//2 mod 4 (F2->pair1, F8->pair0)': dict(model=dict(), run=dict(pair_map={2: 1, 8: 0})),
        'backward Y_T coefficient 11->12': dict(model=dict(c11=12), run=dict()),
        'backward X_T coefficient 12->13': dict(model=dict(c12=13), run=dict()),
    }
    for name, c in cases.items():
        M = Model(R, 3, 3, PRIMES[0], kind='mat', **c['model'])
        inp = random_inputs(M, 2, rng)
        try:
            o = run(M, recipe, inp, **c['run'])
            r = all_decoders(M, o, inp)
        except (KeyError, ZeroDivisionError) as exc:
            r = {'simulation incompatible': False, 'exc': str(exc)}
        res[name] = r
        need(not all(v for v in r.values() if isinstance(v, bool)), 'simulator mutation not detected: ' + name)
    print(f'C. simulator mutations: {len(res)} detected')
    return res


# ---------------------------------------------------------------- D. census
def identity_inputs(M):
    m, n = M.m, 896 * M.m
    inp = []
    for k in range(4):
        sides = []
        for side in (0, 1):
            d = {}
            for p in range(112):
                base = ((k * 2 + side) * 112 + p) * m
                ent = [0] * (m * n)
                for i in range(m):
                    ent[i * n + base + i] = 1
                d[p] = flint.nmod_mat(m, n, ent, M.P)
            sides.append(d)
        inp.append(sides)
    return inp


def rows(arrs):
    out = []
    for a in arrs:
        e = [int(v) for v in a.entries()]
        n = a.ncols()
        out += [e[i * n:(i + 1) * n] for i in range(a.nrows())]
    return out


def census(R, recipe, kind, m, P):
    M = Model(R, m, 11, P, kind=kind)
    o = run(M, recipe, identity_inputs(M))
    n = 896 * m
    blk = {}
    blk['FINAL'] = rows([o['Xf'][k][p] for k in range(4) for p in S] + [o['Yf'][k][p] for k in range(4) for p in S]
                        + [o['C'][k][p] for k in range(4) for p in U] + [q for s in sorted(o['f']) for q in o['f'][s]])
    blk['INIT64'] = rows([q for key in o['INIT'] for q in o['INIT'][key]])
    for k in range(4):
        blk[f'dA{k}'] = rows([o['dA'][k][p] for p in A])
    for k in range(3):
        blk[f'XB{k}'] = rows(o['XB'][k])
    cx, cy = o['cut'][0]
    blk['PR2174sums'] = rows(o['INIT'][(0, 0)] + o['INIT'][(0, 2)] + [cx[p] + cy[p] for p in S]
                             + [cy[p] for p in A] + o['f'][0] + o['f'][2] + [o['C'][0][p] for p in U])
    blk['PR2174sums-noC'] = rows(o['INIT'][(0, 0)] + o['INIT'][(0, 2)] + [cx[p] + cy[p] for p in S]
                                 + [cy[p] for p in A] + o['f'][0] + o['f'][2])
    J = {}
    for k in range(4):
        idx = [((k * 2 + s) * 112 + p) * m + i for s in (0, 1) for p in S for i in range(m)]
        J[k] = [[int(c == j) for c in range(n)] for j in idx]

    def rank(rr):
        rr = [r for r in rr if any(r)]
        return flint.nmod_mat(len(rr), n, [v for r in rr for v in r], P).rank() if rr else 0

    contracts = [
        ('FINAL+INIT64', ['FINAL', 'INIT64']),
        ('FINAL+INIT64+dA0', ['FINAL', 'INIT64', 'dA0']),
        ('FINAL+INIT64+dA0,dA1,dA2', ['FINAL', 'INIT64', 'dA0', 'dA1', 'dA2']),
        ('FINAL+INIT64+dA1,dA2', ['FINAL', 'INIT64', 'dA1', 'dA2']),
        ('FINAL+INIT64+XB (B2X1U,B4X2U,B6X3U)', ['FINAL', 'INIT64', 'XB0', 'XB1', 'XB2']),
        ('PR2174 contract, sums p_A,p_T + d_A', ['PR2174sums']),
        ('PR2174 sums contract without final C0', ['PR2174sums-noC']),
    ]
    out = {}
    for name, parts in contracts:
        O = [r for b in parts for r in blk[b]]
        rO = rank(O)
        lost = [Fraction(rank(O + J[k]) - rO, m) for k in range(4)]
        out[name] = {'rank_O': rO, 'lost_selected_arrays_per_owner': [str(x) for x in lost]}
    return out


EXPECT = {  # lost selected arrays per owner, generic matrix model / permutation model on m points
    'FINAL+INIT64': lambda g: [g, g, g, 0],
    'FINAL+INIT64+dA0': lambda g: [0, g, g, 0],
    'FINAL+INIT64+dA0,dA1,dA2': lambda g: [0, 0, 0, 0],
    'FINAL+INIT64+dA1,dA2': lambda g: [g, 0, 0, 0],
    'FINAL+INIT64+XB (B2X1U,B4X2U,B6X3U)': lambda g: [0, 0, 0, 0],
}


def part_D(R, recipe):
    res = {}
    for kind, m, gen in (('mat', 2, Fraction(8)), ('perm', 3, Fraction(16, 3))):
        per = {}
        for P in PRIMES:
            per[P] = census(R, recipe, kind, m, P)
        need(per[PRIMES[0]] == per[PRIMES[1]], f'census primes agree ({kind})')
        c = per[PRIMES[0]]
        for name, f in EXPECT.items():
            need([Fraction(x) for x in c[name]['lost_selected_arrays_per_owner']] == f(gen),
                 f'census {kind} {name}: {c[name]}')
        need(c['PR2174 contract, sums p_A,p_T + d_A']['lost_selected_arrays_per_owner'][0] == '0',
             'PR2174 sums contract recovers bank 0')
        need(c['PR2174 sums contract without final C0']['lost_selected_arrays_per_owner'][0] != '0',
             'control: dropping final C0 must lose bank 0')
        res[f'{kind} m={m}'] = c
        for name, row in c.items():
            print(f'D. census {kind} m={m}: {name:42s} rankO={row["rank_O"]:5d} lost/owner={row["lost_selected_arrays_per_owner"]}')
    return res


# ---------------------------------------------------------------- E. witnesses
def witness_inputs(M, k, z):
    """Bank-k (k = 1, 2) loss witness from z in ker B_s' cap ker B_t' (bank k+1's INIT64 stages)."""
    s, t = 2 * k + 1, 2 * k + 2
    zero = M.zero(1)
    inp = [[{p: zero for p in range(112)} for _ in (0, 1)] for _ in range(4)]
    for p in U:
        inp[k + 1][0][p] = z[p] * 11
    v = M.B(t, z)                                    # B_t z
    u = {p: w * 11 for p, w in M.PA(t, v).items()}   # = d_A of bank k
    xT = {p: -w for p, w in M.PT(s, M.lA(s, M.PA(t, v))).items()}
    yT = {p: -xT[p] for p in T}
    corr = M.PA(t, M.lT(t, yT))
    yA = {p: u[p] - corr[p] for p in A}
    for p in A:
        inp[k][0][p], inp[k][1][p] = -yA[p], yA[p]
    for p in T:
        inp[k][0][p], inp[k][1][p] = xT[p], yT[p]
    return inp, u


def qblock(M, s, ports, sign=1):
    """8m x (len(ports) m) matrix of x -> sign * sum_p R[:,p] D_s(p) x_p (row-major ints)."""
    m, P = M.m, M.P
    out = []
    for i in range(8):
        for a in range(m):
            row = []
            for p in ports:
                c = sign * M.R[i][p] % P
                row += [int(M.D[s][p][a, b]) * c % P for b in range(m)]
            out.append(row)
    return out


def bank0_witness(M):
    """PR2346-type bank-0 witness (FINAL+INIT64 blind, d_A0 != 0): unknowns (tt, z0, z1) with
    ell_1A tt = B1 z0, ell_2A tt = -B2 z1, B0 z0 = B2 z0 = 0, B3 z1 = B4 z1 = 0, tt != 0."""
    m, P = M.m, M.P
    nA, nU = 8 * m, 96 * m
    Z = lambda r, c: [[0] * c for _ in range(r)]

    def hcat(*bl):
        return [sum((b[i] for b in bl), []) for i in range(len(bl[0]))]
    rws = (hcat(qblock(M, 1, A), qblock(M, 1, U, -1), Z(nA, nU))
           + hcat(qblock(M, 2, A), Z(nA, nU), qblock(M, 2, U))
           + hcat(Z(nA, nA), qblock(M, 0, U), Z(nA, nU)) + hcat(Z(nA, nA), qblock(M, 2, U), Z(nA, nU))
           + hcat(Z(nA, nA), Z(nA, nU), qblock(M, 3, U)) + hcat(Z(nA, nA), Z(nA, nU), qblock(M, 4, U)))
    N = flint.nmod_mat(len(rws), nA + 2 * nU, [v for r in rws for v in r], P)
    X, nul = N.nullspace()
    for c in range(nul):
        col = [int(X[j, c]) for j in range(nA + 2 * nU)]
        if any(col[:nA]):
            vec = lambda off, j: flint.nmod_mat(m, 1, col[off + j * m: off + (j + 1) * m], P)
            tt = {p: vec(0, j) for j, p in enumerate(A)}
            z0 = {p: vec(nA, j) for j, p in enumerate(U)}
            z1 = {p: vec(nA + nU, j) for j, p in enumerate(U)}
            zero = M.zero(1)
            inp = [[{p: zero for p in range(112)} for _ in (0, 1)] for _ in range(4)]
            xT = M.PT(0, M.lA(0, tt))
            for p in U:
                inp[0][0][p], inp[1][0][p] = z0[p], z1[p]
            for p in A:
                inp[0][0][p], inp[0][1][p] = -tt[p], tt[p]
            for p in T:
                inp[0][0][p], inp[0][1][p] = xT[p], -xT[p]
            return inp, {p: -tt[p] for p in A}, nul
    return None, None, nul


def kernel_z(M, k):
    s2, t2 = INIT_STAGES[k + 1]
    m = M.m
    rws = []
    for st in (s2, t2):
        for i in range(8):
            for a in range(m):
                row = []
                for p in U:
                    Dp = M.D[st][p]
                    row += [int(Dp[a, b]) * (M.R[i][p] % M.P) % M.P for b in range(m)]
                rws.append(row)
    N = flint.nmod_mat(len(rws), 96 * m, [v for r in rws for v in r], M.P)
    X, nul = N.nullspace()
    t = 2 * k + 2
    for c in range(nul):
        z = {p: flint.nmod_mat(m, 1, [int(X[j * m + a, c]) for a in range(m)], M.P) for j, p in enumerate(U)}
        if any(not q == M.zero(1) for q in M.B(t, z)):
            return z, nul
    return None, nul


def part_E(R, recipe):
    res = {}
    for P in PRIMES:
        for kind, m in (('mat', 3), ('perm', 3)):
            M = Model(R, m, 5, P, kind=kind)
            for k in (0, 1, 2):
                if k == 0:
                    inp, u, nul = bank0_witness(M)
                    need(inp is not None, 'bank-0 kernel vector')
                else:
                    z, nul = kernel_z(M, k)
                    need(z is not None, 'kernel vector with B_t z != 0')
                    inp, u = witness_inputs(M, k, z)
                o = run(M, recipe, inp)
                Z = M.zero(1)
                zero_obs = (all(o['Xf'][j][p] == Z and o['Yf'][j][p] == Z for j in range(4) for p in S)
                            and all(o['C'][j][p] == Z for j in range(4) for p in U)
                            and all(q == Z for s in o['f'] for q in o['f'][s])
                            and all(q == Z for key in o['INIT'] for q in o['INIT'][key])
                            and (k == 0 or all(o['dA'][0][p] == Z for p in A)))
                sel = sum(1 for side in (0, 1) for p in S if not inp[k][side][p] == Z)
                dAk = any(not o['dA'][k][p] == Z for p in A)
                need(zero_obs and sel > 0, f'witness bank {k} {kind} p={P}')
                need(dAk and all(o['dA'][k][p] == u[p] for p in A), 'witness is seen by its own d_A record')
                res[f'{kind} m={m} p={P} bank{k}'] = {'nonzero_selected_input_arrays': sel,
                                                      'kernel_dim': nul, 'dA_k_nonzero': dAk}
    print(f'E. witnesses: {len(res)} assembled; banks 1,2: all FINAL+INIT64+dA0 observations zero; bank 0: all '
          'FINAL+INIT64 observations zero; bank-k selected inputs nonzero; each detected by its own d_A_k')
    return res


def main(argv):
    quick = '--quick' in argv
    R, recipe, _ = load()
    result = {'A_decoders': part_A(R, recipe, quick), 'B_decoder_mutations': part_B(R, recipe),
              'C_simulator_mutations': part_C(R, recipe)}
    if not quick:
        result['D_census'] = part_D(R, recipe)
    result['E_witnesses'] = part_E(R, recipe)
    print('F. actual cover:')
    result['F_actual_route_bound'] = routes_bound.main([])
    Path(__file__).with_name('CHECK.json' if not quick else 'CHECK-quick.json').write_text(
        json.dumps(result, indent=1, default=str))
    print('PASS_QQ16_FINAL_CONTRACT_DECODER')


if __name__ == '__main__':
    main(sys.argv[1:])
