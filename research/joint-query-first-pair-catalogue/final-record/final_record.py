#!/usr/bin/env python3
"""QQ16 final-record (d18) retention: exact constant-subspace kernel census.

Question (PR2470 open question 1; PR2346 open question 2).  For the pinned
19-event QQ16 recipe with the conserved F/B maps of CONSERVED-QQ16 (1)/(3),
which retained records must be added to the FINAL observation contract of
PR2346 (final selected fields of all pairs, final Y ordinary fields =
completed-X descendants, the five forward histories f0,f2,f4,f6,f8, and
optionally the 64 initial bank-local query coordinates) so that the
original selected inputs are determined?  In particular: does retaining the
post-B1 pair-0 A record d_A = Y_A (PR2346 Q2 / Discussion 2192) restore
full selected recovery?

Model.  Exact rational linear forms in all 896 initial port arrays on the
owner- and payload-constant subspace V_c, on which every actual cover route
D_s acts as the identity (inherited premise, as in PR2470).  Consequences:
  * every kernel vector found here is an exact input of the actual routed
    maps whose permitted observations all vanish (a counterexample);
  * dim J(ker O|V_c) is a lower bound on the number q of owner-covariant
    rational-linear retained complete arrays (each restricts to one
    functional on V_c) that can repair the contract.  It is NOT a bound
    for arbitrary non-covariant records, and constant-subspace sufficiency
    is NOT an actual-cover sufficiency proof.

Witnesses are re-verified by an independent direct Fraction simulation of
all 19 events (vector code path, no matrices).

usage: final_record.py [--controls]
Needs python-flint.  Prepared by Philippe Van Velzen with substantial Claude
Code (Claude Opus 5.5) assistance.  Apache-2.0.  No predecessor executable is
imported; the map formulas are transcribed from CONSERVED-QQ16.md (1)/(3)
and the bridge/terminal semantics from PR2346 section 1.
"""
import gzip
import json
import sys
import tarfile
from fractions import Fraction as Q
from hashlib import sha256
from io import BytesIO
from pathlib import Path

import flint

ROOT = Path(__import__('os').environ.get('REPO', Path(__file__).resolve().parents[3]))
PACKET = 'research/joint-query/certificate.json.gz'
ARCHIVE = 'research/community-frontier1042-audit/selected-replay.tar.gz'
TEMPLATE = 't16/global/base/split-centres/GLOBAL-TEMPLATE.json.gz'
PINS = {PACKET: 'a5ed85fac5ebf3c091482a113ab0d351eb5a7892fe2841a35e16017ec1ff49f7',
        ARCHIVE: '22c3f4dce4b897a1d423d1388801933e288b6a6b1906cb15e1b16a5eb68db754',
        TEMPLATE: 'a179dad87651aba340a5a8ddbcd86a8e757cad87b8a1f8bad48106fbecd38518'}
A = (0, 7, 20, 27, 36, 46, 57, 111)
T = (21, 37, 47, 76, 81, 87, 92, 110)
S = A + T
U = tuple(p for p in range(112) if p not in S)
NP, N = 4, 896


def need(ok, msg):
    if not ok:
        raise SystemExit('FAIL ' + msg)


def fid(k, side, p):
    return k * 224 + side * 112 + p


def load():
    raw = {}
    for name in (PACKET, ARCHIVE):
        raw[name] = (ROOT / name).read_bytes()
        need(sha256(raw[name]).hexdigest() == PINS[name], 'pin ' + name)
    pkt = json.loads(gzip.decompress(raw[PACKET]))
    with tarfile.open(fileobj=BytesIO(raw[ARCHIVE]), mode='r:gz') as tar:
        data = tar.extractfile(TEMPLATE).read()
    need(sha256(data).hexdigest() == PINS[TEMPLATE], 'template pin')
    recipe = json.loads(gzip.decompress(data))['recipe']
    R = [[Q(row.get(str(p), '0')) for p in range(112)] for row in pkt['R_rows']]
    need((tuple(pkt['A_ports']), tuple(pkt['T_ports']), tuple(pkt['U_ports'])) == (A, T, U), 'partition')
    for ports, g in ((U, 48), (A, 4), (T, 4)):
        need(all(sum(R[i][p] * R[j][p] for p in ports) == g * (i == j) for i in range(8) for j in range(8)),
             'Gram')
    need(pkt['F_contract'].startswith('f=R_s x. X_U=x_U-K_s,U f; X_A=P_s,A ell_s,T x; X_T=P_s,T f.'), 'F contract')
    need(pkt['B_contract'].startswith('u=E_t y; a=ell_s,T x; b=R_s u. Y_A=u_A; Y_U=y_U+K_s,U a; '
                                      'Y_T=u_T-P_s,T(b+11a). X_U=x_U-y_U; X_A=x_A+y_A-u_A; X_T=P_s,T(b+12a).'),
         'B contract')
    need(len(recipe) == 19 and recipe[-1] == {'kind': 'terminal', 'pairs': 4}, 'recipe')
    return R, recipe


# ------------------------------------------------------------- vector maps (Fractions)
class Maps:
    def __init__(self, R, c11=11, c12=12):
        self.R, self.c11, self.c12 = R, c11, c12

    def q(self, v, ports):
        return [sum((self.R[i][p] * v[p] for p in ports), Q(0)) for i in range(8)]

    def sec(self, qv, ports):          # P_V q = R_V^T q / 4 (and K_V = R^T/4 on V)
        return {p: sum((self.R[i][p] * qv[i] for i in range(8)), Q(0)) / 4 for p in ports}

    def forward(self, x, y):
        f = self.q(x, range(112))
        k = self.sec(f, U)
        nx = {p: x[p] - k[p] for p in U}
        nx.update(self.sec(self.q(x, T), A))
        nx.update(self.sec(f, T))
        ny = {p: y[p] + nx[p] for p in U}
        ny.update({p: x[p] + y[p] - nx[p] for p in S})
        return nx, ny, f

    def backward(self, x, y):
        corr = self.sec([g + h for g, h in zip(self.q(y, U), self.q(y, T))], A)
        u = dict(y)
        for p in A:
            u[p] = y[p] + corr[p]
        a = self.q(x, T)
        b = self.q(u, range(112))
        ka = self.sec(a, U)
        t11 = self.sec([bi + self.c11 * ai for bi, ai in zip(b, a)], T)
        t12 = self.sec([bi + self.c12 * ai for bi, ai in zip(b, a)], T)
        nx, ny = {}, {}
        for p in U:
            ny[p], nx[p] = y[p] + ka[p], x[p] - y[p]
        for p in A:
            ny[p], nx[p] = u[p], x[p] + y[p] - u[p]
        for p in T:
            ny[p], nx[p] = u[p] - t11[p], t12[p]
        return nx, ny


def simulate(maps, recipe, v):
    """Independent direct simulation of all 19 events on one rational input vector.
    Returns per-event states (lists of 896 Fractions) and forward histories {event: f}."""
    st = list(v)
    states, hist = [], {}
    for ev, e in enumerate(recipe):
        new = list(st)
        if e['kind'] == 'local':
            k = e['pair']
            x = {p: st[fid(k, 0, p)] for p in range(112)}
            y = {p: st[fid(k, 1, p)] for p in range(112)}
            if e['reverse']:
                nx, ny = maps.backward(x, y)
            else:
                nx, ny, f = maps.forward(x, y)
                hist[ev] = f
            for p in range(112):
                new[fid(k, 0, p)], new[fid(k, 1, p)] = nx[p], ny[p]
        elif e['kind'] == 'bridge':
            i, j, c = e['src'], e['dst'], e['coefficient']
            for p in U:
                new[fid(j, 0, p)] = st[fid(j, 0, p)] + c * st[fid(i, 0, p)]
                new[fid(i, 1, p)] = st[fid(i, 1, p)] - c * st[fid(j, 1, p)]
        else:
            for k in range(NP):
                for p in U:
                    new[fid(k, 0, p)], new[fid(k, 1, p)] = st[fid(k, 1, p)], -st[fid(k, 0, p)]
        st = new
        states.append(st)
    return states, hist


# ------------------------------------------------------------- matrix trajectory (flint)
def local_matrix(maps, reverse):
    cols = []
    for j in range(224):
        x = {p: Q(int(j == p)) for p in range(112)}
        y = {p: Q(int(j == 112 + p)) for p in range(112)}
        nx, ny = maps.backward(x, y) if reverse else maps.forward(x, y)[:2]
        cols.append([nx[p] for p in range(112)] + [ny[p] for p in range(112)])
    return [[cols[j][i] for j in range(224)] for i in range(224)]


def to_fmpq(rows):
    return flint.fmpq_mat(len(rows), len(rows[0]),
                          [flint.fmpq(v.numerator, v.denominator) for r in rows for v in r])


def trajectory(maps, recipe):
    """state_t (896 x 896): row i = linear form of field i after event t in the initial arrays."""
    Lf, Lb = local_matrix(maps, False), local_matrix(maps, True)
    need(to_fmpq(Lf).rank() == 224 and to_fmpq(Lb).rank() == 224, 'local maps invertible')
    Rq = to_fmpq(maps.R)
    states, hist = [], {}
    cur = flint.fmpq_mat(N, N, [flint.fmpq(int(i == j)) for i in range(N) for j in range(N)])
    for ev, e in enumerate(recipe):
        E = [[0] * N for _ in range(N)]
        for i in range(N):
            E[i][i] = 1
        if e['kind'] == 'local':
            k, L = e['pair'], (Lb if e['reverse'] else Lf)
            for a in range(224):
                for b in range(224):
                    E[k * 224 + a][k * 224 + b] = L[a][b]
            if not e['reverse']:   # forward history f = R x of this pair before the call
                X = flint.fmpq_mat(112, N, [cur[fid(k, 0, p), c] for p in range(112) for c in range(N)])
                hist[ev] = Rq * X
        elif e['kind'] == 'bridge':
            i, j, c = e['src'], e['dst'], e['coefficient']
            for p in U:
                E[fid(j, 0, p)][fid(i, 0, p)] = c
                E[fid(i, 1, p)][fid(j, 1, p)] = -c
        else:
            for k in range(NP):
                for p in U:
                    E[fid(k, 0, p)][fid(k, 0, p)] = 0
                    E[fid(k, 0, p)][fid(k, 1, p)] = 1
                    E[fid(k, 1, p)][fid(k, 1, p)] = 0
                    E[fid(k, 1, p)][fid(k, 0, p)] = -1
        Em = flint.fmpq_mat(N, N, [flint.fmpq(v.numerator, v.denominator) if isinstance(v, Q) else flint.fmpq(v)
                                   for r in E for v in r])
        cur = Em * cur
        states.append(cur)
    return states, hist


def rows_of(M, idx):
    return [[M[i, c] for c in range(N)] for i in idx]


def mat_rows(rowlists):
    flat = [v for r in rowlists for v in r]
    return flint.fmpq_mat(len(rowlists), N, flat)


def defect(Orows, Jidx):
    """dim J(ker O) = rank[O;J] - rank O on V_c."""
    O = mat_rows(Orows)
    Jr = [[flint.fmpq(int(c == j)) for c in range(N)] for j in Jidx]
    return mat_rows(Orows + Jr).rank() - O.rank(), O.rank()


def kernel_witness(Orows, Jidx):
    """A rational x with O x = 0 and J x != 0 (or None)."""
    O = mat_rows(Orows)
    num, den = O.numer_denom()
    Z = flint.fmpz_mat(num)
    K, nul = Z.nullspace()
    for c in range(nul):
        v = [Q(int(K[r, c])) for r in range(N)]
        if any(v[j] for j in Jidx):
            return v
    return None


def main(argv):
    R, recipe = load()
    maps = Maps(R)
    states, hist = trajectory(maps, recipe)
    fin = states[18]
    fwd = sorted(hist)
    need(fwd == [0, 3, 7, 11, 14], 'forward events')
    first_b = {0: 1, 1: 5, 2: 9, 3: 13}           # each bank's first backward call
    ident = lambda idx: [[flint.fmpq(int(c == i)) for c in range(N)] for i in idx]
    Rq = to_fmpq(R)

    def rq(idx_ports, k):                           # R_U x_kU at the initial cut (constant chart)
        X = flint.fmpq_mat(112, N, [flint.fmpq(int(c == fid(k, 0, p))) if p in idx_ports else flint.fmpq(0)
                                    for p in range(112) for c in range(N)])
        M = Rq * X
        return [[M[i, c] for c in range(N)] for i in range(8)]

    FINAL = (rows_of(fin, [fid(k, s, p) for k in range(NP) for s in (0, 1) for p in S])   # 128 selected
             + rows_of(fin, [fid(k, 1, p) for k in range(NP) for p in U])                # 384 completed-X
             + [r for ev in fwd for r in [[hist[ev][i, c] for c in range(N)] for i in range(8)]])  # 40
    INIT64 = [r for k in range(NP) for _ in range(2) for r in rq(U, k)]   # B_s X_kU (collapsed charts)
    dA = {k: rows_of(states[first_b[k]], [fid(k, 1, p) for p in A]) for k in range(NP)}
    cut0 = rows_of(states[1], [fid(0, s, p) for s in (0, 1) for p in S])     # all 32 post-B1 selected

    J0S = [fid(0, s, p) for s in (0, 1) for p in S]
    JS = [fid(k, s, p) for k in range(NP) for s in (0, 1) for p in S]
    J0 = list(range(224))
    JALL = list(range(N))
    Js = {'pair0 selected (32)': J0S, 'all selected (128)': JS, 'pair0 all (224)': J0, 'all inputs (896)': JALL}
    cases = [
        ('FINAL', FINAL),
        ('FINAL+INIT64', FINAL + INIT64),
        ('FINAL+INIT64+dA0 (PR2346 Q2)', FINAL + INIT64 + dA[0]),
        ('FINAL+INIT64+dA0..3', FINAL + INIT64 + [r for k in range(NP) for r in dA[k]]),
        ('FINAL+INIT64+post-B1 pair-0 selected cut (32)', FINAL + INIT64 + cut0),
    ]
    result = {'model': 'owner/payload-constant subspace, 896 arrays', 'cases': {}}
    for name, O in cases:
        row = {}
        for jn, J in Js.items():
            d, r = defect(O, J)
            row[jn] = d
        row['rank_O'] = r
        result['cases'][name] = row
        print(f'{name:48s} rankO={r:4d} ' + '  '.join(f'{jn}: {row[jn]}' for jn in Js))

    # positive control: PR2174's recorded-cut contract recovers pair-0 selected inputs
    pc = (INIT64[:16] + rows_of(states[1], [fid(0, 0, p) for p in S] + [fid(0, 1, p) for p in A])
          + [r for ev in (0, 3) for r in [[hist[ev][i, c] for c in range(N)] for i in range(8)]]
          + rows_of(fin, [fid(0, 1, p) for p in U]))
    d_pc = defect(pc, J0S)[0]
    print('control PR2174 recorded-cut contract, pair0 selected defect =', d_pc)
    need(d_pc == 0, 'PR2174 positive control')
    result['control_PR2174_defect'] = d_pc

    # explicit witness for PR2346 Q2, verified by independent direct simulation
    O = FINAL + INIT64 + dA[0]
    w = kernel_witness(O, J0S)
    if w is not None:
        st, hs = simulate(maps, recipe, w)
        f = st[18]
        obs_zero = (all(f[fid(k, s, p)] == 0 for k in range(NP) for s in (0, 1) for p in S)
                    and all(f[fid(k, 1, p)] == 0 for k in range(NP) for p in U)
                    and all(v == 0 for ev in hs for v in hs[ev])
                    and all(v == 0 for k in range(NP) for v in maps.q({p: w[fid(k, 0, p)] for p in range(112)}, U))
                    and all(st[1][fid(0, 1, p)] == 0 for p in A))
        sel = {fid(0, s, p): w[fid(0, s, p)] for s in (0, 1) for p in S if w[fid(0, s, p)]}
        need(obs_zero and sel, 'independent simulation of the Q2 witness')
        supp = sum(1 for v in w if v)
        print(f'Q2 witness: support {supp} arrays, all permitted observations zero (direct simulation), '
              f'{len(sel)} nonzero pair-0 selected inputs')
        result['Q2_witness'] = {'support': supp, 'nonzero_pair0_selected': len(sel),
                                'vector': {str(i): str(v) for i, v in enumerate(w) if v}}
    if '--controls' in argv:
        try:
            trajectory(Maps(R, c11=12), recipe)
            raise SystemExit('FAIL control: Y_T coefficient 11->12 not rejected')
        except SystemExit as exc:
            need(str(exc) == 'FAIL local maps invertible', 'control outcome ' + str(exc))
            print('control 11->12 in the backward Y_T row: rejected (local map no longer invertible)')
        cut_only = FINAL + dA[0]
        d_c = defect(cut_only, J0S)[0]
        print('scope note: FINAL+dA0 without INIT64 already gives pair0 selected defect', d_c,
              'on constants (charts coincide; every finite-model loss below is chart-dependent)')
        need(d_c == 0 and result['cases']['FINAL']['pair0 selected (32)'] == 8, 'constant-subspace scope facts')
    Path(__file__).with_name('FINAL-RECORD.json').write_text(json.dumps(result, indent=1))
    print('PASS_QQ16_FINAL_RECORD_CENSUS')


if __name__ == '__main__':
    main(sys.argv[1:])
