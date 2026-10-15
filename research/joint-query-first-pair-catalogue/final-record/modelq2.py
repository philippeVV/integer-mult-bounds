#!/usr/bin/env python3
"""Finite noncommuting-model kernel census for retained QQ16 records under the FINAL contract.

Tests PR2346's open question 2 (retain the post-B1 A record d_A and test FULL recovery, not only
PR2346's family), and the cheapest per-bank extensions of it, against the final-record contract:
  FINAL  = final selected X/Y of all four pairs, final Y ordinary fields (completed-X descendants)
           of all pairs, forward histories f0,f2,f4,f6,f8;
  INIT64 = the 64 initial bank-local query coordinates of PR2346 (pair 0: B0,B2 on X_0U; pair k>=1:
           B_(2k+1), B_(2k+2) on X_kU).
Records (each a complete array set, retained to the end):
  dA_k   = Y_A of pair k right after its first backward call (events 1, 5, 9, 13);
  cut_k  = all 32 selected X/Y of pair k at that cut;
  pT_k   = X_T of pair k at that cut, pA_k = X_A.
For each contract we report dim J_k(ker O) / m for each bank k, i.e. the number of lost selected
arrays per owner in the model (selected inputs of bank k; all four banks have arbitrary inputs).

Model (evidence only, not the actual cover): payload/owner space Q^m; chart 0 is the identity and
every (stage>=1, port) route is an independent random permutation of m points; all maps are the
CONSERVED-QQ16 (1)/(3) words with these routes (R_s = R D_s, K_s = D_s^-1 K, P_s, ell_s, E_t with
t = s+1), bridges/terminal as in PR2346 section 1.  Ranks are computed modulo two primes with
python-flint nmod_mat.  The recipe and R are bound to the pinned sources.

usage: modelq2.py M SEED
Prepared by Philippe Van Velzen with substantial Claude Code (Claude Opus 5.5) assistance. Apache-2.0.
"""
import gzip
import json
import random
import sys
import tarfile
from hashlib import sha256
from io import BytesIO
from pathlib import Path

import flint
import numpy as np

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
PRIMES = (2147483647, 2147483629)


def need(ok, msg):
    if not ok:
        raise SystemExit('FAIL ' + msg)


def load():
    raw = {n: (ROOT / n).read_bytes() for n in (PACKET, ARCHIVE)}
    for n in raw:
        need(sha256(raw[n]).hexdigest() == PINS[n], 'pin ' + n)
    pkt = json.loads(gzip.decompress(raw[PACKET]))
    with tarfile.open(fileobj=BytesIO(raw[ARCHIVE]), mode='r:gz') as tar:
        data = tar.extractfile(TEMPLATE).read()
    need(sha256(data).hexdigest() == PINS[TEMPLATE], 'template pin')
    recipe = json.loads(gzip.decompress(data))['recipe']
    R = [[int(row.get(str(p), '0')) for p in range(112)] for row in pkt['R_rows']]
    need((tuple(pkt['A_ports']), tuple(pkt['T_ports']), tuple(pkt['U_ports'])) == (A, T, U), 'partition')
    for ports, g in ((U, 48), (A, 4), (T, 4)):
        need(all(sum(R[i][p] * R[j][p] for p in ports) == g * (i == j) for i in range(8) for j in range(8)), 'Gram')
    need(len(recipe) == 19 and [e['pair'] for e in recipe if e['kind'] == 'local'] == [0, 0, 0, 1, 1, 2, 2, 3, 3],
         'recipe')
    need(all(e['reverse'] == bool(e['stage'] % 2) for e in recipe if e['kind'] == 'local'), 'orientation')
    return R, recipe


class Model:
    def __init__(self, R, m, seed, P, c11=11, c12=12):
        self.R, self.m, self.P = R, m, P
        self.i4 = pow(4, P - 2, P)
        rng = random.Random(seed)
        # D[s][p]: permutation; apply(v) puts v[i] at position D[i]
        self.D = {0: [list(range(m))] * 112}
        for s in range(1, 10):
            self.D[s] = [rng.sample(range(m), m) for _ in range(112)]
        self.inv = {s: [np.argsort(perm) for perm in self.D[s]] for s in self.D}
        self.c11, self.c12 = c11, c12

    # an array value is an (m, ncols) block
    def ap(self, s, p, v):          # D_s(p) v : w[D[i]] = v[i]  ->  w = v[inv]
        return v[self.inv[s][p]]

    def api(self, s, p, v):         # D_s(p)^-1 v : w[i] = v[D[i]]
        return v[self.D[s][p]]

    def query(self, s, x, ports):
        P = self.P
        out = []
        for i in range(8):
            acc = 0
            for p in ports:
                c = self.R[i][p]
                if c:
                    acc = (acc + c * self.ap(s, p, x[p])) % P
            out.append(acc)
        return out

    def sect(self, s, q, ports):    # D_s^-1 R_V^T q / 4 on V
        P = self.P
        res = {}
        for p in ports:
            acc = 0
            for i in range(8):
                c = self.R[i][p]
                if c:
                    acc = (acc + c * q[i]) % P
            res[p] = self.api(s, p, acc * self.i4 % P)
        return res

    def forward(self, s, x, y):
        P = self.P
        f = self.query(s, x, range(112))
        k = self.sect(s, f, U)
        nx = {p: (x[p] - k[p]) % P for p in U}
        nx.update(self.sect(s, self.query(s, x, T), A))
        nx.update(self.sect(s, f, T))
        ny = {p: (y[p] + nx[p]) % P for p in U}
        ny.update({p: (x[p] + y[p] - nx[p]) % P for p in S})
        return nx, ny, f

    def backward(self, s, t, x, y):
        P = self.P
        # u = E_t y = y + P_(t,A)(R_t y - ell_(t,A) y): A rows change only
        qa = [(g + h) % P for g, h in zip(self.query(t, y, U), self.query(t, y, T))]
        corr = self.sect(t, qa, A)
        u = dict(y)
        for p in A:
            u[p] = (y[p] + corr[p]) % P
        a = self.query(s, x, T)
        b = self.query(s, u, range(112))
        ka = self.sect(s, a, U)
        t11 = self.sect(s, [(bi + self.c11 * ai) % P for bi, ai in zip(b, a)], T)
        t12 = self.sect(s, [(bi + self.c12 * ai) % P for bi, ai in zip(b, a)], T)
        nx, ny = {}, {}
        for p in U:
            ny[p], nx[p] = (y[p] + ka[p]) % P, (x[p] - y[p]) % P
        for p in A:
            ny[p], nx[p] = u[p], (x[p] + y[p] - u[p]) % P
        for p in T:
            ny[p], nx[p] = (u[p] - t11[p]) % P, t12[p]
        return nx, ny


def run(model, recipe):
    m, P = model.m, model.P
    n = 896 * m
    eye = np.eye(n, dtype=np.int64)
    # bank k, side s (0 X, 1 Y), port p -> rows of the identity
    st = [[{p: eye[((k * 2 + s) * 112 + p) * m:((k * 2 + s) * 112 + p + 1) * m] for p in range(112)}
           for s in (0, 1)] for k in range(4)]
    init = [dict(st[k][0]) for k in range(4)]
    hist, cuts = {}, {}
    first_b = {}
    for e in recipe:
        if e['kind'] == 'local':
            k, s = e['pair'], e['stage']
            x, y = st[k]
            if e['reverse']:
                nx, ny = model.backward(s, s + 1, x, y)
                if k not in first_b:
                    first_b[k] = s
                    cuts[k] = ({p: nx[p] for p in S}, {p: ny[p] for p in S})
            else:
                nx, ny, f = model.forward(s, x, y)
                hist[s] = f
            st[k] = [nx, ny]
        elif e['kind'] == 'bridge':
            i, j, c = e['src'], e['dst'], e['coefficient']
            for p in U:
                st[j][0][p] = (st[j][0][p] + c * st[i][0][p]) % P
            for p in U:
                st[i][1][p] = (st[i][1][p] - c * st[j][1][p]) % P
        else:
            for k in range(4):
                for p in U:
                    st[k][0][p], st[k][1][p] = st[k][1][p], (-st[k][0][p]) % P
    need(sorted(hist) == [0, 2, 4, 6, 8] and first_b == {0: 1, 1: 3, 2: 5, 3: 7}, 'schedule shape')
    blocks = {}
    blocks['FINAL'] = ([st[k][s][p] for k in range(4) for s in (0, 1) for p in S]
                       + [st[k][1][p] for k in range(4) for p in U]
                       + [q for s in sorted(hist) for q in hist[s]])
    init_stages = {0: (0, 2), 1: (3, 4), 2: (5, 6), 3: (7, 8)}
    blocks['INIT64'] = [q for k in range(4) for s in init_stages[k] for q in model.query(s, init[k], U)]
    for k in range(4):
        cx, cy = cuts[k]
        blocks[f'dA{k}'] = [cy[p] for p in A]
        blocks[f'pA{k}'] = [cx[p] for p in A]
        blocks[f'pT{k}'] = [cx[p] for p in T]
        blocks[f'dT{k}'] = [cy[p] for p in T]
    blocks['INIT0'] = blocks['INIT64'][:16]
    blocks['F02'] = list(hist[0]) + list(hist[2])
    blocks['Y0FIN'] = [st[0][1][p] for p in U]
    blocks['YFIN'] = [st[k][1][p] for k in range(4) for p in U]
    blocks['FALL'] = [q for s in sorted(hist) for q in hist[s]]
    blocks['IN0'] = [eye[i * m:(i + 1) * m] for i in range(224)]     # all pair-0 inputs revealed
    return {kname: np.vstack(v) % P for kname, v in blocks.items()}


def nrank(M, P):
    M = M[np.any(M != 0, axis=1)]
    if M.shape[0] == 0:
        return 0
    return flint.nmod_mat(M.shape[0], M.shape[1], [int(v) for v in M.ravel()], P).rank()


def census(blocks, contracts, m, P):
    n = 896 * m
    out = {}
    Jrows = {}
    for k in range(4):
        idx = [((k * 2 + s) * 112 + p) * m + i for s in (0, 1) for p in S for i in range(m)]
        J = np.zeros((len(idx), n), dtype=np.int64)
        J[np.arange(len(idx)), idx] = 1
        Jrows[k] = J
    for name, parts in contracts:
        O = np.vstack([blocks[b] for b in parts])
        rO = nrank(O, P)
        lost = [Q_(nrank(np.vstack([O, Jrows[k]]), P) - rO, m) for k in range(4)]
        out[name] = {'rank_O': rO, 'lost_selected_arrays_per_bank': lost}
    return out


def Q_(a, m):
    from fractions import Fraction
    v = Fraction(a, m)
    return int(v) if v.denominator == 1 else str(v)


CONTRACTS = [
    ('FINAL', ['FINAL']),
    ('FINAL+INIT64 (PR2346 contract)', ['FINAL', 'INIT64']),
    ('+dA0 (PR2346 Q2 / #2192 record)', ['FINAL', 'INIT64', 'dA0']),
    ('+dA0..dA3 (per-bank d_A, 32)', ['FINAL', 'INIT64', 'dA0', 'dA1', 'dA2', 'dA3']),
    ('+dA0,pA0,pT0 (PR2174 cut of pair 0)', ['FINAL', 'INIT64', 'dA0', 'pA0', 'pT0']),
    ('+per-bank dA,pA,pT (96)', ['FINAL', 'INIT64'] + [f'{x}{k}' for k in range(4) for x in ('dA', 'pA', 'pT')]),
    ('+per-bank full selected cut (128)', ['FINAL', 'INIT64'] + [f'{x}{k}' for k in range(4) for x in ('dA', 'pA', 'pT', 'dT')]),
    ('FINAL+dA0..dA3 (no INIT64)', ['FINAL', 'dA0', 'dA1', 'dA2', 'dA3']),
    ('FINAL+INIT64+all pair-0 inputs revealed', ['FINAL', 'INIT64', 'IN0']),
    ('+dA0,dA1,dA2 (24)', ['FINAL', 'INIT64', 'dA0', 'dA1', 'dA2']),
    ('+dA1,dA2 only', ['FINAL', 'INIT64', 'dA1', 'dA2']),
    ('control: PR2174 cut contract, literal (final Y_0U only)', ['INIT0', 'dA0', 'pA0', 'pT0', 'dT0', 'F02', 'Y0FIN']),
    ('control: PR2174 cut contract + all final Y_U', ['INIT0', 'dA0', 'pA0', 'pT0', 'dT0', 'F02', 'YFIN']),
    ('PR2174 cut contract + all final Y_U + all f', ['INIT0', 'dA0', 'pA0', 'pT0', 'dT0', 'FALL', 'YFIN']),
]


def main(argv):
    m = int(argv[0]) if argv else 3
    seed = int(argv[1]) if len(argv) > 1 else 20261015
    R, recipe = load()
    res = {'m': m, 'seed': seed, 'primes': {}}
    for P in PRIMES:
        model = Model(R, m, seed, P)
        blocks = run(model, recipe)
        out = census(blocks, CONTRACTS, m, P)
        res['primes'][str(P)] = out
        for name, row in out.items():
            print(f'p={P} m={m} seed={seed} {name:40s} rankO={row["rank_O"]:5d} lost/bank={row["lost_selected_arrays_per_bank"]}',
                  flush=True)
    a, b = (res['primes'][str(P)] for P in PRIMES)
    need(a['control: PR2174 cut contract, literal (final Y_0U only)']['lost_selected_arrays_per_bank'][0] == 0, 'PR2174 control')
    need(a['+dA0 (PR2346 Q2 / #2192 record)']['lost_selected_arrays_per_bank'][0] == 0, 'bank-0 repair')
    need(a == b, 'two primes agree')
    Path(__file__).with_name(f'MODELQ2-m{m}-s{seed}.json').write_text(json.dumps(res, indent=1))
    print('PASS_QQ16_MODEL_Q2_CENSUS')


if __name__ == '__main__':
    main(sys.argv[1:])
