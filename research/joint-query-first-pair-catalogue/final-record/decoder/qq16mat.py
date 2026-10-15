"""Noncommutative operator model of the pinned 19-event QQ16 recipe (library).

Every (chart s, port p) route D_s(p) is an operator on the payload space F_P^m.  The default is a
uniformly random invertible m x m matrix (noncommuting), so the CONSERVED-QQ16 (1)/(3) words are
evaluated with genuinely ordered products.  Optional route families: 'perm' (random permutation
matrices, the modelq2 family), 'diag' (random invertible diagonal matrices: commuting control), and
'ident' (D = I: the owner/payload-constant subspace).  Chart 0 is the identity unless d0=True.

An array value is an m x K nmod_mat (K simultaneous input columns).  The R coefficients are exact
integers; 1/4 and 1/11 are modular inverses (both primes are odd and != 11).

Prepared by Philippe Van Velzen with substantial Claude Code (Claude Opus 5.5) assistance. Apache-2.0.
The map formulas are transcribed from CONSERVED-QQ16.md (1)/(3) and the bridge/terminal semantics from
PR2346 section 1 / PR2174 section 2; no predecessor executable is imported.
"""
import gzip
import json
import os
import random
import tarfile
from hashlib import sha256
from io import BytesIO
from pathlib import Path

import flint

ROOT = Path(os.environ.get('REPO', str(Path(__file__).resolve().parents[4])))
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
PRIMES = (2147483647, 4611686018427387847)      # 2^31-1 and a 62-bit prime
FIRST_B = {0: (1, 2), 1: (3, 4), 2: (5, 6), 3: (7, 8)}   # (s, t) of each bank's backward call
INIT_STAGES = {0: (0, 2), 1: (3, 4), 2: (5, 6), 3: (7, 8)}  # PR2346 INIT64 reads B_s X_kU


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
        need(all(sum(R[i][p] * R[j][p] for p in ports) == g * (i == j) for i in range(8) for j in range(8)),
             'Gram')
    need(len(recipe) == 19 and recipe[-1] == {'kind': 'terminal', 'pairs': 4}, 'recipe length')
    need([e['pair'] for e in recipe if e['kind'] == 'local'] == [0, 0, 0, 1, 1, 2, 2, 3, 3], 'pair sequence')
    need(all(e['reverse'] == bool(e['stage'] % 2) for e in recipe if e['kind'] == 'local'), 'orientation')
    return R, recipe, pkt


class Model:
    """Operator model.  c11/c12 are the backward Y_T / X_T coefficients (mutation hooks)."""

    def __init__(self, R, m, seed, P, kind='mat', d0=False, c11=11, c12=12):
        self.R, self.m, self.P, self.kind = R, m, P, kind
        self.c11, self.c12 = c11, c12
        self.i4 = pow(4, P - 2, P)
        self.i11 = pow(11, P - 2, P)
        rng = random.Random(f'{seed}:{m}:{kind}:{P}')
        I = flint.nmod_mat(m, m, [int(i == j) for i in range(m) for j in range(m)], P)
        self.I = I
        self.D, self.Di = {}, {}
        for s in range(9):
            self.D[s], self.Di[s] = [], []
            for p in range(112):
                M = I if (s == 0 and not d0) or kind == 'ident' else self._rand(rng, kind)
                self.D[s].append(M)
                self.Di[s].append(M.inv())

    def _rand(self, rng, kind):
        m, P = self.m, self.P
        if kind == 'perm':
            pi = rng.sample(range(m), m)
            return flint.nmod_mat(m, m, [int(pi[j] == i) for i in range(m) for j in range(m)], P)
        if kind == 'diag':
            return flint.nmod_mat(m, m, [rng.randrange(1, P) if i == j else 0 for i in range(m) for j in range(m)], P)
        while True:
            M = flint.nmod_mat(m, m, [rng.randrange(P) for _ in range(m * m)], P)
            if M.rank() == m:
                return M

    def zero(self, K):
        return flint.nmod_mat(self.m, K, self.P)

    # ---- ordered primitives (chart s) ------------------------------------------------------
    def query(self, s, x, ports):
        """sum_p R[:,p] D_s(p) x_p  -> list of 8 arrays."""
        K = next(iter(x.values())).ncols()
        Dx = {p: self.D[s][p] * x[p] for p in ports}
        out = []
        for i in range(8):
            acc = self.zero(K)
            for p in ports:
                c = self.R[i][p]
                if c:
                    acc = acc + Dx[p] * (c % self.P)
            out.append(acc)
        return out

    def sect(self, s, q, ports):
        """p -> D_s(p)^-1 (R[:,p] . q)/4 for p in ports (K_s on U, P_sA on A, P_sT on T)."""
        K = q[0].ncols()
        res = {}
        for p in ports:
            acc = self.zero(K)
            for i in range(8):
                c = self.R[i][p]
                if c:
                    acc = acc + q[i] * (c % self.P)
            res[p] = self.Di[s][p] * (acc * self.i4)
        return res

    # named operators
    def B(self, s, z):
        return self.query(s, z, U)

    def Kop(self, s, q):
        return self.sect(s, q, U)

    def PA(self, s, q):
        return self.sect(s, q, A)

    def PT(self, s, q):
        return self.sect(s, q, T)

    def lA(self, s, x):
        return self.query(s, x, A)

    def lT(self, s, x):
        return self.query(s, x, T)

    # ---- CONSERVED-QQ16 words ---------------------------------------------------------------
    def forward(self, s, x, y):
        f = self.query(s, x, range(112))
        k = self.sect(s, f, U)
        nx = {p: x[p] - k[p] for p in U}
        nx.update(self.sect(s, self.query(s, x, T), A))
        nx.update(self.sect(s, f, T))
        ny = {p: y[p] + nx[p] for p in U}
        ny.update({p: x[p] + y[p] - nx[p] for p in S})
        return nx, ny, f

    def backward(self, s, t, x, y):
        qa = [g + h for g, h in zip(self.query(t, y, U), self.query(t, y, T))]
        corr = self.sect(t, qa, A)
        u = dict(y)
        for p in A:
            u[p] = y[p] + corr[p]
        a = self.query(s, x, T)
        b = self.query(s, u, range(112))
        ka = self.sect(s, a, U)
        t11 = self.sect(s, [bi + ai * self.c11 for bi, ai in zip(b, a)], T)
        t12 = self.sect(s, [bi + ai * self.c12 for bi, ai in zip(b, a)], T)
        nx, ny = {}, {}
        for p in U:
            ny[p], nx[p] = y[p] + ka[p], x[p] - y[p]
        for p in A:
            ny[p], nx[p] = u[p], x[p] + y[p] - u[p]
        for p in T:
            ny[p], nx[p] = u[p] - t11[p], t12[p]
        return nx, ny


def random_inputs(model, K, rng):
    m, P = model.m, model.P
    return [[{p: flint.nmod_mat(m, K, [rng.randrange(P) for _ in range(m * K)], P) for p in range(112)}
             for _ in (0, 1)] for _ in range(4)]


def run(model, recipe, inputs, bridge_sign=1, pair_map=None):
    """Simulate all 19 events.  Returns the observation dictionary.
    bridge_sign / pair_map are mutation hooks for controls."""
    st = [[dict(inputs[k][0]), dict(inputs[k][1])] for k in range(4)]
    init = [dict(inputs[k][0]) for k in range(4)]
    hist, cut = {}, {}
    for e in recipe:
        if e['kind'] == 'local':
            k, s = e['pair'], e['stage']
            if pair_map:
                k = pair_map.get(s, k)
            x, y = st[k]
            if e['reverse']:
                nx, ny = model.backward(s, s + 1, x, y)
                if k not in cut:
                    cut[k] = ({p: nx[p] for p in S}, {p: ny[p] for p in S})
            else:
                nx, ny, f = model.forward(s, x, y)
                hist[s] = f
            st[k] = [nx, ny]
        elif e['kind'] == 'bridge':
            i, j, c = e['src'], e['dst'], e['coefficient'] * bridge_sign
            for p in U:
                st[j][0][p] = st[j][0][p] + st[i][0][p] * (c % model.P)
            for p in U:
                st[i][1][p] = st[i][1][p] - st[j][1][p] * (c % model.P)
        else:
            for k in range(4):
                for p in U:
                    st[k][0][p], st[k][1][p] = st[k][1][p], -st[k][0][p]
    obs = {
        'Xf': [{p: st[k][0][p] for p in S} for k in range(4)],      # final selected X
        'Yf': [{p: st[k][1][p] for p in S} for k in range(4)],      # final selected Y
        'C': [{p: st[k][1][p] for p in U} for k in range(4)],       # final Y ordinary (completed-X descendants)
        'f': hist,                                                  # forward histories f0,f2,f4,f6,f8
        'INIT': {(k, s): model.B(s, init[k]) for k in range(4) for s in INIT_STAGES[k]},
        'dA': {k: {p: cut[k][1][p] for p in A} for k in cut},       # Y_A right after each bank's first B
        'cut': cut,
        # cross-bank initial reads used by the alternative contract (pre-run, like INIT64)
        'XB': {0: model.B(2, init[1]), 1: model.B(4, init[2]), 2: model.B(6, init[3])},
        'final_state': st,
    }
    return obs


# ---- small vector helpers ----------------------------------------------------------------------
def vadd(a, b):
    return [x + y for x, y in zip(a, b)]


def vsub(a, b):
    return [x - y for x, y in zip(a, b)]


def dadd(a, b):
    return {p: a[p] + b[p] for p in a}


def dsub(a, b):
    return {p: a[p] - b[p] for p in a}


def dscale(a, c):
    return {p: a[p] * c for p in a}
