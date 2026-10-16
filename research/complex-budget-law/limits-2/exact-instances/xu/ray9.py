"""expXU: expQR's ray test (all scalings t of a frame direction N0, any denominator) on the exact lane-free q = 9 schedule.
Rational-frame (contract) model: source frames M_in(k) = N0_k, sink frames M_out(k) = I + N0_k (rho = id on sessions).
Z_A(tN0) = (I - t K_A) A~, K_A = A~ N0 A~^{-1} - N0.  f_A = charpoly(K_A) mod p, x-power stripped; trivial gcd over 3
independent random char-0 realisations A => no scaling of N0 survives GRS (s >= Wm for every interior frame, every t).
Also records the fixed-scale (t = 1) nullity of Z for draw 0.
Families: expQR's signed norm-2 label families (F1..F8, R) and expXN's 23 F2 families lifted to {0,1} rational matrices.
usage: ray9.py m [qr|xn|both]"""
import os, sys, json, time
import numpy as np, flint
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import circ9 as c
sys.path.insert(1, os.path.join(HERE, '..', 'expQR'))
import families as QR                      # expQR/families.py (uses expQR/circ.py helpers only; P identical)
assert QR.P == c.P
P = c.P; W = c.W


class Adapter:                             # what expQR's family() / conflict_graph() read from a Circuit
    W = c.W
    src = list(c.R)
    arcs = c.ARCS
    rho = {a: a for a in c.R}


def ray(N0, As):
    m = N0.shape[1]; polys = []; I = np.eye(m, dtype=np.int64)
    for A in As:
        Ai = np.array(flint.nmod_mat(A.tolist(), P).inv().entries(), dtype=np.int64).reshape(W, W)
        AN = (A[:, :, None, None] * N0[None, :, :, :]) % P
        AN = AN.transpose(0, 2, 1, 3).reshape(W * m, W * m)
        K = flint.nmod_mat(AN.tolist(), P) * flint.nmod_mat(np.kron(Ai, I).tolist(), P)
        K = np.array(K.entries(), dtype=np.int64).reshape(W * m, W * m)
        for k in range(W): K[k*m:(k+1)*m, k*m:(k+1)*m] = (K[k*m:(k+1)*m, k*m:(k+1)*m] - N0[k]) % P
        co = [int(x) for x in flint.nmod_mat(K.tolist(), P).charpoly().coeffs()]; z = 0
        while co[z] == 0: z += 1
        polys.append(flint.nmod_poly(co[z:], P))
    g = polys[0]
    for f in polys[1:]: g = g.gcd(f)
    pairs = [polys[i].gcd(polys[j]).degree() for i in range(3) for j in range(i + 1, 3)]
    return [f.degree() for f in polys], g.degree(), pairs


def null_t1(N0, A):
    m = N0.shape[1]; I = np.eye(m, dtype=np.int64)
    Mout = (I[None] + N0) % P
    D = (Mout[:, None] - N0[None, :]) % P
    Zm = ((A[:, :, None, None] * D) % P).transpose(0, 2, 1, 3).reshape(W * m, W * m)
    return W * m - flint.nmod_mat(Zm.tolist(), P).rank()


if __name__ == '__main__':
    m = int(sys.argv[1]); which = sys.argv[2] if len(sys.argv) > 2 else 'both'
    As = [c.transfer_p(np.random.default_rng(5000 + g)) for g in range(3)]
    print(json.dumps(dict(ranks=[flint.nmod_mat(A.tolist(), P).rank() for A in As],
                          nonzeros=[int((A != 0).sum()) for A in As])), flush=True)
    fams = []
    if which in ('qr', 'both'):
        ad = Adapter(); nb = QR.conflict_graph(ad)
        names = QR.FAMS + ['R_random_control'] + (QR.FAMS8 if m == 8 else [])
        for f in names:
            if f in ('F3_Droot_rankhalf', 'F6_herm_realified') and m % 2: continue
            Min, Mout, info = QR.family(f, m, ad, np.random.default_rng(7 * 0 + m), nb)
            fams.append(('QR:' + f, Min % P, info))         # N0 = M_in (= -P_k); M_out = I + N0 as in expQR
    if which in ('xn', 'both'):
        import screen9
        F = screen9.families(m)
        for name, betas in F.items():
            if isinstance(betas, list): fams.append(('XN:' + name, np.array(betas, dtype=np.int64) % P, {}))
    for name, N0, info in fams:
        t = time.time(); degs, gdeg, pairs = ray(N0, As); n1 = null_t1(N0, As[0])
        rec = dict(m=m, fam=name, W=W, deg=degs, gcd_deg_all=gdeg, gcd_deg_pairs=pairs, nullity_t1_draw0=n1,
                   sec=round(time.time() - t, 1), **{k: v for k, v in info.items() if isinstance(v, (int, float, str))})
        print(json.dumps(rec), flush=True)
        with open(os.path.join(HERE, 'ray_results.jsonl'), 'a') as fh: fh.write(json.dumps(rec) + '\n')
