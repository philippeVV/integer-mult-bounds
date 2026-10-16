"""(3) Signed norm-2 label families on the PG(2,16)-M 784-arc schedule with 12 C/J lanes: nullity of Z_G mod p.
Frames: source M_in(k) = -P_k, sink M_out(rho k) = I - P_k, P_k = l l^T H / nrm (denominator 2 or 4); rank-r: sum over
r mutually H-orthogonal labels. usage: families.py gseed mode m_list fam_list"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import sys, json, time, itertools, numpy as np
from circ import *

def droots(m):
    R = []
    for a, b in itertools.combinations(range(m), 2):
        for sa, sb in [(1, 1), (1, -1)]:
            v = np.zeros(m, dtype=np.int64); v[a] = sa; v[b] = sb; R.append(v)
    return R

def e8roots():
    R = droots(8)
    for signs in itertools.product([1, -1], repeat=8):
        if signs.count(-1) % 2 == 0: R.append(np.array(signs, dtype=np.int64))  # (1/2)^8 scaled by 2: norm 8
    return R  # D8 part: norm 2 (H=I, den 2); half part stored x2: norm 8 -> P = l l^T / 8

def hamming8():
    # extended Hamming [8,4,4]: 14 weight-4 words; labels +-1 on support, H = I/2 -> norm 2, P = l l^T/4
    G = np.array([[1,1,1,1,0,0,0,0],[0,0,1,1,1,1,0,0],[0,0,0,0,1,1,1,1],[0,1,0,1,0,1,0,1]])
    words = set()
    for c in itertools.product([0, 1], repeat=4):
        w = tuple(np.array(c) @ G % 2)
        if sum(w) == 4: words.add(w)
    return sorted(words)

def proj(ls, p=P):
    """sum of l l^T / (l.l) for mutually orthogonal integer labels -> mod p matrix (rational idempotent)."""
    m = len(ls[0]); M = np.zeros((m, m), dtype=np.int64)
    for l in ls:
        n = int(l @ l); inv = pow(n, p - 2, p)
        M = (M + (np.outer(l, l) % p) * inv) % p
    return M

def realify_c(z):
    """complex vector z (list of complex ints) -> real 2r vector; complex line proj realified = rank-2 real proj
    spanned by realify(z), realify(i z) (orthogonal, same norm)."""
    a = np.array([int(c.real) for c in z] + [int(c.imag) for c in z], dtype=np.int64)
    iz = [1j * c for c in z]
    b = np.array([int(c.real) for c in iz] + [int(c.imag) for c in iz], dtype=np.int64)
    return [a, b]

def conflict_graph(c):
    """sessions sharing a source gate or a sink gate."""
    W = c.W; gs = {}; 
    for k, a in enumerate(c.src):
        u, w = c.arcs[a]; gs.setdefault(('s', w), []).append(k)
        u2, _ = c.arcs[c.rho[a]]; gs.setdefault(('t', u2), []).append(k)
    nb = [set() for _ in range(W)]
    for ks in gs.values():
        for x in ks:
            nb[x].update(y for y in ks if y != x)
    return nb

def greedy_orth(c, roots, rng, nb):
    W = c.W; lab = [None] * W; bad = 0
    for k in rng.permutation(W):
        best = None; bc = 1e9
        for t in rng.permutation(len(roots))[:200]:
            l = roots[t]; cost = sum(1 for y in nb[k] if lab[y] is not None and int(lab[y] @ l) != 0)
            if cost < bc: bc, best = cost, l
            if cost == 0: break
        lab[k] = best; bad += (bc > 0)
    return lab, bad

def family(name, m, c, rng, nb):
    W = c.W; info = {}
    if name == 'F1_Droot_r1':
        R = droots(m); Ps = [proj([R[rng.integers(len(R))]]) for _ in range(W)]
    elif name == 'F2_Droot_gateorth':
        lab, bad = greedy_orth(c, droots(m), rng, nb); Ps = [proj([l]) for l in lab]; info['conflicts'] = bad
    elif name.startswith('F3_Droot_rank'):
        r = m // 2
        # r mutually orthogonal D_m roots: pair up coordinates, choose sign/partner per pair
        Ps = []
        for _ in range(W):
            perm = rng.permutation(m); ls = []
            for t in range(r):
                a, b = perm[2 * t], perm[2 * t + 1]; v = np.zeros(m, dtype=np.int64); v[a] = 1; v[b] = rng.choice([1, -1]); ls.append(v)
            Ps.append(proj(ls))
        info['prank'] = r
    elif name == 'F4_E8_r1':
        assert m == 8; R = e8roots(); Ps = [proj([R[rng.integers(len(R))]]) for _ in range(W)]
    elif name == 'F4b_E8_rank4':
        assert m == 8; R = e8roots(); Ps = []
        for _ in range(W):
            ls = []
            for t in rng.permutation(len(R)):
                l = R[t]
                if all(int(l @ x) == 0 for x in ls): ls.append(l)
                if len(ls) == 4: break
            Ps.append(proj(ls))
    elif name == 'F5_F16_hamming':
        assert m == 8; H = hamming8()
        # F16-derived: session arc (P_i, L_j) -> a = (i*j + j) mod 15 picks a Hamming word (+ sign from i parity)
        Ps = []
        for k, a in enumerate(c.src):
            u, w = c.arcs[a]; x = (u * w + w) % 14; sg = 1 if (u + w) % 2 == 0 else -1
            l = np.array(H[x], dtype=np.int64); idx = np.nonzero(l)[0]; l[idx[-1]] *= sg
            Ps.append(proj([l]))
    elif name == 'F5b_hamming_gateorth':
        assert m == 8; Hw = hamming8(); roots = []
        for w in Hw:
            idx = np.nonzero(w)[0]
            for sg in itertools.product([1, -1], repeat=3):
                v = np.array(w, dtype=np.int64); v[idx[1:]] *= np.array(sg); roots.append(v)
        lab, bad = greedy_orth(c, roots, rng, nb); Ps = [proj([l]) for l in lab]; info['conflicts'] = bad
    elif name == 'F6_herm_realified':
        assert m % 2 == 0; h = m // 2
        Ps = []
        for _ in range(W):
            z = [0j] * h
            if h == 1 or rng.random() < 0.3:
                z[rng.integers(h)] = complex(1, rng.choice([1, -1]))   # (1+-i) e_a, hermitian norm 2
            else:
                a, b = rng.choice(h, 2, replace=False); z[a] = 1; z[b] = [1, -1, 1j, -1j][rng.integers(4)]
            Ps.append(proj(realify_c(z)))
    elif name == 'F7_twoclass':
        R = droots(m); l = R[0]; Pm = proj([l]); I = np.eye(m, dtype=np.int64)
        Ps = [Pm if rng.random() < 0.5 else (I - Pm) % P for _ in range(W)]
    elif name == 'R_random_control':
        Ps = [rng.integers(-2, 3, size=(m, m)) for _ in range(W)]
    elif name == 'F8_common_control':
        R = droots(m); Pm = proj([R[0]]); Ps = [Pm] * W
    else:
        raise ValueError(name)
    Ps = np.array(Ps) % P
    I = np.eye(m, dtype=np.int64)
    Min = (-Ps) % P; Mout = (I[None] - Ps) % P
    return Min, Mout, info

FAMS = ['F1_Droot_r1', 'F2_Droot_gateorth', 'F3_Droot_rankhalf', 'F6_herm_realified', 'F7_twoclass', 'F8_common_control']
FAMS8 = ['F4_E8_r1', 'F4b_E8_rank4', 'F5_F16_hamming', 'F5b_hamming_gateorth']

if __name__ == '__main__':
    gseed = int(sys.argv[1]); mode = sys.argv[2]
    ms = [int(x) for x in sys.argv[3].split(',')]
    fams = sys.argv[4].split(',') if len(sys.argv) > 4 else None
    c = Circuit(mode=mode, seed=0); nb = conflict_graph(c)
    A = c.transfer(np.random.default_rng(1000 + gseed))
    out = []
    for m in ms:
        fl = fams or (FAMS + (FAMS8 if m == 8 else []))
        for f in fl:
            if f in ('F3_Droot_rankhalf', 'F6_herm_realified') and m % 2: continue
            rng = np.random.default_rng(7 * gseed + m)
            Min, Mout, info = family(f, m, c, rng, nb)
            t = time.time(); r = rank_mod(Zmat(A, Min, Mout))
            rec = dict(gseed=gseed, mode=mode, m=m, fam=f, rank=r, Wm=c.W * m, nullity=c.W * m - r, sec=round(time.time() - t, 1), **info)
            print(json.dumps(rec), flush=True); out.append(rec)
    json.dump(out, open(_os.path.join(HERE, f'fam_{mode}_g{gseed}_m{sys.argv[3].replace(",","-")}.json'), 'w'))
    print('PASS_EXPQR_FAMILIES' if out and max(r['nullity'] for r in out) == 0 else 'FAIL_EXPQR_FAMILIES')
