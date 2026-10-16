"""expXH part 2: Hermitian (Gaussian / Eisenstein) port systems with the one-valuation condition.
Ports v with H(v,v)=2; labels L_v = O-span of (v,1) mod 2 (F2-dim u=2); free iff H_ab in 1+2O_(2).
Decoder N = H/2 (complex), valid iff every nonzero off-diagonal H_ab lies in 1+2O_(2).
Checks: (i) TENSOR LAW: the explicit 2V x 2V F2 label Gram B_C equals (gauge) B_R (x) T, so h_C = 2 h_R and every term's
pairing rank doubles; (ii) max valid sets in Gaussian / Eisenstein E8, complex rank, gamma_C with complex terms.
usage: part2.py"""
import itertools, math
MAXV = []
import numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from lib import rank_f2_rows, bits_of, max_clique

# ---- E8 roots (real, standard coords) ----
R = []
for i, j in itertools.combinations(range(8), 2):
    for si, sj in itertools.product((1, -1), repeat=2):
        v = [0.0] * 8; v[i] = si; v[j] = sj; R.append(v)
for s in itertools.product((0.5, -0.5), repeat=8):
    if sum(1 for t in s if t < 0) % 2 == 0: R.append(list(s))
R = np.array(R); assert len(R) == 240
lines = R[[next(t for t in r if t != 0) > 0 for r in R]]
assert len(lines) == 120


def f2rank(M):
    return rank_f2_rows(bits_of(np.asarray(M, dtype=np.uint8)))


def analyse(name, X, Hfun, resid_pairs, u_basis, trace_form):
    """X: candidate real vectors; Hfun(x,y)->complex H; resid: maps complex H -> (is_zero, in_1+2O, residue pair mod 2)."""
    V0 = len(X)
    H = np.array([[Hfun(x, y) for y in X] for x in X])
    assert np.allclose(np.diag(H), 2)
    isz = np.isclose(np.abs(H), 0)
    good = np.zeros_like(isz)
    for a in range(V0):
        for b in range(V0):
            if a != b and not isz[a, b]:
                good[a, b] = resid_pairs(H[a, b])[1]
    adj = isz | good; np.fill_diagonal(adj, False)
    cl, opt, ub = max_clique(adj, seconds=60)
    S = np.array(cl); V = len(S)
    Hs = H[np.ix_(S, S)]
    rC = np.linalg.matrix_rank(Hs, tol=1e-8)
    rRe = np.linalg.matrix_rank(Hs.real, tol=1e-8)
    # real pattern labels: B_R = 1 iff H=0 or diagonal
    BR = (np.isclose(np.abs(Hs), 0) | np.eye(V, dtype=bool)).astype(np.uint8)
    hR = f2rank(BR)
    # explicit F2 label Gram: for ports a,b and basis elements lam, mu of O/2 (u_basis), entry = trace_form(lam*conj(mu)*z_ab)
    # with z = H + 1 (label (v,1)), reduced mod 2
    BC = np.zeros((2 * V, 2 * V), dtype=np.uint8)
    for a in range(V):
        for b in range(V):
            z = Hs[a, b] + 1
            for i, lam in enumerate(u_basis):
                for j, mu in enumerate(u_basis):
                    BC[2 * a + i, 2 * b + j] = trace_form(lam * np.conj(mu) * z)
    hC = f2rank(BC)
    assert hC == 2 * hR, (hC, hR)   # tensor law (limits-1 replay)
    MAXV.append(V)
    print(f"{name}: candidates={V0} max valid set V={V} (optimal={opt}, ub={ub}); complex rank={rC}, rank Re H={rRe}; "
          f"h_R={hR}, h_C={hC} (tensor law 2h_R: {hC == 2 * hR})")
    return S, Hs, BR, BC, V, hR, hC, rC


# Gaussian structure on E8: i acts on coordinate pairs (x0,x1)->(-x1,x0)
def Jg(x):
    y = np.empty(8); y[0::2] = -x[1::2]; y[1::2] = x[0::2]; return y
def Hg(x, y):
    return complex(x @ y, x @ Jg(y))   # Re = b(x,y); sesquilinear in the i-structure (convention-free for the checks)
for x in lines[:20]:
    for y in lines[:20]:
        assert abs(Hg(Jg(x), y) - 1j * Hg(x, y)) < 1e-9 or abs(Hg(Jg(x), y) + 1j * Hg(x, y)) < 1e-9
def resid_g(Hv):
    a, b = round(Hv.real), round(Hv.imag); assert abs(Hv - complex(a, b)) < 1e-9
    return (a == 0 and b == 0), (a % 2 == 1 and b % 2 == 0), None
def tr_g(z):   # F2 form on Z[i]/2: Re(z) mod 2
    a = round(z.real); assert abs(z.real - a) < 1e-9; return a % 2
Sg = analyse('Gaussian E8', lines, Hg, resid_g, [1, 1j], tr_g)

# complex terms for Gaussian: complex coordinates z_j = x_{2j} + i x_{2j+1} -> conj-consistent with Hg up to convention
def zc(x):
    return x[0::2] - 1j * x[1::2]
S, Hs, BR, BC, V, hR, hC, rC = Sg
Z = np.array([zc(x) for x in lines[S]])
Hrec = Z @ Z.conj().T
assert np.allclose(Hrec, Hs) or np.allclose(Hrec, Hs.conj()), 'coordinate form mismatch'
CR = CC = 0
for c in range(4):
    Sc = np.nonzero(~np.isclose(np.abs(Z[:, c]), 0))[0]
    r = f2rank(BR[np.ix_(Sc, Sc)])
    rc = f2rank(BC[np.ix_(np.repeat(2 * Sc, 2) + np.tile([0, 1], len(Sc)), np.repeat(2 * Sc, 2) + np.tile([0, 1], len(Sc)))])
    assert rc == 2 * r, (rc, r)
    CR += r; CC += rc
print(f"  Gaussian E8 complex coordinate decoder (4 terms): C_C={CC} = 2*{CR} (term ranks double: OK); u=2, "
      f"gamma_C = C_C/(uV) = {CC / (2 * V):.3f}; merit gamma*(h-1)/u = {CC / (2 * V) * (hC - 1) / 2:.2f}")

# Eisenstein structure: sigma = product of A2 Coxeter elements on 4 mutually orthogonal A2's in E8
def refl(a):
    return np.eye(8) - np.outer(a, a)  # roots have norm 2: s_a = I - a a^T (2/(a.a)=1)
A2s = []
used = []
for a in lines:
    if any(abs(a @ u_) > 1e-9 for u_ in used): continue
    for b in R:
        if abs(a @ b + 1) < 1e-9 and all(abs(b @ u_) < 1e-9 for u_ in used):
            A2s.append((a, b)); used += [a, b]; break
    if len(A2s) == 4: break
assert len(A2s) == 4
sig = np.eye(8)
for a, b in A2s:
    sig = sig @ refl(a) @ refl(b)
assert np.allclose(np.linalg.matrix_power(sig, 3), np.eye(8)) and np.allclose(sig @ sig + sig + np.eye(8), 0)
s3 = math.sqrt(3)
def He(x, y):   # Hermitian with Re = b, H(sig x, y) = omega H(x,y)  (verified below)
    return complex(x @ y, (2 * (x @ (sig @ y)) + x @ y) / s3)
om = complex(-0.5, s3 / 2)
for x in lines[:15]:
    for y in lines[:15]:
        h0 = He(x, y)
        assert abs(He(sig @ x, y) - om * h0) < 1e-9 or abs(He(sig @ x, y) - np.conj(om) * h0) < 1e-9
def to_po(Hv):    # Hv = p + q*omega, p,q rational with denominators dividing 3
    q = 2 * Hv.imag / s3; p = Hv.real + q / 2
    P, Q = round(3 * p), round(3 * q); assert abs(3 * p - P) < 1e-7 and abs(3 * q - Q) < 1e-7
    return P, Q   # 3H = P + Q omega
def resid_e(Hv):
    P, Q = to_po(Hv)
    # H in 1+2O_(2)  <=>  3H - 3 in 2 O  <=>  P-3 even and Q even
    return (P == 0 and Q == 0), ((P - 3) % 2 == 0 and Q % 2 == 0), None
def tr_e(z):     # F2 form on O/2 = F4: Tr(z) = 2Re z mod 2, with 2-adic units (1/3) allowed
    t = 2 * z.real * 3; T = round(t); assert abs(t - T) < 1e-7; return (T * 1) % 2   # 3*Tr(z), 3 odd
Se = analyse('Eisenstein E8', lines, He, resid_e, [1, om], tr_e)
S, Hs, BR, BC, V, hR, hC, rC = Se
# complex coordinates: an orthonormal-ish basis of the omega-eigenspace of sig; generic basis -> full supports
w, U = np.linalg.eig(sig)
E = U[:, np.isclose(w, om) | np.isclose(w, np.conj(om))]
best = None
# terms: support-minimising choice among coordinate bases is not attempted; report generic (all-support) cost bound
r_full = f2rank(BR)
print(f"  Eisenstein E8: complex rank {rC} -> generic {rC} full-support terms, C_C = {rC} * 2 * {r_full} = {2 * rC * r_full}, "
      f"gamma_C = {rC * r_full / V:.3f}")

assert MAXV == [12, 4], MAXV   # Gaussian / Eisenstein E8: max valid sets 12 and 4 (CP-SAT optimal)
print('PASS_EXPXH_PART2')
