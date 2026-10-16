"""Analytic gauge Hessian of Theta_3 at a two-sided flat code (girth >= 8):
  Theta_3(D_y) - 1/8 = (4/E) y^T Gam y + O(|y|^3),
  Gam = (1/4) sym[ (3/16)(I-X)(I-Y) + (1/8)(I-X)Y(I-X) + (1/8)(I-Y)X(I-Y) + (3/16)(I-Y)XY(I-X) ],
  X = 4(|A|^{o2} - I/4) (point blocks), Y = 4(|B|^{o2} - I/4) (line blocks): symmetric, zero diagonal, row sums 1.
lambda_min(Gam | 1^perp) < 0  =>  every code near C0 has Theta_3^* < 1/8 (gauge descent beats any code ascent at 2nd order)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'xg'))
import xg_core as D
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
from xk_sparse import load_any, orth, Struct, gauge_value_grad

def XY_from_codes(V, nP, E, Qs):
    X = np.zeros((E, E)); Y = np.zeros((E, E))
    for i, (es, q) in enumerate(zip(V, Qs)):
        P = q @ q.conj().T; O = 4 * (abs(P) ** 2); np.fill_diagonal(O, 0)
        (X if i < nP else Y)[np.ix_(es, es)] = O
    return X, Y

def Gam(X, Y):
    I = np.eye(len(X)); A = I - X; B = I - Y
    G = (3/16) * A @ B + (1/8) * A @ Y @ A + (1/8) * B @ X @ B + (3/16) * B @ X @ Y @ A
    return (G + G.T) / 8

def lam_min(X, Y):
    E = len(X); G = Gam(X, Y); one = np.ones(E) / np.sqrt(E); Pp = np.eye(E) - np.outer(one, one)
    w, U = np.linalg.eigh(Pp @ G @ Pp + 10 * np.outer(one, one)); return w[0], U[:, 0]

if __name__ == '__main__':
    name = sys.argv[1]; V, nP, E, g = load_any(name); S = Struct(V, nP, E); rng = np.random.default_rng(0)
    Qf = D.fourier_codes(V, rng); X, Y = XY_from_codes(V, nP, E, Qf); l, v = lam_min(X, Y)
    # finite-difference check against the true Theta_3 along v
    for tau in [1e-2, 3e-3]:
        th = gauge_value_grad(S, Qf, 3, tau * v)[0]; print(f'{name} Fourier: analytic (4/E) lam tau^2 = {4/E*l*tau**2:.3e}; true Theta(tau v) - 1/8 = {th-0.125:.3e}')
    print(name, 'Fourier lambda_min(Gam) =', l)
