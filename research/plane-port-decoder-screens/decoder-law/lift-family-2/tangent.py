"""Tangent dimension of the masked-completion variety {(X, Y): x_a.y_b = delta_ab on required pairs (a = b or a perp b)}
at each saved completion (sols_0_20.npz), minus the GL9 x torus orbit dimension = moduli tangent dimension.
Prints singular-value gaps so that numerical ranks are visible."""
import sys, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from e8 import *
D = np.load('./sols_0_20.npz')
REQ = [(a, b) for a in range(NL) for b in range(NL) if a == b or IP[a][b] == 0]
def jac(X, Y):
    J = np.zeros((len(REQ), 2 * NL * 9))
    for r, (a, b) in enumerate(REQ):
        J[r, a * 9:(a + 1) * 9] += Y[b]
        J[r, NL * 9 + b * 9: NL * 9 + (b + 1) * 9] += X[a]
    return J
def orbit(X, Y):
    V = []
    for i in range(9):
        for j in range(9):
            E = np.zeros((9, 9)); E[i, j] = 1
            V.append(np.concatenate([(X @ E).ravel(), (-Y @ E.T).ravel()]))
    for a in range(NL):
        dX = np.zeros_like(X); dY = np.zeros_like(Y); dX[a] = X[a]; dY[a] = -Y[a]
        V.append(np.concatenate([dX.ravel(), dY.ravel()]))
    return np.array(V)
def nrank(M, tag):
    sv = np.linalg.svd(M, compute_uv=False); sv = sv / sv[0]
    r = int((sv > 1e-8).sum()); gap = (sv[r - 1], sv[r] if r < len(sv) else 0.0)
    return r, gap
if __name__ == '__main__':
  keys = sys.argv[1:] or sorted(D.keys(), key=lambda k: tuple(map(int, k.split('_'))))
  for k in keys:
      V, W = D[k]
      # normalise scale for conditioning: GL9 balance
      J = jac(V, W); rJ, gJ = nrank(J, 'J'); tdim = 2 * NL * 9 - rJ
      O = orbit(V, W); rO, gO = nrank(O, 'O')
      assert np.abs(J @ O.T).max() < 1e-6 * np.abs(J).max() * np.abs(O).max()
      rX = [np.linalg.matrix_rank(M, 1e-7) for M in (V, W)]
      print(k, 'tangent', tdim, f'(rank J {rJ}, sv gap {gJ[0]:.1e}/{gJ[1]:.1e})', 'orbit', rO, f'(gap {gO[0]:.1e}/{gO[1]:.1e})', 'moduli', tdim - rO, flush=True)
