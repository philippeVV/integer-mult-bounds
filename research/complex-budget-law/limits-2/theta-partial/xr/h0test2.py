import numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(1)
def H0(X, Y):
    I = np.eye(len(X)); return np.block([[4*Y, 3*(I+Y@X)], [3*(I+X@Y), 4*X]])
def mineig(M): return abs(np.linalg.eigvalsh((M+M.T)/2)).min()
n = 2
def mk(p, cap):
    th1, th2, a1, a2, b1, b2 = p
    R = lambda t: np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    X = R(th1) @ np.diag(cap*np.tanh([a1, a2])) @ R(th1).T; Y = R(th2) @ np.diag(cap*np.tanh([b1, b2])) @ R(th2).T; return X, Y
for cap in [1.0, 0.999, 0.99, 0.9]:
    best = min((minimize(lambda p: np.log(mineig(H0(*mk(p, cap)))), rng.standard_normal(6)*2, method='Nelder-Mead', options=dict(maxiter=6000, xatol=1e-12, fatol=1e-14)) for _ in range(20)), key=lambda r: r.fun)
    X, Y = mk(best.x, cap); print('cap', cap, 'min|eig|', np.exp(best.fun), 'eigX', np.linalg.eigvalsh(X).round(4), 'eigY', np.linalg.eigvalsh(Y).round(4), 'angle', round((best.x[0]-best.x[1]) % np.pi, 4))
