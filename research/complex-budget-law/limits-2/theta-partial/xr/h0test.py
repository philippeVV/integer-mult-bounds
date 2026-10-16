import numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(0)
def H0(X, Y):
    I = np.eye(len(X)); return np.block([[4*Y, 3*(I+Y@X)], [3*(I+X@Y), 4*X]])
def inert(M, tol=1e-10):
    ev = np.linalg.eigvalsh((M+M.T)/2); return (int((ev>tol).sum()), int((ev<-tol).sum())), abs(ev).min()
def symc(p, n):   # symmetric contraction via eigen-parametrisation
    A = p[:n*n].reshape(n, n); Q = np.linalg.qr(A)[0]; e = np.tanh(p[n*n:n*n+n]) ; return (Q*e)@Q.T
worst = 1e9
for n in [2, 3, 4, 6]:
    for t in range(300):
        p = rng.standard_normal(2*(n*n+n))*2
        X = symc(p[:n*n+n], n); Y = symc(p[n*n+n:], n)
        (ip, im), mn = inert(H0(X, Y))
        worst = min(worst, mn)
        if (ip, im) != (n, n): print('INERTIA CHANGE', n, ip, im)
    # adversarial: minimise smallest |eig|
    f = lambda p: np.log(inert(H0(symc(p[:n*n+n], n), symc(p[n*n+n:], n)))[1])
    best = min((minimize(f, rng.standard_normal(2*(n*n+n)), method='Nelder-Mead', options=dict(maxiter=4000)) for _ in range(5)), key=lambda r: r.fun)
    print(n, 'random min|eig|', worst, ' adversarial min|eig|', np.exp(best.fun))
