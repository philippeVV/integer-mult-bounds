"""Smallest x=log2 p at which each 'eventually' side condition of the support-local H call holds at PR2044's
exponents (eps=1-2u, beta=1-6u, u=1e-30010), unit constants (Theta(.)=1). Fixed point x = log2(f(x))/delta,
computed in log-domain with mpmath. Also the reservation threshold R(h)<=h/16 and the deferred-phase geometric sum."""
import mpmath as mp, json, sys, math
mp.mp.dps = 60
u = mp.mpf(10)**-30010
def solve(delta, f):
    x = mp.mpf(1)
    for _ in range(200): x = mp.log(f(x), 2) / delta
    return x
G = lambda x: 4*x + 6; K = lambda x: 128*G(x)
rows = {
 'axis fits K-neighbourhood above local core: p^(1-eps) > K + log2 L (L=p^128)': solve(2*u, lambda x: K(x) + 128*x),
 'router reach: n=p >= 9kG, k=2gK': solve(6*u, lambda x: 9*2*K(x)*G(x)),
 'short-record gap: n > 2(k+1)ceil(log2 p)^2': solve(6*u, lambda x: 2*(2*K(x)+1)*x*x),
 'all g neighbourhoods fit address: gK <= n': solve(6*u, lambda x: K(x)),
 'group below dimension: g <= d (no log factor)': mp.mpf(0),
}
# reservation: R(h)= ceil(2hG/K)+ceil(hG/K)+q0, K=128G, q0=1 chunk (row stock <= ~60 log2 h bits << K)
th = next(h for h in range(1, 10**6) if all(math.ceil(2*H/128)+math.ceil(H/128)+1 <= H/16 for H in range(h, h+2000)))
# geometric deferred sum with exponent 1-q, q ~ 0.00168
q = 0.0016797549561255
ratio = max((math.ceil(2*h/128)+math.ceil(h/128)+1)/h for h in range(th, 10**6))
geo = 1/(1 - ratio**(1-q))
out = {k: mp.nstr(v, 8) for k, v in rows.items()}
out['reservation R(h)<=h/16 holds for all h >='] = th
out['max R(h)/h above that'] = ratio
out['sum_j h_j^(1-q) / g^(1-q) <='] = geo
json.dump(out, open(sys.argv[1], 'w'), indent=1); print(json.dumps(out, indent=1))
