"""For fixed L: min over lattice gauges (8-term decoders N = (G/2)|_L in basis M / dual) of the cost.
Starts: best of 2025 orthogonal frames, simple-root basis, then anneal over rays."""
import sys, random, math
from catalog import *
from gsearch import rand_ray, indep
def best_frame(L):
    return min((gauge_cost([ROOT[a] for a in F], L), F) for F in frames())
def anneal(L, M, iters, rnd):
    M = list(M); cur = gauge_cost(M, L); best = (cur, list(M)); T = 1.0
    Ls = list(bits(L))
    for it in range(iters):
        k = rnd.randrange(8); old = M[k]
        u = rnd.random()
        if u < 0.4: M[k] = ROOT[rnd.randrange(NL)]
        elif u < 0.7: M[k] = nullvec_from(L, rnd)
        else: M[k] = rand_ray(rnd)
        if M[k] is None or not indep(M): M[k] = old; continue
        new = gauge_cost(M, L)
        if new <= cur or rnd.random() < math.exp((cur - new) / T): cur = new
        else: M[k] = old
        T = max(0.05, T * 0.999)
        if cur < best[0]: best = (cur, list(M))
    return best
from gsearch import nullvec
def nullvec_from(L, rnd):
    Ls = list(bits(L)); return nullvec([ROOT[rnd.choice(Ls)] for _ in range(7)])
if __name__ == '__main__':
    which = sys.argv[1:] if len(sys.argv) > 1 else None
    C = catalog(); rnd = random.Random(7)
    for nm, L in C.items():
        if which and nm not in which: continue
        n = bin(L).count('1')
        fc, F = best_frame(L)
        sc = gauge_cost(list(SIMPLE), L); wc = gauge_cost(list(OMEGA), L)
        best = min(fc, sc, wc)
        for s in range(3):
            b, M = anneal(L, [ROOT[a] for a in F], 1500, rnd); best = min(best, b)
        print(f'{nm:22s} |L|={n:3d} frame {fc:3d} simple {sc:3d} weights {wc:3d} anneal-best {best:3d}  gain {n - best:3d}  |L|-48={n - 48}', flush=True)
