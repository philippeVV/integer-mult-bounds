"""Block decoders: N = (+)_i (G/2 in gauge M_i)|_{L_i}  (+) I on leftover lines.  Block-diagonal N is always a
valid masked decoder (cross entries 0 are allowed everywhere).  Maximise total gain = sum_i (|L_i| - cost_i)."""
import sys, random, math
from e8 import *
from gsearch import rand_ray, indep, nullvec
def gainb(M, L): return bin(L).count('1') - gauge_cost(M, L) if L else 0
if __name__ == '__main__':
    seed, parts, iters = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rnd = random.Random(seed)
    Fs = frames()
    Ms = [[ROOT[a] for a in Fs[rnd.randrange(2025)]] for _ in range(parts)]
    asg = [0] * NL if len(sys.argv) <= 4 else [rnd.randrange(parts + 1) for _ in range(NL)]   # start: all in block 0 (gain 48)
    def masks(): 
        m = [0] * (parts + 1)
        for a in range(NL): m[asg[a]] |= 1 << a
        return m
    mk = masks(); g = [gainb(Ms[i], mk[i]) for i in range(parts)]
    cur = sum(g); best = cur; T = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    for it in range(iters):
        if rnd.random() < 0.15:
            i = rnd.randrange(parts); k = rnd.randrange(8); old = Ms[i][k]
            Ls = list(bits(mk[i])) or list(range(NL))
            u = rnd.random()
            Ms[i][k] = ROOT[rnd.randrange(NL)] if u < 0.4 else (nullvec([ROOT[rnd.choice(Ls)] for _ in range(7)]) if u < 0.8 else rand_ray(rnd))
            if Ms[i][k] is None or not indep(Ms[i]): Ms[i][k] = old; continue
            ng = gainb(Ms[i], mk[i]); new = cur - g[i] + ng
            if new >= cur or rnd.random() < math.exp((new - cur) / T): cur = new; g[i] = ng
            else: Ms[i][k] = old
        else:
            a = rnd.randrange(NL); o = asg[a]; p = rnd.randrange(parts + 1)
            if p == o: continue
            asg[a] = p; mk[o] ^= 1 << a; mk[p] ^= 1 << a
            go = gainb(Ms[o], mk[o]) if o < parts else 0; gp = gainb(Ms[p], mk[p]) if p < parts else 0
            new = cur - (g[o] if o < parts else 0) - (g[p] if p < parts else 0) + go + gp
            if new >= cur or rnd.random() < math.exp((new - cur) / T):
                cur = new
                if o < parts: g[o] = go
                if p < parts: g[p] = gp
            else: asg[a] = o; mk[o] ^= 1 << a; mk[p] ^= 1 << a
        T = max(0.15, T * 0.9997)
        if cur > best:
            best = cur; print('it', it, 'gain', best, 'blocks', [(bin(mk[i]).count('1'), g[i]) for i in range(parts)], flush=True)
    print('final best', best, [(bin(mk[i]).count('1'), g[i]) for i in range(parts)])
