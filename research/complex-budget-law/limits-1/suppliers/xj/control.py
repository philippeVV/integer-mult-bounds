"""Control: can the phase-frame (beta) model express a saving at m <= 3 at all?
Exhaustive min s on expNJ's flow-deficient (unrealisable) W = 6 instances, where rational frames give 10-11 < 12 at m = 2."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # python -I omits the script dir
import beta as B

RES = []
INST = [([(0, 1), (0, 2, 3, 4), (0, 1, 2, 5)], [3, 4, 0, 1, 2, 5]),
        ([(0, 1, 2, 3), (0, 1, 4, 5), (0, 2)], [1, 3, 4, 5, 0, 2])]
for gates, rho in INST:
    W = 6; first, last, segs = {}, {}, []
    for gi, g in enumerate(gates):
        for w in g:
            if w in last: segs.append((last[w], gi))
            else: first[w] = gi
            last[w] = gi
    for m in (1, 2, 3):
        nb, rk, cost, I = B.tables(m)
        for name, tab in (('rank', rk), ('cost', cost)):
            best = None
            for fr in itertools.product(range(1 << nb), repeat=len(gates) - 1):
                F = (0,) + fr
                s = sum(tab[F[a] ^ F[b]] for a, b in segs)
                s += sum(tab[F[first[w]] ^ F[last[rho[w]]] ^ I] for w in range(W))
                if best is None or s < best[0]: best = (s, F)
            print(f"{gates} rho={rho} m={m} {name}: min s = {best[0]}  Wm = {W*m}  frames={best[1]}")
            RES.append((m, best[0]))
# limits-1 replay: on both unrealisable instances the model expresses s = 4m < Wm = 6m (down to the flow bound m F_*)
assert all(s == 4 * m for m, s in RES), RES
print('PASS_EXPXJ_CONTROL')
