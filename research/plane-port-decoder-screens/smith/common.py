import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gain'))
from base import *
EV = [set(b for b in range(n) if b != a and dot(P[a], P[b]) == 1) for a in range(n)]
def cliques():
    cl = []
    def bk(Rr, Pp, X):
        if not Pp and not X: cl.append(frozenset(Rr)); return
        u = max(Pp | X, key=lambda v: len(EV[v] & Pp))
        for v in list(Pp - EV[u]):
            bk(Rr | {v}, Pp & EV[v], X & EV[v]); Pp = Pp - {v}; X = X | {v}
    bk(set(), set(range(n)), set())
    return cl
CL = cliques()
CL8 = sorted([sorted(c) for c in CL if len(c) == 8])
def linerank(R, A):
    return pair_rank([P[a] for a in R], [P[b] for b in A])
