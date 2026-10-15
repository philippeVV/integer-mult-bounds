"""Enumerate all nonempty orthogonal cliques of P112 (bitmasks over live positions)."""
import pickle, sys
from collections import Counter
from model import NL, ONB, popc

def all_cliques():
    out = []
    def ext(cur, cand):
        out.append(cur)
        c = cand
        while c:
            b = c & -c
            i = b.bit_length() - 1
            c ^= b
            # only extend with higher indices to avoid duplicates
            ext(cur | b, c & ONB[i])
    for a in range(NL):
        higher = ONB[a] & ~((1 << (a + 1)) - 1)
        ext(1 << a, higher)
    return out

if __name__ == "__main__":
    cl = all_cliques()
    dist = Counter(popc(q) for q in cl)
    print(len(cl), sorted(dist.items()))
    pickle.dump(cl, open("cliques.pkl", "wb"))
