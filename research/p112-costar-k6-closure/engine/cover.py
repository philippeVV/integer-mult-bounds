"""Covering sets of normals: sets A of root lines whose co-stars cover all 112 live lines."""
import pickle
from collections import Counter
from model import N, COSTAR, FULL, popc, F0

def covering_sets(k):
    out = []
    def rec(start, chosen, union):
        if len(chosen) == k:
            if union == FULL:
                out.append(tuple(chosen))
            return
        for a in range(start, N):
            rec(a + 1, chosen + [a], union | COSTAR[a])
    rec(0, [], 0)
    return out

if __name__ == "__main__":
    c3 = covering_sets(3)
    print("covering 3-sets:", len(c3))
    # max coverage by 3 co-stars
    best = 0
    for a in range(N):
        for b in range(a + 1, N):
            ab = COSTAR[a] | COSTAR[b]
            for c in range(b + 1, N):
                best = max(best, popc(ab | COSTAR[c]))
    print("max live lines covered by 3 co-stars:", best)
    # 4-sets: prune by requiring union of 3 to leave a remainder coverable by one co-star
    out = []
    for a in range(N):
        for b in range(a + 1, N):
            ab = COSTAR[a] | COSTAR[b]
            for c in range(b + 1, N):
                miss = FULL & ~(ab | COSTAR[c])
                if popc(miss) > 60:
                    continue
                for d in range(c + 1, N):
                    if COSTAR[d] & miss == miss:
                        out.append((a, b, c, d))
    print("covering 4-sets:", len(out))
    print("number of F0 normals per covering 4-set:", sorted(Counter(sum(1 for x in s if x in F0) for s in out).items()))
    pickle.dump(out, open("cover4.pkl", "wb"))
