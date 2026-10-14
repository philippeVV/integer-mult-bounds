import os, sys, pickle
_ENG = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'engine')
sys.path.insert(0, _ENG)
from model import L, GM, F0, LIVE, N
OPEN = pickle.load(open(os.path.join(_ENG, 'k6_open.pkl'),'rb'))
def split(A, c):
    two = [a for a, k in zip(A, c) if k == 2]
    one = [a for a, k in zip(A, c) if k == 1]
    return two, one
def classes(two):
    E6 = [r for r in range(N) if all(GM[r][a] == 0 for a in two)]
    cl = {a: [r for r in range(N) if GM[r][a] == 0 and r not in E6 and r != a and r not in two] for a in two}
    return E6, cl
if __name__ == '__main__':
    for A, c in OPEN:
        two, one = split(A, c)
        E6, cl = classes(two)
        info = []
        for a in two:
            X = cl[a]
            live = [r for r in X if r not in F0]
            n0 = [r for r in live if all(GM[r][d] != 0 for d in one)]
            info.append((a, len(X), len(live), len(n0)))
        print(A, c, 'E6', len(E6), 'E6live', len([r for r in E6 if r not in F0]), 'A2live', [a not in F0 for a in two], info, 'one in F0', [d in F0 for d in one])
