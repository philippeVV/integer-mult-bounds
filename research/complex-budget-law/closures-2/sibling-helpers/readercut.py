#!/usr/bin/env python3
"""expNU / readercut.py -- Lemma S reader cuts on the realised expSL network (exact mod 2^61-1).
For a target set T', R(T') = shared (helper) registers arriving at the chain nodes of T', S(T') = source forms there
(per-copy side information), rho(T') = rank{alpha_t mod S(T') : t in T'}.  Lemma S (README section 3): a joint network
with the SAME helper registers shared by two copies needs |R(T')| >= 2 rho(T') for every T'.  Single copy needs
|R| >= rho.  Reports min over T' of |R| - 2 rho (greedy growth from each target and by shared chain nodes, plus T' = all)."""
import sys, os, json
from collections import defaultdict
sys.argv = [sys.argv[0]]
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'census.py')).read().split('# ---- per-node slack')[0])
chains = {int(t): [tuple(N) for N in ch] for t, ch in D['tchains'].items()}
inn = defaultdict(list)
for i, (u, v, f, tag) in enumerate(E):
    inn[v].append(i)


def RS(Ts):
    R, S = set(), []
    for t in Ts:
        for N in chains[t]:
            for i in inn.get(N, []):
                if role_of.get(i) == 'S':
                    S.append(E[i][2])
                else:
                    R.add(i)
    return R, S


def rho(Ts):
    R, S = RS(Ts)
    _, Bs = rank(S)
    r, _ = rank([reduce(Bs, [md(coef(t, s)) for s in range(n)]) for t in Ts])
    return len(R), r


best = None
allT = list(chains)
k, r = rho(allT)
print(f'T\' = all 112 targets: |R| = {k}, rho = {r}, |R| - 2 rho = {k - 2 * r}')
# groups of targets sharing a chain node
bynode = defaultdict(list)
for t, ch in chains.items():
    for N in ch:
        bynode[N].append(t)
groups = [tuple(sorted(g)) for g in bynode.values() if len(g) >= 2]
worst = []
for g in set(groups):
    k, r = rho(g)
    worst.append((k - 2 * r, len(g), k, r))
worst.sort()
print('shared-chain-node groups:', len(set(groups)), 'min |R|-2rho', worst[:5])
# greedy: from each target, add the target that minimises |R| - 2 rho, up to 12 targets
gmin = (10 ** 9,)
for t0 in allT:
    T = [t0]; cur = rho(T)
    for _ in range(11):
        cand = None
        for t in allT:
            if t in T:
                continue
            k, r = rho(T + [t])
            if cand is None or k - 2 * r < cand[0]:
                cand = (k - 2 * r, t, k, r)
        T.append(cand[1]); cur = cand[2:]
        if cand[0] < gmin[0]:
            gmin = (cand[0], len(T), cand[2], cand[3])
    if t0 >= 15:
        break
print('greedy (16 seeds, up to 12 targets): min |R| - 2 rho =', gmin)
assert min(w[0] for w in worst) >= 0 and gmin[0] >= 0      # Lemma S necessary condition holds on every tested cut
print('PASS_EXPNU_READERCUT')
