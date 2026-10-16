"""Per-orbit idempotent-rank caps rho'_K.  ub: exact upper bounds (rank_Q default; rpub.py scans).  lb: lower bounds
(feasible idempotent supports found by rpgrow/lball; the true rho' lies in [lb, ub])."""
import pickle
UB = {6: {26: 8, 32: 8, 16: 4, 18: 4, 20: 4}, 5: {7: 4, 43: 15}}          # rpub.py: exact, no idempotent support of larger rank
LBX = {6: {26: 8, 32: 8, 77: 16}, 5: {}}                    # rpgrow runs (best found)
def caps(n, mode, NK):
    import core6
    from rhoprime import rkQ
    D = core6.load(n); korb = sorted(set(D['corb']))
    out = []
    for k, ko in enumerate(korb):
        i = D['corb'].index(ko); Wm, Qm = D['cols'][i]
        r = rkQ(core6.bits(Wm), core6.bits(Qm)) if D['beta'][i] else 0
        if mode == 'ub': out.append(min(r, UB[n].get(k, r)))
        else:
            lb = pickle.load(open('lb%d.pkl' % n, 'rb'))
            out.append(max(lb.get(k, 0), LBX[n].get(k, 0)))
    return out
