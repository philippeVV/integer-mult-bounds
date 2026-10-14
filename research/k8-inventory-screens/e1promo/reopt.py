"""Exact-domain frame re-optimisation with optional E1 promotions (expB discovery tool, needs scipy).

Uses PR1148's frame model (research/k8-eight-ends/discovery/frame_model.py, downloaded read-only to pr1148/ee/) for the
complete frame domains [forward join, backward meet] and the per-component local-polytope LP; MILP fallback if fractional.
Objective: caller's per-child excess phi(r)-r (EXC in e1lib), i.e. the same ranking as dRC at fixed mass.
usage: reopt.py WORD.json.gz OUT.json.gz [--promote H:x ...] [--only-changed]
  --promote H:x  raise helper H's start frame by vector x (E' = E + <x>)
"""
import sys, os, json, gzip, hashlib, importlib.util, argparse, time
from collections import Counter
import numpy as np
from scipy.optimize import linprog, milp, Bounds, LinearConstraint
from scipy.sparse import csc_matrix
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import e1lib as L
s = importlib.util.spec_from_file_location('frame_model', os.environ.get('FRAME_MODEL', os.path.join(HERE, 'pr1148/ee/frame_model.py')))
M = importlib.util.module_from_spec(s); s.loader.exec_module(M)
COST = [0.0] + [L.EXC[r] for r in range(1, 10)]


def solve_component(labels, unary, pairs):
    rr, rhs, cols = M.lp(labels, unary, pairs)
    ri, ci, val = [], [], []
    for j, (_, _, es) in enumerate(cols):
        for i, a in es: ri.append(i); ci.append(j); val.append(a)
    A = csc_matrix((np.array(val, float), (ri, ci)), shape=(len(rr), len(cols)))
    c = np.array([sum(n * COST[d] for d, n in enumerate(vec)) for _, vec, _ in cols])
    sol = linprog(c, A_eq=A, b_eq=rhs, bounds=(0, None), method='highs'); assert sol.success, sol.message
    frac = any(1e-7 < x < 1 - 1e-7 for x in sol.x)
    if frac:
        sol = milp(c, integrality=np.ones(len(cols)), bounds=Bounds(0, 1), constraints=LinearConstraint(A, rhs, rhs),
                   options={'time_limit': 60}); assert sol.x is not None, sol.message
    choose = {key[1]: key[2] for j, (key, _, _) in enumerate(cols) if key[0] == 'x' and sol.x[j] > .9}
    return choose, frac, float(sol.fun)


def optimize(nw, comps_filter=None):
    m = M.Model(nw, relax_given=True); new = m.baseline(); nfrac = 0; ncomp = 0
    given = {}
    for comp in m.components():
        labels, unary, pairs = m.component_problem(comp)
        choose, frac, val = solve_component(labels, unary, pairs); nfrac += frac; ncomp += 1
        assert len(choose) == len(comp)
        new.update(choose)
    for i in range(len(m.gates)):
        fr = m.frame_of(('G', i), new)
        if i < m.nA: nw['A'][i][-1] = list(fr)
        elif i < m.nA + m.nB: nw['B'][i - m.nA][-1] = list(fr)
        else: nw['early_cleanup'][i - m.nA - m.nB][1] = list(fr)
    prof = m.profile(new)
    return nw, prof, ncomp, nfrac


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('word'); ap.add_argument('out'); ap.add_argument('--promote', nargs='*', default=[])
    a = ap.parse_args(); t = time.time()
    w, sha = L.load(a.word)
    for p in a.promote:
        h, x = map(int, p.split(':'))
        old = L.basis(w['start_extra'][str(h)]); new = L.basis(list(old) + [x]); assert len(new) == len(old) + 1
        w['start_extra'][str(h)] = list(new)
    nw, prof, ncomp, nfrac = optimize(w)
    # recompute report profile/widths exactly with the audit-equivalent forward pass
    W = L.Word(nw); paid = W.prepare() + W.centre_paid()
    assert dict(paid) == {int(k): c for k, c in prof.items()}, (dict(paid), prof)
    v, R = nw['report']['v'], nw['report']['R']
    inv = Counter(L.dim(W.end[q]) - L.dim(W.start[q]) for q in range(2 * v, 2 * v + R))
    nw['report']['local_profile'] = {str(k): c for k, c in sorted(paid.items()) if c}
    nw['report']['local_rank_mass'] = sum(k * c for k, c in paid.items())
    nw['report']['width_inventory'] = {str(k): c for k, c in sorted(inv.items()) if c}
    data = json.dumps(nw, sort_keys=True, separators=(',', ':')).encode()
    open(a.out, 'wb').write(gzip.compress(data, mtime=0))
    print(json.dumps(dict(input=sha, promote=a.promote, profile=nw['report']['local_profile'], widths=nw['report']['width_inventory'],
                          components=ncomp, fractional=nfrac, out_sha=hashlib.sha256(data).hexdigest(), seconds=round(time.time() - t, 1))), flush=True)
