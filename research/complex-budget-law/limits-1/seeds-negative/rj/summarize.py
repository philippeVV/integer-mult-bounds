#!/usr/bin/env python3
"""expRJ: tables from big_*.jsonl (a_H), rr.jsonl (exact LP), tanner.jsonl. stdlib only."""
import json, math, os, sys
def load(f):
    return [json.loads(l) for l in open(f)] if os.path.exists(f) else []

print("== A. far-pair, uniform-hop ceiling a_H = 1 - ln(E/(r k))/ln b_min, best b0 per graph (r = min(sigma_hat, 16-d))")
big = load('big_rr.jsonl') + load('big_struct.jsonl')
best = {}
for r in big:
    key = (r['name'], r['n'])
    if key not in best or r['a_H'] > best[key]['a_H']: best[key] = r
print(f"{'seed':22s} {'n':>7s} {'d':>3s} {'diam':>4s} {'b0':>3s} {'k':>6s} {'sigma':>6s} {'r':>6s} {'roles':>6s} {'arity':>5s} {'E/rk':>6s} {'a_H':>7s}")
first = {}
for (nm, n), r in sorted(best.items(), key=lambda kv: (kv[1]['d'], kv[0][0], kv[0][1])):
    print(f"{nm:22s} {n:7d} {r['d']:3d} {r['diam']:4d} {r['b0']:3d} {r['k']:6d} {r['sigma']:6.3f} {r['rate']:6.3f} {r['roles']:6d} {r['arity']:5.2f} {r['E_over_rk']:6.3f} {r['a_H']:7.4f}")
    for th in (0.317, 0.5):
        if r['a_H'] >= th and r['roles'] <= 1e5 and r['arity'] <= 16:
            if th not in first or r['roles'] < first[th]['roles']: first[th] = r
for th, r in sorted(first.items()):
    print(f"smallest (by roles) with a_H >= {th}: {r['name']} n={r['n']} b0={r['b0']} k={r['k']} roles={r['roles']} arity={r['arity']} a_H={r['a_H']}")

print("\n== B. far-pair, exact HiGHS max-concurrent-flow (<=512 roles): gap_cap = sigma_hat/lambda*, a_G = ln(gap)/ln b0, a_H")
rr = load('rr.jsonl')
print(f"{'seed':10s} {'n':>5s} {'k':>4s} {'b0':>3s} {'bmean':>6s} {'sigma':>6s} {'lam*':>7s} {'lam_vol':>7s} {'lamECMP':>7s} {'gap':>6s} {'a_G':>7s} {'a_H':>7s} {'arity':>5s} {'tLP':>6s}")
for r in rr:
    aH = math.log(r['rate'] * r['k'] * r['b0'] / r['E']) / math.log(r['b0'])
    print(f"{r['name']:10s} {r['n']:5d} {r['k']:4d} {r['b0']:3d} {r['bmean']:6.2f} {r['sigma']:6.3f} {r['lam_lp']:7.4f} {r['lam_vol']:7.4f} {r['lam_ecmp']:7.4f} {r['gap']:6.3f} {r['a_cap']:7.4f} {aH:7.4f} {r.get('arity', 0):5.2f} {r['tlp']:6.1f}")

print("\n== C. Tanner hybrid (sessions = balanced edge set R, rate <= 1, arity = degree), exact LP")
for r in load('tanner.jsonl'):
    print(f"{r['name']:22s} n={r['n']:4d} S={r['k']:4d} S/E={r['S_over_E']:.3f} b0={r['b0']} sigma={r['sigma']:.3f} lam*={r['lam_lp']:.4f} "
          f"gap={r['gap']:.3f} a_G={r['a_cap']:.4f} a_expQK(uniform)={r['a_expQK']:.4f}")
print("PASS_EXPRJ_SUMMARY")
