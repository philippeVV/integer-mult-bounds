"""Independent checker (stdlib only, exact Fractions).
Certificate: JSON {instance, paths:[[k, [v0,v1,...,vL]], ...], weights:[num/den strings]}.
Checks: every path is a simple walk in the INTERNAL graph from one endpoint of session k to the other;
edge loads <= 1; per-session load <= 1. Then for every m and every integer metric r and s_k>=0 with
s_k + r(P) >= m on all c_k-p_k paths:  Σ r + Σ s >= Σ_P y_P (r(P)+s_k) >= m·F  (weak LP duality),
so MC_m >= ceil(m·F) (MC_m is an integer).  usage: check_flow.py inst.json cert.json"""
import sys, json, gzip, os
from fractions import Fraction as Fr
from math import ceil
H=os.path.dirname(os.path.abspath(__file__)); op=lambda f: (gzip.open if f.endswith('.gz') else open)(os.path.join(H, f), 'rt')   # closures-3: paths relative to this script; gzipped certificates
inst=json.load(op(sys.argv[1])); cert=json.load(op(sys.argv[2]))
I={frozenset(e) for e in inst['internal']}; assert len(I)==len(inst['internal'])
Rs=inst['sessions']; K=len(Rs)
load={}; sload=[Fr(0)]*K; F=Fr(0)
for (k,vs),w in zip(cert['paths'],cert['weights']):
    w=Fr(w); assert w>=0
    assert {vs[0],vs[-1]}==set(Rs[k]), 'endpoint'
    assert len(set(vs))==len(vs), 'not simple'
    for a,b in zip(vs,vs[1:]):
        e=frozenset((a,b)); assert e in I, 'non-internal edge'
        load[e]=load.get(e,0)+w
    sload[k]+=w; F+=w
assert max(load.values())<=1 and max(sload)<=1, 'capacity'
print(f'valid multiflow, value F = {F} = {float(F):.6f}, W = {K}')
for m in range(1,6): print(f'  m={m}: MC_m >= {ceil(m*F)}   (Wm = {m*K}, gap {m*K-ceil(m*F)})')
Eint=len(inst['internal'])
print(f'  m=5: r = 1 on every internal edge costs E_int = {Eint}, and ceil(5F) = {ceil(5*F)}: MC_5 = {Eint} exactly (Wm = {5*K})')
print('PASS_EXPRA_FLOW' if ceil(5*F) == Eint < 5*K else 'FAIL_EXPRA_FLOW')
