"""Independent exact reproduction of PR2044's outer kappa from new/SPEC.json (no PR code imported),
plus a rigorous mpmath interval enclosure of the native moment at b, b_up, cl and cl+10^-DIG."""
import json, sys, time
from fractions import Fraction as F
sys.set_int_max_str_digits(0)
P = sys.argv[1].rstrip('/') + '/'
spec = json.load(open(P + 'new/SPEC.json'))
oc = json.load(open(P + 'outer/OUTER-CLAIM.json'))
bj = json.load(open(P + 'outer/B.json'))
led = json.loads(open(P + 'outer/LEDGER.json', 'rb').read().replace(b'\r', b''))
DIG = spec['outer_refinement']['endpoint_digits']
cl = F(spec['expected']['complex_lower'])
u = F(1, 10**(DIG + 10))
b = cl - u; q = b - u; eps = 1 - 2*u
kappa = F((eps*q*10**(DIG+46)).__floor__(), 10**(DIG+46)) - F(1, 10**(DIG+22))
out = {}
out['DIG'] = DIG
out['b_matches_B.json'] = b == F(bj['b'])
out['b_up_is_cl_plus_1e-DIG'] = F(bj['b_above_endpoint']) == cl + F(1, 10**DIG)
out['kappa_matches_OUTER-CLAIM_rational'] = kappa == F(oc['kappa_rational'])
row = led['rows']['struct_boundary_radix_union']
out['kappa_matches_LEDGER_rational'] = kappa == F(row['kappa']['rational'])
out['kappa_matches_OUTER-CLAIM_decimal'] = kappa == F(oc['kappa'])
k1897 = None
out['q_minus_kappa'] = '%.6e' % float((q - kappa)*10**30013) + 'e-30013'
out['b_over_1_plus_b_equals_supplier_kappa'] = cl/(1+cl) == F(spec['expected']['kappa']) if 'kappa' in spec['expected'] else None
out['lift_kappa_minus_supplier'] = '%.12e' % float(kappa - F(spec['expected']['kappa']))
out['b2_over_1pb'] = '%.12e' % float(cl*cl/(1+cl))
prof = {int(k): int(v) for k, v in spec['expected']['global_profile'].items()}
m = spec['geometry']['dimension']; W = spec['expected']['live']
out['m_W'] = (m, W); out['M0'] = str(sum(F(n*r, m*W) for r, n in prof.items()))
print(json.dumps(out, indent=1)); sys.stdout.flush()
if len(sys.argv) > 2:
    import mpmath
    from mpmath import iv
    iv.prec = int(sys.argv[2])
    def M(x):
        X = iv.mpf([x.numerator, x.numerator]) / x.denominator
        s = iv.mpf(0)
        for r, n in prof.items():
            s += iv.mpf(n*r) / (m*W) * iv.exp(X * iv.log(iv.mpf(m) / r))
        return s
    for name, x, want in [('b', b, '<1'), ('b_up', F(bj['b_above_endpoint']), '>1'), ('cl', cl, '<1'), ('cl+1e-DIG', cl + F(1, 10**DIG), '>1')]:
        t = time.time(); v = M(x); one = iv.mpf(1)
        ok = (v.b < 1) if want == '<1' else (v.a > 1)
        gap = (v.a - 1) if want == '>1' else (1 - v.b)
        print(name, want, 'PASS' if ok else 'FAIL', 'gap', mpmath.nstr(mpmath.mpf(gap), 6), 'width', mpmath.nstr(mpmath.mpf(v.b - v.a), 3), '%.1fs' % (time.time()-t)); sys.stdout.flush()
