#!/usr/bin/env python3
"""expQS price.py: what laminar (helper-free) decoders would be worth, in the expL2 pair model (s = 3 via expNE eta3,
s = 9 via expNB cal9 / nohelp) and in a Lambda-share model (a ~ nu/Lambda, helper share of Lambda).
Multipliers are over 0.0017.  Rows:
 (A) realised: no laminar decoder with gamma < 1 was found -> x1.00 (identity decoder: nu < 0).
 (B) partial laminarity on the shipped E8 word: fraction f of corrections direct (maxlam.json), helper mass * (1 - f).
 (C) hypothetical helper-free laminar E8-type unit at gamma (rho = parking floor 1/(2a)).
 (D) small units h_c = 5, 6 at the expN2 LP floors with rho = 1/(2a) (the s = 3 rescue), from expL2's model."""
import sys, json, math
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)   # rho_min (expL2), kappa_pair (expNB), transfer (expNE)
import rho_min as RM
import kappa_pair as KP

REF = 0.0017
def pair(a, b): A = min(a, b); return A / (1 + A) / REF

def e8_row(x_b, x_c, rho_b, rho_c):
    s3 = pair(RM.side('bit', 8, x_b, rho_b, False), RM.side('cx', 9, x_c, rho_c, False))
    s9 = pair(KP.side_saving('bit', 8, x_b, 'cal9'), KP.side_saving('cx', 9, x_c, 'cal9'))
    nh = pair(KP.side_saving('bit', 8, x_b, 'nohelp'), KP.side_saving('cx', 9, x_c, 'nohelp'))
    return s3, s9, nh

if __name__ == '__main__':
    ml = json.load(open(os.path.join(HERE, 'maxlam.json')))
    xb, xc = 63 / 111, 64 / 112; rb, rc = RM.RHO_T['bit'], RM.RHO_T['cx']
    s3, s9, nh = e8_row(xb, xc, rb, rc)
    print('control: E8 today  s=3 x%.3f  s=9 cal9 x%.3f (expNB: 1.009)  s=9 nohelp x%.3f' % (s3, s9, nh))
    assert abs(s9 - 1.0088) < 0.003
    print('\n(B) partial laminarity on E8 (greedy direct fraction f of corrections; helper mass scaled by 1-f)')
    for name in ('P120', 'P112', 'P104'):
        f = ml[name]['frac']
        s3f = pair(RM.side('bit', 8, xb, rb * (1 - f), False), RM.side('cx', 9, xc, rc * (1 - f), False))
        # Lambda-share model at s = 9: a ~ 1/Lambda, helper share phi of Lambda; phi from the proved x1.88 cap and from #1094's 65%
        lam = {phi: 1 / (1 - phi * f) for phi in (1 - 1 / 1.88, 0.65)}
        print('  %s f=%.3f: s=3 x%.3f (was %.3f); s=9 Lambda-share x%.3f (phi=.468) / x%.3f (phi=.65) over today'
              % (name, f, s3f, s3, s9 * lam[1 - 1 / 1.88], s9 * lam[0.65]))
    print('\n(C) hypothetical helper-free laminar E8-type unit (bit h=8, complex h=9, x = gamma both sides, rho = 1/16)')
    print('  gamma    s=3(rho=1/16)  s=9 nohelp  best')
    rows = []
    for g in (64 / 112, 0.60, 0.65, 0.70, 0.80, 8 / 9):
        s3g = pair(RM.side('bit', 8, g, 1 / 14, False), RM.side('cx', 9, g, 1 / 16, False))
        nhg = pair(KP.side_saving('bit', 8, g, 'nohelp'), KP.side_saving('cx', 9, g, 'nohelp'))
        rows.append((g, s3g, nhg)); print('  %.3f    x%-12.2f  x%-9.2f  x%.2f' % (g, s3g, nhg, max(s3g, nhg)))
    print('\n(D) small units at the expN2 LP floors, rho = 1/(2a) (s = 3 rescue), from expL2 mult()')
    for h, x in ((5, 0.5485), (6, 0.4862), (5, 0.45), (6, 0.40)):
        print('  h_c=%d x=%.4f: rho=1/(2a) x%.2f | rho=today x%.2f' % (h, x, RM.mult(h, x, 1 / (2 * (h - 1))), RM.mult(h, x, RM.RHO_T['cx'])))
    json.dump(dict(control=(s3, s9, nh), C=rows), open(os.path.join(HERE, 'price.json'), 'w'), indent=1)
    print('PASS_EXPQS_PRICE')
