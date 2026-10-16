"""Heuristic pricing of the certified sparse patterns (expQF price.py): gamma ~ 4/7 + 0.075 (K - 8) (expNO extrapolation).
Run: python3 -I -B price_rank.py   (stdlib; uses score.py in this directory)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score import rho, mult, lam_model
# hypothetical: structured sparse pattern with certified min-rank K>=12..16; gamma per expNO extrapolation 4/7 + 0.075 (K-8)
for K, d in ((12, 30), (15, 28), (16, 17), (16, 26)):
    g = 4 / 7 + 0.075 * (K - 8)
    print('K=%d d=%d gamma~%.3f  pair x%.2f  Lambda x%.2f' % (K, d, g, mult(9, g, rho(d)), lam_model(g, d)))
    assert mult(9, g, rho(d)) < 1 and lam_model(g, d) < 1
print('PASS_PRICE')
