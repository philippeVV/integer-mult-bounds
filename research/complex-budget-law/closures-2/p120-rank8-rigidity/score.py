"""Scores for a unit with copy fraction gamma and transfers per port d at h=9 (labels F2^9, a=8).
rho(d) = max(1/(2a), EFF*d/(2a)) with EFF = 4142/(112*52) (expL2 E8 sharing heuristic; d=52 -> 2.31 today).
(1) pair model (expL2 mult: both sides, best of s3 and cal9; cal9 is rho-blind);
(2) s3-only pair model; (3) expL4 sec.4 Lambda model: (nu/nu_E8)/(0.35+0.65 d/52), nu = 1-9/8 gamma."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rho_min import mult, EFF
def rho(d,h=9): a=h-1; return max(1/(2*a), EFF*d/(2*a))
def lam_model(g,d): nu=1-9/8*g; return (nu/(5/14))/(0.35+0.65*d/52)
if __name__=='__main__':
    print('gamma   d   rho   pair(best)  pair(s3)  Lambda-model')
    for g in (4/7,0.6,0.675,0.75,0.5,0.4):
        for d in (52,26,15,8,4,1):
            print('%.3f %3d  %.3f   x%.2f      x%.2f     x%.2f'%(g,d,rho(d),mult(9,g,rho(d)),mult(9,g,rho(d),False),lam_model(g,d)))
    assert abs(mult(9,64/112,rho(52))-1)<0.05
    print('PASS_SCORE')
