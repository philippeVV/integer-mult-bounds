"""expRF markov.py method at unequal traces: free law of PQP, tr P = a = 8/17 (points), tr Q = b = 9/17 (lines).
nu = (1-min(a,b)) delta_0 + max(a+b-1,0) delta_1 + density sqrt((x-l-)(l+-x))/(2 pi x(1-x)) on [l-,l+],
l+- = a+b-2ab +- sqrt(4ab(1-a)(1-b)).  Girth 2g=6: moments j<=g-1=2 equal free ones (tree-like), so
S0/E = nu_true({1}) <= M = min{ int p dnu_free : deg p <= 2, p >= 0 on [0,1], p(1)=1 }.  Also checks half-rate -> 2/11."""
import numpy as np
from math import log, sqrt
from scipy.optimize import linprog
from scipy.integrate import quad
def law(a,b):
    lm=a+b-2*a*b-sqrt(4*a*b*(1-a)*(1-b)); lp=a+b-2*a*b+sqrt(4*a*b*(1-a)*(1-b))
    dens=lambda x: sqrt(max(0,(x-lm)*(lp-x)))/(2*np.pi*x*(1-x))
    at0=1-min(a,b); at1=max(a+b-1,0)
    def mom(k):
        v=quad(lambda x: x**k*dens(x),lm,lp,limit=200)[0]
        return v+at1+(at0 if k==0 else 0)
    return lm,lp,at0,at1,mom
def M(a,b,g=3,grid=20001):
    lm,lp,at0,at1,mom=law(a,b); d=g-1
    mo=np.array([mom(k) for k in range(d+1)])
    xs=np.linspace(0,1,grid); Vm=np.vstack([xs**k for k in range(d+1)]).T
    r=linprog(mo,A_ub=-Vm,b_ub=np.zeros(len(xs)),A_eq=np.ones((1,d+1)),b_eq=[1.0],bounds=[(None,None)]*(d+1),method='highs')
    return r.fun,mo,(lm,lp,at0),r.x
thr=lambda a,b=5: 1/(1+b**(1-a))
for (a,b),name in [((.5,.5),'half-rate check (expRF 2/11)'),((8/17,9/17),'PG(2,16) full [17,8]/[17,9]'),((8/16,8/16),'PG(2,16)-M [16,8]')]:
    m,mo,sup,p=M(a,b)
    amax=1-log((1-m)/m)/log(5)
    print(f'{name}: traces ({a:.4f},{b:.4f}) moments {np.round(mo,5)} support {np.round(sup,4)} M={m:.5f} '
          f'thr(a=0)={thr(0):.4f} thr(a=1/9)={thr(1/9):.4f} a_max={amax:+.4f} kappa_max={amax/(1+amax):.4f} p={np.round(p,4)}')
# limits-2 check: the half-rate value is expRF's 2/11 and the full PG(2,16) cap is below the a = 1/9 threshold
mh=M(.5,.5)[0]; mf=M(8/17,9/17)[0]; af=1-log((1-mf)/mf)/log(5)
ok=abs(mh-2/11)<1e-4 and abs(mf-0.18147)<5e-5 and mf<thr(1/9) and abs(af/(1+af)-0.0602)<5e-4
print('PASS_EXPXX_MARKOV' if ok else 'FAIL_EXPXX_MARKOV')
