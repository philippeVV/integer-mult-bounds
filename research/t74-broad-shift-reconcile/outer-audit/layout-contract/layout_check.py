"""Exact profile-dependent quantities of the compact-layout contract for PR2044's T74 supplier (and the selected T16 for comparison).
Usage: python3 -I layout_check.py SPEC.json T16_PROFILE.json OUT.json  (stdlib + mpmath)"""
import json, sys, math
from fractions import Fraction as Q
import mpmath as mp
sys.set_int_max_str_digits(0)
spec=json.load(open(sys.argv[1])); t16=json.load(open(sys.argv[2]))
e=spec['expected']; m=spec['geometry']['dimension']; W=int(e['live']); phys=int(e['physical'])
H={int(r):int(n) for r,n in e['global_profile'].items()}
H16={int(r):int(n) for r,n in t16['profile'].items()}; W16,phys16=14336,14826
out={}
def prof(H,W,phys,tag):
    s=sum(H.values()); mass=sum(r*n for r,n in H.items()); tmax=max(H)
    hd=1
    while m**hd<=2*tmax**hd: hd+=1
    group=2**(36*36)
    for j in range(1,37): group*=2**(2*j)-1
    stock=group*phys
    mp.mp.dps=40
    d=dict(children_s=s, rank_mass=mass, deficit=m*W-mass, s_over_W=str(mp.mpf(s)/W),
           maxchild=tmax, halving_degree=hd, halving_ints=[m**hd,2*tmax**hd],
           physical_minus_live=phys-W, stock_bits=stock.bit_length(),
           naive_uniform_sigma_log_m_s_over_W=str(mp.log(mp.mpf(s)/W)/mp.log(m)),
           rank_one_expansion_rate_1_minus_log_m_mass_over_W=str(1-mp.log(mp.mpf(mass)/W)/mp.log(m)))
    out[tag]=d; return d,stock
T,stock74=prof(H,W,phys,'T74_pr2044'); prof(H16,W16,phys16,'T16_pr1094')
# depth function of endpoint-gauge-complex.tex: D(e)=0 (e<m), 1+D(t*floor(e/m))
def D(x,t=48):
    k=0
    while x>=m: x=t*(x//m); k+=1
    return k
dep={}
for L in [10,100,1000,10**4,10**5]:
    h=2**L; dep[str(L)]=dict(D=D(h),D_over_log2h=str(Q(D(h),L)),two_log2h=2*L)
out['depth_D_of_2^L']=dep
out['depth_asymptotic_slope_1_over_log2(73/48)']=str(1/mp.log(mp.mpf(73)/48,2))
# row stock: chain's coefficient (2*stock_bits + 9909 + 252 + 2*payload)
payload=spec['geometry']['charged_prefix_bits']
row=2*stock74.bit_length()+9909+252+2*payload
out['row_coefficient_T74']=row; out['row_limit']=20161
# chunks q0 = ceil(coef*log2 h / K), h<=p, K=128(4ceil(log2 p)+6)
def q0(coef,X): return -(-coef*X//(128*(4*X+6)))
out['q0_chunks_at_h=p']={str(X):dict(row=q0(row,X),conservative_51_25_x_20161=q0(-(-51*20161//25),X),degree_1e6=q0(10**6,X)) for X in [64,10**4,10**30]}
def thresh(c):
    h=1
    while True:
        if all(-(-2*x//128)-(-x//128)+c<=Q(x,16) for x in range(h,h+256)):
            # verify monotone afterwards: for x>=h, 3x/128+2+c<=x/16 iff x>=(2+c)*128/5
            assert Q(3*h,128)+2+c<=Q(h,16) or h>=(2+c)*128//5-256
            return h
        h+=1
out['R(h)<=h/16_threshold']={str(c):max(thresh(c),-(-(2+c)*128//5)) for c in [1,39,40,79,1954]}
# stop fraction
b=Q(1)-0; certified=Q(e['complex_lower']); DIG=spec['outer_refinement']['endpoint_digits']
R0=Q(1,10**(DIG+10)); b=certified-R0; q=b-R0; a=(q+b)/2; STOP=Q(1,10**(DIG+16))
def sci(x): mp.mp.dps=20; return mp.nstr(mp.mpf(x.numerator)/x.denominator,12)
out['stop']={'STOP':'1e-%d'%(DIG+16),'leaf_minus_a':sci((1-STOP)*b-a),'leaf_minus_q':sci((1-STOP)*b-q),'with_stop_1e-9_leaf_minus_a':sci((1-Q(1,10**9))*b-a),'onset_log2_h_for_h^stop>=2m':mp.nstr(mp.log(2*m,2),8)+'e%d'%(DIG+16)}
out['tau_minus_sigma(=b-a)']=sci(b-a)
json.dump(out,open(sys.argv[3],'w'),indent=1); print(json.dumps(out,indent=1))
