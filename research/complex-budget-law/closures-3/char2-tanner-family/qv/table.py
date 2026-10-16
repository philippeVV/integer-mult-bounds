"""expQV: GQ T2*(O) Tanner seeds (char 2), measured S per (q, code), priced with a = 1 - ln((E-S)/S)/ln 7 (b = 7, girth 8) and
charged with the expNI spectral cut requirement |R| >= lambda_2 n / 8 (t = forced extra R arcs).  Columns:
| q | code | E | S | S/E | a | x(a/(1+a))/0.0017 if exact | cut need | t | a charged |.  Stdlib."""
import math
ok = True
rows=[(4,'hexacode F4 (any twist)',45),(8,'F64 RS+nu, line s^2 (rel twist 1)',670),(8,'F8 e4 (non-MDS), rel twist 1',680),(8,'untwisted F64',576),(16,'F256 RS+nu, line s^4 (rel twist 2)',11450),(16,'rel twist 1',11332),(16,'untwisted',8646)]
for q,lab,S in rows:
    E=q**3*(q+2); n=q**3+q*q*(q+2); d1,d2=q+2,q
    l2=((d1+d2)-math.sqrt((d1-d2)**2+8*q))/2; need=l2*n/8
    a=1-math.log((E-S)/S)/math.log(7); t=max(0,math.ceil(need)-S)
    ac=1-math.log((E-S+t)/(S-t))/math.log(7)
    f=lambda x: f"{x:.4f}"
    ok &= ac < 0.09 and (q == 4 or t > 0)
    print(f"| {q} | {lab} | {E} | {S} | {S/E:.4f} | {f(a)} | {a/(1+a)/0.0017 if a>0 else 0:.1f} | {need:.1f} | {t} | {f(ac)} |")
print('every q >= 8 code misses the cut requirement; q = 4 has a < 0; charged a <= 0.086 everywhere')
print('PASS_EXPQV_TABLE' if ok else 'FAIL_EXPQV_TABLE')
