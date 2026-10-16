"""min z subject to l <= c * min w (sufficient form: every hop-shortest path <= c*omega), cutting planes, HiGHS."""
import sys; sys.path.insert(0, __file__.rsplit('/',1)[0])
from variants import full, minimal
from cp import solve
from math import log
a = 147/157
for f in (full, minimal):
    nv, E, ss = f()
    for c in (7, 7.5, 8, 10):
        res, w, y, d, om = solve(nv, E, ss, c=c)
        z = res.fun; lw = d.max()/w.min()
        extra = f"  analytic a+2(1-a)/c = {a+2*(1-a)/c:.10f}" if f is full else ""
        print(f"{f.__name__:8s} c<={c:<4} z*={z:.10f} gap={1/z:.6f} actual l/w={lw:.6f}  allN(ln/4ln c)={log(1/z)/(4*log(lw)):.6f} allN(c1=c+1)={log(1/z)/(4*log(lw+1)):.6f}{extra}", flush=True)
