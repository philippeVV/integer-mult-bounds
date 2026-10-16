"""expQL pricing: hypothetical gamma = 1/phi_eff priced with expNM's kappa.lp (framework-ceiling engine, best odd s).
The E8-40 rows are hypothetical realisations (none exists in PG(3,4) or PG(3,8)); the control is today's E8 P112."""
import os, sys, io, contextlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
with contextlib.redirect_stdout(io.StringIO()):
    from kappa import lp            # closures-2/field-line-large-ports/kappa.py, unchanged
res = {}
for nm, h, u, phi in [('E8 P112 control', 8, 1, 1.75), ('E8 P112 control h9', 9, 1, 1.75), ('F4-18 in PG(3,8) [hyp]', 12, 3, 1.125),
                      ('E8-40 in PG(3,8) [hyp]', 12, 3, 1.25), ('E8-40 in PG(3,16) [hyp]', 16, 4, 1.25), ('phi=1.3 h12', 12, 3, 1.3),
                      ('phi=1.45 h12 (break-even)', 12, 3, 1.45), ('phi=1.75 h12', 12, 3, 1.75), ('phi=1.75 h16', 16, 4, 1.75)]:
    kk, s = lp(h, u, 1 / phi); res[nm] = kk
    print(f"{nm:28s} h={h:2d} u={u} k={u/(h-1):.3f} gamma={1/phi:.3f}: " + (f"ceiling kappa {kk:.6f} = x{kk/0.0017:.2f} over 0.0017 (s{s})" if s else "nu<=0: no gain"))
ctl = min(res['E8 P112 control'], res['E8 P112 control h9'])
assert res['E8-40 in PG(3,8) [hyp]'] < ctl and res['E8-40 in PG(3,16) [hyp]'] < ctl
print('E8-40 (phi = 5/4), even if realised, prices below the E8 control: idea 1 realised x1.0')
print('PASS_EXPQL_PRICE')
