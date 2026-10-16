"""expXH price2: tensor law h_C = u*h_R, gamma_C = gamma of the real pattern under a complex decoder.
Price lp (expNM kappa.lp, data floor) and the selected-word formula (expXA price.py sel) at (u*h_R, u, gamma) vs (h_R, 1, gamma)."""
import sys, os, io, contextlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'xa'))
with contextlib.redirect_stdout(io.StringIO()):
    import price as P
ref_lp = P.lp(9, 1, 4 / 7)[0]; ref_sel = P.sel(9, 4 / 7)
print(f"E8 P112 control: lp {ref_lp:.5f} (x{ref_lp/0.0017:.2f} over 0.0017), sel {ref_sel:.5f} (x{ref_sel/0.0017:.2f})")
print(f"{'h_R':>3s} {'u':>2s} {'h':>3s} {'k':>6s} {'gamma':>6s} | lp (x over 0.0017, x vs E8) | sel (x over 0.0017, x vs E8)")
for hR in (5, 7, 9, 10, 13):
    for g in (4 / 7, 0.4, 0.2, 0.0):
        for u in (1, 2, 3, 4, 8):
            h = u * hR
            kk, s = P.lp(h, u, g)
            se = P.sel(h, g, u=u)
            print(f"{hR:3d} {u:2d} {h:3d} {u/(h-1):6.3f} {g:6.3f} | {kk:.5f} (x{kk/0.0017:5.2f}, x{kk/ref_lp:4.2f}) s{s} | {se:.5f} (x{se/0.0017:5.2f}, x{se/ref_sel:4.2f})")
# limits-1 replay: under the tensor law u > 1 never beats u = 1 at the same (h_R, gamma), in either regime
for hR in (5, 7, 9, 10, 13):
    for g in (4 / 7, 0.4, 0.2, 0.0):
        lp1, sel1 = P.lp(hR, 1, g)[0], P.sel(hR, g, u=1)
        for u in (2, 3, 4, 8):
            assert P.lp(u * hR, u, g)[0] <= lp1 + 1e-12 and P.sel(u * hR, g, u=u) <= sel1 + 1e-12, (hR, g, u)
print('PASS_EXPXH_PRICE2')
