import sys; sys.argv = ['inc.py', '0', '1']
exec(open('inc.py').read().split('p = rng.normal(size=2*nv); m = 0')[0])
Lall = np.exp(-2j*np.pi*np.outer(np.arange(N), D)/N)
for s in range(3):
    p = np.random.default_rng(s).normal(size=2*nv); J = jac(p, Lall)
    sv = np.linalg.svd(J, compute_uv=False); print('real rank', (sv > sv[0]*1e-10).sum(), 'of', J.shape, 'sv tail', sv[95*2-4:95*2+4] if len(sv)>200 else '')
    print(np.round(np.log10(sv/sv[0]),1)[::10])
