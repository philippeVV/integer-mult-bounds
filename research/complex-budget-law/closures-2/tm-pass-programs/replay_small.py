#!/usr/bin/env python3
"""Replay the small exact cells of the TM-native tables (step model and pass model) and check route == xor.

Step model (tmsat.py): exact minimum steps of a k-tape oblivious binary TM with c control bits that transposes a 2x2 matrix
of b-bit records, searched upward from T = 1 (so the minimum is certified by UNSAT at T - 1), workspace L = V + 2.
Pass model (passsat.py): exact minimum whole-volume masked-lag passes interchanging two u-bit fields.
Expected values are the scouting tables (README).  Needs python-sat (CaDiCaL), e.g.
    uv run -q --no-project --python 3.12 --with python-sat python -I -B replay_small.py
Prints PASS_TM_REPLAY.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tmsat import min_T, rec_transpose   # noqa: E402
from passsat import minP                 # noqa: E402

STEP = {  # (b, k, c): expected minT for both route and xor
    (1, 1, 1): 4, (1, 2, 0): 4, (1, 2, 1): 4,
    (2, 1, 1): 11, (2, 2, 0): 9, (2, 2, 1): 9,
    (3, 2, 0): 15, (3, 2, 1): 14, (3, 3, 0): 14,
}
PASS = {(1, 2): 3, (1, 3): 3, (2, 2): 6}   # (u, k): expected minP for both route and xor

ok = True
print('step model, 2x2 transpose of b-bit records (L = V + 2, searched from T = 1):')
for (b, k, c), want in STEP.items():
    perm = rec_transpose(b, 2); V = len(perm)
    got = {}
    for mode in ('route', 'xor'):
        t0 = time.time()
        got[mode] = min_T(V, perm, k, c, V + 2, mode, Tlo=1, log=False)
        got[mode + '_s'] = time.time() - t0
    good = got['route'] == got['xor'] == want
    ok &= good
    print(f'  b={b} V={V:2d} k={k} c={c}: route {got["route"]:3d}  xor {got["xor"]:3d}  expected {want:3d}  '
          f'{"OK" if good else "MISMATCH"}  ({got["route_s"]:.1f}s + {got["xor_s"]:.1f}s)', flush=True)
print('pass model, 2^u x 2^u field interchange:')
for (u, k), want in PASS.items():
    got = {}
    for mode in ('route', 'xor'):
        t0 = time.time()
        got[mode] = minP(u, k, mode, log=False)
        got[mode + '_s'] = time.time() - t0
    good = got['route'] == got['xor'] == want
    ok &= good
    print(f'  u={u} k={k}: route {got["route"]}  xor {got["xor"]}  expected {want}  {"OK" if good else "MISMATCH"}  '
          f'({got["route_s"]:.1f}s + {got["xor_s"]:.1f}s)', flush=True)
print('PASS_TM_REPLAY' if ok else 'FAIL_TM_REPLAY')
