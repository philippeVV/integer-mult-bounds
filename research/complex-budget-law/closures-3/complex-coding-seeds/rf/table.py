"""Digest of run.py logs: per layout, best S0/E (char 0), best S2/E (saturated mod-2), count of rows with S2 < S0."""
import json, glob, math, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))   # logs live next to this script
OK = True
print('| layout | n | E | girth | rows | best S0 (S0/E) | config | S2 of that config | best S2/E | rows S2<S0 | rows S2>S0 |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for f in ['pg3', 'pgm4', 'biaffine4', 'w3', 'wm4', 'pg7', 'pgm8', 'biaffine8', 'w7', 'wm8', 'biaffine16', 'pgm16']:
    try: L = [json.loads(l) for l in open(f + '.log')]
    except FileNotFoundError: continue
    h = L[0]; rows = [r for r in L if 'S0' in r and 'TP' in r]
    if not rows: continue
    b = max(rows, key=lambda r: (r['S0'], r['S2']))
    OK &= b['S0'] / h['E'] < 1 / 6 and all(r['S2'] >= r['S0'] for r in rows)
    print(f"| {h['layout']} | {h['n']} | {h['E']} | {h['girth']} | {len(rows)} | {b['S0']} ({b['S0']/h['E']:.4f}) | {b['TP']}/{b['TL']}/{b['mode']} | {b['S2']} | "
          f"{max(r['S2'] for r in rows)/h['E']:.4f} | {sum(r['S2'] < r['S0'] for r in rows)} | {sum(r['S2'] > r['S0'] for r in rows)} |")
print('PASS_EXPRF_TABLE' if OK else 'FAIL_EXPRF_TABLE')   # every S0/E < 1/6 and S2 >= S0 in every row
