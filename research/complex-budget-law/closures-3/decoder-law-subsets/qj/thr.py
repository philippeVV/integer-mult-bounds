import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import e8m, sys
rows = [list(map(int, l.split())) for l in open(sys.argv[1] if len(sys.argv) > 1 else 'sub_out.txt')]
out = {}
for (m, s), b in zip(e8m.CLASSES, rows):
    n = bin(s).count('1')
    T = tuple(n - x for x in b)   # deletions for span <= 8,7,...,1
    out.setdefault(n, set()).add(T)
for n in sorted(out):
    for T in sorted(out[n]): print(n, 120 - n, T)
