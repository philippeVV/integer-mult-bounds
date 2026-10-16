import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # -I drops the script dir
import e8m, sys
sets = [e8m.supp(e8m.LINES[i]) for i in range(120)]
with open(sys.argv[1], 'w') as f:
    f.write("%d\n%s\n%d\n" % (120, " ".join(map(str, e8m.LAB)), len(sets)))
    for s in sets: f.write("%d %d\n" % (s >> 64, s & ((1 << 64) - 1)))
