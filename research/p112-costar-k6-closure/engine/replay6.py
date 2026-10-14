import pickle, sys, glob, os
from multiprocessing import Pool
import verify_cert as V
def run(path):
    out = path + ".ok"
    if os.path.exists(out): return open(out).read()
    jobs = pickle.load(open(path, "rb"))
    for A, c, fin, steps in jobs:
        V.replay(A, c, fin, steps)
    open(out, "w").write(str(len(jobs)))
    return str(len(jobs))
if __name__ == "__main__":
    files = sorted(glob.glob("k6chunks/c*.pkl"))
    with Pool(4) as p:
        tot = sum(int(x) for x in p.imap_unordered(run, files))
    print("k=6 certificates replayed:", tot)
