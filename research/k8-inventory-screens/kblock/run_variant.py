"""Rebuild the k8 leaver word up to the cold44 input (HPE word) for a declared variant, without pins, and report the
exact per-block pairing (|I|, rank C, dim L, dim K, rank(LK)) with pairlib (same definitions as cold44x.py).
Steps are reproduce_leaver.py's 1-4 (gauge-forward, port0-3, regroup where the pinned `before` gates occur, greedy HPE);
the repository files are used read-only, outputs go to --out.
usage: python3 run_variant.py --main REPO --pr635 PR635 --word-dir NOTRI_BASE --out DIR
                              [--xparam X.json] [--cells CODED-CELLS.json] [--cgx CGX.py]
"""
import argparse, gzip, hashlib, importlib.util, json, os, shutil, subprocess, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pairlib import Word, rd


def enc(x): return json.dumps(x, sort_keys=True, separators=(',', ':')).encode()


def main():
    ap = argparse.ArgumentParser()
    for k in ('--main', '--pr635', '--word-dir', '--out'): ap.add_argument(k, type=Path, required=True)
    ap.add_argument('--xparam', type=Path); ap.add_argument('--cells', type=Path); ap.add_argument('--cgx', type=Path)
    a = ap.parse_args()
    L = a.main.resolve() / 'research/gauge44-leaver8/construction'; HERE = L / 'leaver'; NOTRI = L / 'notri'
    out = a.out.resolve(); assert not out.exists(); out.mkdir(parents=True)
    up = a.pr635.resolve() / 'research/complex-571-cold14-st2'
    src = a.main.resolve() / 'research/complex-k492-owner-repair/vendor/multicold'
    base = src / 'vendor/pr558'; parent = base / 'vendor/pr552'; gk = parent / 'complex-local/static44/gk'
    structural = a.main.resolve() / 'research/structural-gain/frontier-composition'
    wd = a.word_dir.resolve()
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONUTF8='1', PYTHONINTMAXSTRDIGITS='0', MULTICOLD_SOURCE=str(src), GAIN_PR564_SOURCE=str(src))
    for k in ('PYTHONOPTIMIZE', 'XPARAM', 'EXTRAH', 'PAIRING_RANKS', 'DIAG', 'HYBBUDGET'): env.pop(k, None)

    def run(label, script, *args):
        with (out / (label + '.log')).open('w') as f:
            rc = subprocess.call([sys.executable, '-B', str(src / 'code/run_portable.py'), str(script), *map(str, args)], env=env, stdout=f, stderr=subprocess.STDOUT)
        assert rc == 0, (label, (out / (label + '.log')).read_text()[-2000:])

    shutil.copy(HERE / 'TRI-LEAVER-k8.json', out / 'TRI.json')
    shutil.copy(a.cells or HERE / 'CODED-CELLS-k8.json', out / 'CODED-CELLS.json')
    env['XPARAM'] = str(a.xparam or NOTRI / 'PRIO6.json')
    raw = out / 'raw'
    run('gauge-forward', a.cgx or HERE / 'cgx_leaver.py', wd / 'hamming', raw, out / 'TRI.json')
    word = raw / 'FULL-WORD.json.gz'
    chain = [parent / 'complex-local/static44/inputs/RAW504.json.gz', structural / 'global-replay/base/split-centres/FULL-WORD.json.gz',
             wd / 'shadow-early/FULL-WORD.json.gz', structural / 'multicold-eight-word/FULL-WORD.json.gz', wd / 'shadow-prefix/FULL-WORD.json.gz']
    for i, (b, c) in enumerate(zip(chain, chain[1:])):
        t = out / ('port%d' % i)
        run('port%d' % i, NOTRI / 'port_stage_x.py', '--shadow-input', b, '--shadow-output', c, '--input', word, '--map', raw / 'HELPER-MAP.json', '--gk-dir', gk, '--out', t)
        word = t / 'FULL-WORD.json.gz'
    w = rd(word); rg = {json.dumps(c['before']): c['after'] for c in rd(HERE / 'REGROUP-LEAVER-k8.json')['changes']}
    n = 0
    for part in ('A', 'B'):
        for i, e in enumerate(w[part]):
            k = json.dumps(e)
            if k in rg: w[part][i] = rg[k]; n += 1
    sys.path.insert(0, str(gk)); import gk_word as GW, gk_common as GC
    GW.V = GC.V = 112; GW.DATA = GC.DATA = 7168
    H = GW.profile(w)
    w['report'].update(local_profile={str(k): z for k, z in sorted(H.items()) if z}, local_rank_mass=sum(r * z for r, z in H.items()))
    (out / 'regroup').mkdir(); (out / 'regroup/FULL-WORD.json.gz').write_bytes(gzip.compress(enc(w), mtime=0))
    s = importlib.util.spec_from_file_location('upstream_hpe', up / 'code/build_hpe.py'); h = importlib.util.module_from_spec(s); s.loader.exec_module(h)
    used, ch = set(), []
    for c in h.candidates(w):
        rs = set(h.regs(w['cleanup'][c['cleanup_index']]))
        if rs & used: continue
        used |= rs; ch.append(c)
    w = h.apply(w, ch); (out / 'hpe').mkdir(); (out / 'hpe/FULL-WORD.json.gz').write_bytes(gzip.compress(enc(w), mtime=0))
    W = Word(w, wd / 'hamming'); res = []
    for b in range(14):
        d = W.pairing(b); res.append(dict(block=b, I=d['I'], rankC=d['rankC'], dimL=d['dimL'], dimK=d['dimK'], r=d['r']))
    rep = dict(variant=dict(xparam=str(a.xparam), cells=str(a.cells), cgx=str(a.cgx)), regroup_applied=n, hpe_changes=len(ch),
               hpe_sha256=hashlib.sha256(enc(w)).hexdigest(), R=w['report']['R'], pairing=res, total=sum(x['r'] for x in res))
    (out / 'PAIRING.json').write_text(json.dumps(rep, indent=1))
    print('regroup', n, 'hpe', len(ch), 'R', w['report']['R'], 'hpe sha', rep['hpe_sha256'][:12])
    for x in res: print('  block %2d |I| %d rankC %d dimL %d dimK %d rank(LK) %d' % (x['block'], len(x['I']), x['rankC'], x['dimL'], x['dimK'], x['r']))
    print('TOTAL', rep['total'])


if __name__ == '__main__':
    main()
