# k = 7: every normal set with 7 distinct normals is closed

## Scope

**Result (model, notation and lemmas as in `../README.md`).** No nine-term masked P112 co-star
decoder has terms using exactly 7 distinct normals, over any field of characteristic ≠ 2, with zeros
and unequal left/right factors allowed. No new κ, decoder or cost claim.

The proof has two kinds of step. They do not have the same standing.

- **Exact, replayed certificates.** 7,534,138 of the 7,536,480 (orbit, composition) cases are closed by
  certificates: an initial L1 clique, or engine kills plus a final L0/L1/L1Z/L3odd contradiction.
  `verify_cert.replay` replays each one with plain Python sets. It imports neither `model.py` nor any
  engine. The orbit list is checked complete by brute force.
- **Two SAT solvers agreeing.** The other 2,342 cases are closed by Lemma S (`../lemmaS/README.md`)
  with capacity 3 allowed. Each is UNSAT on both CaDiCaL 1.5.3 and Glucose 4. No DRAT proofs were
  produced. This is the same standard as the published `lemmaS/` and Lemma C in `lemmas23/`.

k = 8 and k = 9 are **not** covered here. Of those sizes, only the normal sets containing F0 are closed
(`../f0/`). `../lemmaS/` closes the engine survivors of a sample of k = 9 sets, which is not a full
closure.

## Method

1. **Stage 1 (`stage1.py`, `ext7.py`, `merge1.py`).** Every 7-set is H-equivalent to B ∪ {x}, where
   B is one of the 32,292 6-orbit representatives in `../engine/orbk.pkl` and x ∉ B. Stage 1 runs
   over all (B, x). A non-covering set is recorded only by its canonical code: L0 at full supports
   kills every composition, because some live line lies in no co-star. A covering set gets the initial
   L1 filter over all 28 compositions of 9 into 7 positive parts, and the first A seen per canonical
   code (`canonfast.canon`, the lex-least sorted image under H) is kept. `merge1.py` writes
   `reps7.npz`.
2. **Stage 2 (`stage2.py`, `engineC.py`, `engineC3.py`, `cscan.c`).** For every covering representative
   and every one of the 28 compositions:
   - if the composition was killed by the filter, the certificate is `("L1", clique)` with no steps;
   - otherwise `engineC.propagate` runs (20 s Z-scan limit) and returns a certificate or SURVIVE.

   Every certificate is replayed in-process by `verify_cert.replay`. Stage 2 also records |Stab_H(A)|.
   `engineC` is `../engine/engine4.py` with the clique scan done in C. It gave identical output to
   engine4 on 309 sampled pairs in the original run, and on 46 more here (`cmp_engine4.py 1 30`).
   Correctness does not depend on this, because every certificate is replayed.
3. **Re-replay (`replay_s2.py`).** This reloads all certificates from disk and replays them again. It
   checks that each (rep, composition) appears exactly once and that no engine certificate is left
   over.
4. **Lemma S (`check_S7.py`, `check_S7s.py`, `s3run.sh`).** This runs on the 2,342 survivors
   (`cases_k7_all.json`). `check_S7.py` is `../lemmaS/check_S.py` with one change: capacity 3 is
   allowed, and in (R) a group touching Q contributes e_g = min(c_g, |S_g ∩ Q|) touch literals.
   `check_S7s.py` is the same code with the static (R) slack bound read from `SLACK` (default 2).

`model.py`, `group.py`, `cliques.py`, `canonfast.py`, `engine2.py` and `verify_cert.py` are
byte-identical copies of the files in `../engine/`. They are copied so that everything runs from this
directory. `ext7.py` keeps the run's `init`/`filt`/`filt_fast` unchanged; its unused prototype driver
was removed.

## Completeness of the case list

- `orb_indep.py` uses only numpy, `group_elems.npy` and `reps7.npz`. By brute force over all 172,032
  elements of H, it checks that every one of the 419,952 representatives is the lex-least sorted
  image of itself, so the representatives are pairwise inequivalent. It also counts each stabiliser
  and checks Σ |H|/|Stab| = C(120,7) = 59,487,568,920 (`logs/orb_indep.log`: MATCH). Together these
  show that the list has exactly one representative per H-orbit of 7-sets. This argument does not
  depend on stage 1 or on the completeness of the 6-orbit list.
- `orbsum.py` is a second route that uses the stage-2 stabilisers. It checks that H is a set of
  distinct F0-preserving |Gram|-isometries containing the identity. It confirms the covering /
  non-covering split with a direct union test, and it gets the same sum (`logs/orbsum.log`):
  covering 37,720,632,568 sets, non-covering 21,766,936,352.
- H preserves the mask and the co-star structure. So closing a representative for all 28
  compositions closes its whole orbit. Non-covering sets fail L0, since some diagonal 4 is
  unreachable.

## Counts

| quantity | count |
|---|---:|
| 6-orbit reps B × lines x ∉ B (raw extensions) | 3,681,288 (1,305,898 non-covering, 2,375,390 covering) |
| raw (covering extension, composition) pairs surviving the L1 filter | 32,010,853 (34,500,067 killed) |
| 7-orbits | 419,952 (269,160 covering, 150,792 non-covering) |
| covering orbits with ≥ 1 filter survivor; surviving (orbit, composition) pairs | 261,919; 3,636,308 |
| stage-2 cases (269,160 × 28) | 7,536,480 |

Stage 2 by capacity shape (`logs/summarize.log`):

| shape | cases | initial L1 | L0 | L1 | L1Z | L3odd | survive |
|---|---:|---:|---:|---:|---:|---:|---:|
| (3,1⁶) | 1,884,120 | 1,222,909 | 319,384 | 341,732 | 86 | 0 | **9** |
| (2,2,1⁵) | 5,652,360 | 2,677,263 | 1,249,398 | 1,709,807 | 13,316 | 243 | **2,333** |
| total | 7,536,480 | 3,900,172 | 1,568,782 | 2,051,539 | 13,402 | 243 | **2,342** |

Certificates: 7,534,138. All replayed in stage 2, and again from disk by `replay_s2.py`: 7,534,138
OK, 0 failures, each (rep, composition) exactly once, 269,160 reps (`logs/replay_s2.log`).

The survivors are **not** all of the A2 type seen in the sampled run (`../lemmaS/`):

| survivor type | count |
|---|---:|
| (2,2,1⁵), the two capacity-2 normals span an A2 (non-orthogonal) | 2,111 |
| (2,2,1⁵), the two capacity-2 normals orthogonal | 222 |
| (3,1⁶) | 9 |

Lemma S on the 2,342 survivors:

| run | cases | result |
|---|---:|---|
| `--frames` (1,260 live 8-frames), `logs/s3/k00–k23.log` | 2,342 | 2,332 UNSAT on both solvers, 10 SAT on both |
| the 10 (`cases_k7_satframes.json`), all cliques with slack ≤ 2, `logs/s3/sat10_full.log` | 10 | 8 UNSAT on both, 2 SAT on both |
| the last 2 (`open1.json`, `open2.json`), slack ≤ 4 (192,714 rank cliques, 9,837,878 clauses), `logs/s3/open{1,2}_slack4*.log` | 2 | UNSAT on both, in the original run and in an independent re-run |

The last two cases are the set (0,2,3,14,15,56,80) with capacity 2 on {3,15} or on {2,14}. They are
equivalent: the set's stabiliser in H has order 32, and 16 of its elements map {3,15} to {2,14}.

The earlier estimate in `../lemmaS/README.md` (27.93M raw pairs, about 3.2M after dedup) came from an
incomplete stage-1 run. The true figures are 32.01M raw filter-surviving pairs and 3,636,308 after
orbit dedup.

## Files

| file | content |
|---|---|
| `reps7.npz` (1.6 MB) | `cov_code`, `cov_A` (269,160 covering reps), `cov_mask` (28-bit filter-survivor masks), `non_code` (150,792 non-covering reps, 7 × 7-bit codes) |
| `cases_k7_all.json` | the 2,342 stage-2 survivors `[normals, caps]` |
| `cases_k7_satframes.json` | the 10 that are SAT with `--frames` |
| `open1.json`, `open2.json` | the 2 that need slack ≤ 4 |
| `logs/` | stage 1, merge, stage 2 (several resumed sessions), summarize, orbsum, replay, orb_indep, and the Lemma S logs (`s3/`) |
| `CERTS-MANIFEST.json` | sha256 and size of each uncommitted output (see below) |

**Not committed (size):**
- `s2/e*.xz`: engine certificates, 318 MB.
- `s2/d*.npz`: per-pair table and stabilisers, 31 MB.
- `s2/v*.pkl`: survivors, 0.6 MB.
- `s1/c*.npz`: 9 MB.
- `group_elems.npy`: 20.6 MB.
- `cliques.pkl`: 3.7 MB.

`CERTS-MANIFEST.json` lists each of these files. `python3 group.py` regenerates `group_elems.npy`;
the result was byte-identical here (sha256 `3510f278…54af`). `python3 cliques.py` under Python 3.14
regenerates `cliques.pkl`, also byte-identical. Under uv Python 3.12 the pickle bytes differ, but the
clique list is the same. The stage outputs are regenerated by running the stages again. The bytes of
`e*.xz` depend on the lzma and pickle versions and on the engine's 20 s Z-scan limit. The invariant is
the replayed counts. On a much slower machine the limit could turn some cases into extra survivors,
which would then need Lemma S.

## Reproduce

Run everything from this directory. You need Python ≥ 3.12 with numpy, plus python-sat for
`check_S7*.py`, and a C compiler. The original runs used
`uv run -q --no-project --python 3.12 --with numpy --with python-sat python …`, on an 8-core Ryzen
AI 7 PRO 350. Cap memory with something like `systemd-run --user --scope -p MemoryMax=5G`.

```sh
python3 group.py && python3 cliques.py           # group_elems.npy, cliques.pkl (seconds)
cc -O2 -shared -fPIC -o cscan.so cscan.c
python3 stage1.py 3      # ~40 min with 3 workers -> s1/ (reads ../engine/orbk.pkl; env ORBK overrides)
python3 merge1.py        # -> reps7.npz (compare with the committed copy)
python3 stage2.py 3      # ~11 CPU-h -> s2/ (resumable; prints replayFAIL count, must stay 0)
python3 summarize.py     # tables above; rewrites cases_k7_all.json
python3 replay_s2.py 3   # ~25 min: REPLAY ALL OK [7534138, 0, 2342, 7536480, 269160]
python3 orbsum.py        # orbit sum via stage-2 stabilisers: MATCH
python3 orb_indep.py . 3 # independent completeness: all lex-least True, MATCH
                         # (not timed; a 300-rep slice takes ~3 s, extrapolated ~1-2 CPU-h)
PY=python3 ./s3run.sh 0 & PY=python3 ./s3run.sh 1; wait   # Lemma S --frames, ~1 h on 2 workers
python3 check_S7s.py cases_k7_satframes.json              # slack 2: 8 UNSAT, 2 SAT
SLACK=4 python3 check_S7s.py open1.json                   # 9.8M clauses, UNSAT on both
SLACK=4 python3 check_S7s.py open2.json
```

Spot checks run on the committed copies for this package (in a scratch copy of the package):
- regenerated `group_elems.npy`: sha matches.
- `stage1.py` chunk 0 output identical to the original `s1/c0000.npz`.
- `stage2.py` on the first 20 reps: rows and stabilisers identical to the original `s2/d0000.npz`.
- `replay_s2.py` on original chunk 0: 28,000 OK.
- `orb_indep.py` slices [0, 300) and [300000, 300300): all lex-least.
- `check_S7s.py --frames` on 4 cases, one with caps (3,1⁶): lines identical to the original logs.
- slack-2 run on a frames-SAT case: UNSAT, as in `sat10_full.log`.
- `cmp_engine4.py 1 30`: 46 pairs, 0 differences.

The full runs were not repeated.
