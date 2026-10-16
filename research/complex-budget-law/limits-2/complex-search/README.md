# Complex search (expXP): a spectral Boot window on PG(2,q), and complex MDS codes that fall away from it

Status: one field-free theorem with exact checks, plus a numerical search, author-checked, LOCAL-ONLY until published. No κ
claim. Base: main `0fdd7d82`. Multipliers are over κ = 0.0017. Inputs: `../../limits-1/` (complex-dual, bootstrap,
seeds-negative), expXE, expXG, expXK, expXM.

## Verdict

| Item | Value |
|---|---|
| Question | Problem X(ℂ, a) on the smallest hosts: can **non-equivariant** complex MDS codes, protected against gauge artefacts, follow Boot into the window (E/g, Boot]? |
| **Spectral Boot bound** (proved, any field) | **Boot(PG(2,q)) ≥ λ₂·n/8 = (q+1−√q)(q²+q+1)/4**, i.e. **Boot/E ≥ ¼(1 − √q/(q+1))**. Proof: cost(σ) = ½·(total ascent of the cut profile) ≥ ½·cut(middle prefix) ≥ ½·λ₂·(n/2)(n/2)/n, and λ₂ = q + 1 − √q for the incidence graph (adjacency spectrum ±(q+1), ±√q) |
| Room above E/6 | the certificate gives Boot > E/6 **exactly for q ≥ 7** (√q/(q+1) < 1/3 ⇔ q > 6.85). PG(2,7) is the smallest graph with Boot/E > 1/g at girth 6: Boot ∈ [77, 79] vs E/6 = 76. PG(2,5) is closed: Boot ≤ 31 = E/6 |
| Room above the a = 1/9 threshold | certified above thr = 1/(1+5^{8/9}) = 0.19300 **exactly for q ≥ 19** (PG(2,19): Boot/E ≥ 0.1955). PG(2,17) misses by 1.5 sessions (LB 1065.05 < 1066.54 ≤ UB 1095): still "not excluded" |
| **Complex search** (`cmds.py`) | Riemannian ascent on ∏Gr(k, 2k) with the canonical gauge (P-Barthe-flat, then doubly normalised) re-imposed every step and an MDS barrier per code in its own Barthe gauge. Healthy (f → 0, flattened margin ≥ e^-3) codes reach **r ≈ 16–17 on PG(2,5)** (E/6 = 31, Boot 30–31) and **r ≈ 30 on PG(2,7)** (E/6 = 76, Boot 77–79); r = 40, 50, 60, 77 on PG(2,7) fail (f ≥ 3.7e-3, margin ≤ e^-5.6; at 77 the L-leverages go to {0, 1}). PG(2,9) at r = 152: flat at f = 7.2, stopped |
| **r/E falls with q** | 0.154 (PG(2,3), calibration, = Boot), 0.091 (PG(2,5)), 0.066 (PG(2,7)): generic complex MDS codes track roughly E/(2q + …), not Boot; the numerical ceiling sits below E/6 by ×2.5 on PG(2,7). Same rigidity as expRB's frequency codes (about 1/(4q)) and expXM's Singer codes |
| Gauge artefacts | without the canonical gauge, PG(2,5) at r = 31 reaches f = 4e-5 in its own gauge; re-gauged canonically f = 1.46 (Lemma D on data). Without the barrier, r = 20, 31 "converge" to Plücker margins e^-23, e^-29 (non-MDS). A reverse continuation at r = 20 (barrier 1e-12 → 1e-3) raises f monotonically from 1e-4 to 2.4e-2: no MDS point on the path |
| Girth 10 | degree-4 voltage lifts of K₄,₄ have Boot/E ≤ 0.092 (QC4_107, E = 1712) and ≤ 0.090 (QC4_97): below 1/10, excluded in every field. Degree 6 needs E ≥ 4686 (Moore bound); the greedy QC search found no girth-10 lift with L ≤ 560 in 40 min |
| Multiplier | **×0**. Stop test fired: no certified complex code beats E/6 with a positive MDS margin anywhere, so no GRS screen or schedule check was reached |

## Bootstrap table

| graph | g | d | E | greedy/SA UB | spectral LB | LB/E | E/g | thr(1/9)·E |
|---|---|---|---|---|---|---|---|---|
| PG(2,5) | 6 | 6 | 186 | 31 | 29.2 → 30 | 0.157 | 31 | 35.9 |
| PG(2,7) | 6 | 8 | 456 | 79 (expXE) / 80 | 76.3 → **77** | 0.167 | 76 | 88.0 |
| PG(2,9) | 6 | 10 | 910 | 166 (= BH dim, exact) | 159.3 | 0.175 | 151.7 | 175.6 |
| PG(2,11) | 6 | 12 | 1596 | 300 | 288.7 | 0.181 | 266 | 308.0 |
| PG(2,13) | 6 | 14 | 2562 | 491 | 475.5 | 0.186 | 427 | 494.5 (UB below: excluded) |
| PG(2,17) | 6 | 18 | 5526 | 1095 | 1065.05 | 0.1927 | 921 | 1066.5 (undecided) |
| PG(2,19) | 6 | 20 | 7620 | 1527 | **1489.8** | **0.1955** | 1270 | 1470.7 (**certified above**) |
| PG(2,23) | 6 | 24 | 13272 | 2715 | 2655.0 | 0.2000 | 2212 | 2561.5 (certified above) |
| QC4_97 | 10 | 4 | 1552 | 140 | 69.6 | 0.045 | 155.2 | 192.8 (excluded) |
| QC4_107 | 10 | 4 | 1712 | 158 | 72.9 | 0.043 | 171.2 | 212.7 (excluded) |

## Reading

The window (E/6, Boot] exists from q = 7 on and the a = 1/9 window from q = 19 on, so the field-free ceiling does not close the
complex side on large planes. What closes the small planes is the codes: on the one graph with the cleanest window (PG(2,7),
77 > 76) random non-equivariant complex MDS codes fail at r = 40, about half of E/6, and the reachable fraction falls with q.
This is numerical evidence, not a proof; the proved statements on the complex side are in `../theta-partial/` and
`../../limits-1/complex-dual/`.

## Next step and stop test

1. Do not run larger unstructured complex searches on PG(2,q). **Stop line:** reopen only if a structured, non-flat code family
   reaches r ≥ 0.12·E on PG(2,7) (the best complex value anywhere is 0.156, at PG(2,3)).
2. PG(2,17): a Boot ≥ 1067 certificate needs one unit more than λ₂ gives (for example a bisection SDP with bipartite cuts).
   Low value: complex codes are far from Boot regardless.
3. Girth 10: a degree-6 girth-10 graph with E ≈ 5–10k needs a better construction (PEG or CD(5,q) components).

## Files and reproduction

| Script (run from `xp/`) | Checks | Runtime | Output |
|---|---|---|---|
| `uv run -q --no-project --python 3.12 --with numpy python -I -B spectral_boot.py` | λ₂ = q + 1 − √q on PG(2,p), p prime ≤ 23 (homogeneous coordinates, dense eigvalsh); closed form = `boot_pg.log`; over all prime powers q ≤ 10⁴: Boot > E/6 certified iff q ≥ 7, above 0.19300 iff q ≥ 19; PG(2,7) window 77 > 76 | 2 s | `spectral_boot.out`, `PASS_EXPXP_SPECTRAL` |
| `python3 -I -B ladder_summary.py` (stdlib) | re-reads the logged ladders: best healthy r = 17 on PG(2,5), 30 on PG(2,7); failing rungs [18, 19] and [40, 50, 60, 77]; r/E decreasing | < 1 s | `ladder_summary.out`, `PASS_EXPXP_LADDERS` |
| `PASS_F=1e-10 uv run … python -I -B cmds.py ../../../limits-1/complex-dual/xg/graphs_small.json PG3 8 0` | calibration: PG(2,3), r = 8 = Boot: f = 5e-14, margin e^-0.30, gap 0.293 | 3 s | `cmds_pg3_r8.out`, `PASS_EXPXP_CMDS_PG3_r8` |
| `PASS_F=1e-6 uv run … cmds.py ../../../limits-1/complex-dual/xg/graphs_small.json PG5 16 4 4000` (slow) | one healthy PG(2,5) rung: r = 16, seed 4 | 30–42 s | `cmds_pg5_r16.out`, `PASS_EXPXP_CMDS_PG5_r16` |

`spectral_boot.py` and `ladder_summary.py` are new; `cmds.py` gained only the optional `PASS_F` check at the end (converged
below PASS_F with flattened margin ≥ e^-3). The PG(2,5) rung replays to f = 5.0e-7 and margin e^-1.56 against the logged
f = 1.0e-7 and e^-1.72 (not bit-identical; same verdict). The `graphs_small.json` PG3/PG5 entries are identical to expXP's
`graphs.json`. The upper bounds in the table are expXE/expXP greedy values (`boot_pg.log`); `bootsa.py` (needs numba) and
expXP's 1 MB `graphs.json` are not replayed (the graph file is not copied).

Included but not replayed: `bootsa.py` (greedy plus float-key SA upper bound and the spectral bound), `qc10.py` (girth-≥ G voltage
lifts of K_{d,d}) with `qc4.log`, `qc4b.log`, `qc6.log`; logs `boot_pg.log`, `ladder_pg5.jsonl`, `ladder_pg7.jsonl`,
`ladder_pg7b.jsonl`, `cont_pg5_r20.err` (reverse continuation), `t_pg5_c{0,1}.out` (canonical re-gauge of the own-gauge hit),
`pg9_r152.err`.

Source: expXP. Prepared with Claude Code (Claude Opus 5.5) assistance. The spectral bound is proved; the search is numerical
evidence. Not independently reviewed.
