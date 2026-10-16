# Seeds, negative (expXD, expXM, expRJ; expRB/expRF in closures-3): hexagons, the PG(2,17) Singer complex, and expander ceilings

Status: scoped closures, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`. Multipliers are over
κ = 0.0017. Metric: a = 1 − ln((E−S)/S)/ln b with b = g − 1, κ ≈ a/(1+a).

Cross-reference: the char-0 frequency (ζ) and Walsh–Hadamard seeds (expRB, expRF) are closed in
[`../../closures-3/complex-coding-seeds/`](../../closures-3/complex-coding-seeds/) and are not repeated here; its T1–T3 are
extended to all leverage profiles in `../complex-dual/`.

## Verdict

| Lane | Verdict | Realised | Ceiling left |
|---|---|---|---|
| Hexagons H(q) (expXD) | **closed at constructible size**: natural [4,2] code on H(3) has S = 128 (S/E = 0.0879, a = 0.024); H(5) [6,3] S = 2263 (0.0966, a = 0.068); H(4)−M 163, and the x² twist *lowers* it to 119; random gates ≤ 3. Required for a ≥ ½: S/E ≥ 0.232. Stop test (S/E < 0.15) fires | ×0 | conditional ×14 (H(3)) / ×37 (H(5), ×18 under the κ < 1/33 guard); the first genuine Hermitian twist is H(9), E = 664,300, not recommended |
| PG(2,17) Singer complex (expXM) | **closed for Singer-equivariant codes**: every exact equivariant family gives S ≤ 81 = k² (Chebotarev rigidity at prime N = 307); numerical MDS solutions saturate at 120 rank-1 blocks (≈ parameter count 145); the stop line is M₆·E = 1005, Boot = 1095. Larger "solutions" are provably non-codes (equivariant-degeneracy lemma) | ×0 | non-equivariant complex codes on PG(2,17) (≈ 50k unknowns) are untested |
| Expander ceilings (expRJ) | **cut/flow ceilings do not close the level-scaling question**; what decides it is the existence of codes. Tanner hybrids: **a ≥ ½ ⇔ S/E ≥ 1/(1+√(g−1))** (at S/E = ¼ this means girth ≥ 10). Far pairs: a_H ≤ 0.458 within 10⁵ roles (LPS X^{5,37}); ½ needs ≳ 10⁷ nodes. Lemma C forbids compiled seeds (flow–cut gaps ≤ 1.43 vs > 3 needed) | ×0 | conditional ×185 (far pair), ×196–207 (LPS(5,37), H(3) at S/E = ¼), ×147 under the unfused outer |

## Hexagons (expXD)

**Construction.** Q(6,q): X₀X₄ + X₁X₅ + X₂X₆ = X₃²; lines of H(q) are totally singular lines whose Grassmann coordinates satisfy
p₁₂ = p₃₄, p₅₄ = p₃₂, p₂₀ = p₃₅, p₆₅ = p₃₀, p₀₁ = p₃₆, p₄₆ = p₃₁ (Van Maldeghem 2.4.13). Verified for q = 2…5: (q⁶−1)/(q−1)
points and lines, order (q, q), incidence girth 12. Local codes: homogeneous forms of degree k − 1 in canonical coordinates
(G₂(q)-equivariant, MDS), optionally Frobenius-twisted.

| seed | E | S | S/E | a (b = 11) |
|---|---|---|---|---|
| H(3), natural [4,2] (F₃, F₉, F₂₇, dual variants) | 1456 | 128 | 0.0879 | 0.024 |
| H(3), random gauge / random MDS gates | 1456 | 0–3 | ≤ 0.002 | < 0 |
| H(4)−M, natural [4,2] (F₄, F₁₆) | 5460 | 163 | 0.0299 | −0.45 |
| H(4)−M, x² twist on one side | 5460 | 119 | 0.0218 | −0.59 |
| H(5), natural [6,3] over F₅ | 23436 | 2263 | 0.0966 | 0.068 |

H(3) passes the cut test (SA bisection 254, so S ≥ 127 is needed and 128 is available), but the exact-schedule search reached
|R| = 390 against 128. **Char-0 analogue:** the vertex stabiliser in G₂(q) induces 2-transitive PGL(2,q) on the q + 1 edges,
so equivariant char-0 codes are flat and Thm F of `../complex-dual/` gives S/E ≤ 1/12 at girth 12: a ≤ 0 exactly. The char-p
codes beat 1/12 only marginally (semicontinuity excess, not a mechanism).

## PG(2,17) Singer complex (expXM)

Singer model N = 307, D normalised with 17·D = D (six ×17-orbits of size 3), [18,9] codes, E = 5526. For equivariant pairs
(A, B): **dim W = Σ_{χ∈Z₃₀₇} dim(A ∩ Λ_χB)**, Λ_χ = diag(ζ₃₀₇^{χd})_{d∈D}, each block 18 × 18.

| search | best dim W | S/E | status |
|---|---|---|---|
| exact structured families (112 codes; frequency, Galois-twisted GRS with Gaussian-period weights, pencil GRS lifts, Gaussian-period matrices, semilinear) | **81** = k² | 0.0147 | exact mod two split primes (agree on all rows); MDS checked on all 48620 minors |
| numerical MDS, rank-1 block targets (det-Newton) | **120** | 0.0217 | numerical, spectral gap ≈ 3·10⁻⁴; 130 targets fail (0 of 4 seeds) |
| rank-2/3/4 targets | 20 / 18 / 16 | ≤ 0.004 | numerical |
| degenerate attractors | j·307, j = 2…7 (up to 2149) | — | **not codes** (lemma below) |

The char-0 "Frobenius twist" is either trivial (u ∈ ⟨17⟩: a coordinate permutation, S = 81) or kills everything (other u:
S = 0). Fixed-field rows are constant on ×17-orbits (rank ≤ 6), so fixed-field codes must break the multiplier symmetry and then
behave generically.

**Equivariant-degeneracy lemma (rigorous).** If A ∩ B contains a nonzero Λ-invariant subspace, it is a coordinate subspace F^X
(the ζ^d are distinct), and then dim ≥ |X| in every block, so S ≥ |X|·N. An MDS [18,9] code has minimum distance 10, so it
contains no vector supported on ≤ 9 coordinates: such X cannot occur. Hence every "solution" with a constant block profile
dim = j at all 307 χ (S = 614 … 2149, including 1228 > 1005) is a non-MDS boundary point (measured Plücker ≤ e^{−36}). This is
the equivariant form of Lemma D and of the bootstrap boundary; any hit from a loss without an MDS barrier must be discarded. ∎

Reading: Singer-equivariant complex codes stay ×8.4 (numerical) to ×12 (exact) below the stop line, while char 2 nearly reaches
the field-free ceiling on PG(2,16) (900 vs 914): the char-0 rigidity of expRB, now measured against Boot.

## Expander ceilings (expRJ)

Ceiling forms (W = r·k roles, E internal edges, b_min the minimum own-pair distance): a_H = 1 − ln(E/(rk))/ln b_min (the
repo dictionary made consistent; primary), a_G = ln(σ̂/λ*)/ln b_min (an envelope, not a cost).

- **Far pairs:** a_H ≈ 1 − ln c_d/ln(log_{d−1} n − 1), c_d ≈ 6.5, 4.3, 3.6, 3.3, 2.9, 2.7 for d = 3, 4, 5, 6, 8, 10. a_H ≥ ½ needs
  b_min ≥ c_d² ≥ 7.5, i.e. ≳ 10⁷ nodes. Best within 10⁵ roles: LPS(5,37) 0.458, LPS(5,41) 0.433, rr6 0.431; first ≥ 0.317 at
  ≈ 2k roles (LPS(5,13) 0.336).
- **Tanner hybrids:** cut and flow ceilings are vacuous for designed session sets (S/E up to ½); random balanced sets bind at
  S/E ≈ 0.3. **No graph-theoretic quantity caps S/E in (¼, ½); only the algebra of the local codes does.** At S/E = s,
  a = 1 − ln((1−s)/s)/ln(g−1), so a ≥ ½ ⇔ (1−s)/s ≤ √(g−1) ⇔ **s ≥ 1/(1+√(g−1))**.
- **Lemma C (expNN):** a compiled seed needs a per-copy gap σ/λ* > 3; measured flow–cut gaps are 1.0 (n ≤ 128) to 1.43
  (n = 512). Only single-copy exact-schedule seeds can carry these ceilings.
- Hand-offs: H(3) (girth 12, a = 0.542 at S/E = ¼) and H(4)−M were screened by expXD above (closed); LPS X^{5,37} (girth exactly
  10, 50,616 gates) has no known code mechanism.

## Next step and stop test

None of these lanes should be continued. Reopen triggers: (i) hexagons — a proof or disproof that the natural H(q) code has
s → 1/12 + o(1); (ii) PG(2,17) — an MDS-barrier-protected non-equivariant solver reaching r = 1006 with a positive flat-gauge
Plücker margin (stop if unreachable); (iii) expanders — any linear multiple-unicast code with rate near the sparsest cut on a
bounded-degree expander (expL7's extremal problem).

## Files and reproduction

Runtimes under a 4G cap, one worker.

| Script (run from its subdirectory) | Checks | Runtime | Output |
|---|---|---|---|
| `xd/`: `uv run -q --no-project --python 3.12 --with numpy --with numba python -I -B check_xd.py` | H(3) build (364 + 364, girth 12), natural code S = 128 over F₃ and F₉, random gates ≤ 3; H(4)−M girth 12, S = 163, twist 119 | 5 s | `check_xd.out`, `PASS_EXPXD` |
| `xd/`: `uv run … --with numpy --with numba python -I -B h5.py 5` (slow) | H(5) natural [6,3], blocked rank mod 5: S = 2263 | 43 s | `h5.out`, `PASS_EXPXD_H5` |
| `xm/`: `uv run … --with numpy python -I -B fam.py 0` and `fam.py 1` | 112 exact equivariant families mod P₀ and P₁; max S = 81 | 3 s, 5 s | `fam_P0.out`/`.json`, `fam_P1.out`/`.json`, `PASS_EXPXM_FAM` |
| `xm/`: `python3 -I -B agree.py` | the two primes agree on every row | < 1 s | `agree.out`, `PASS_EXPXM_AGREE` |
| `xm/`: `uv run … --with numpy --with scipy python -I -B lm.py m1:120 3` | det-Newton reaches 120 rank-1 blocks with an MDS Plücker margin (min log ≈ −12) | 1 s | `lm_120_s3.out`, `lm_m1_120_s3.npy`, `PASS_EXPXM_LM` |
| `rj/`: `python3 -I -B summarize.py` | tables from the expRJ sweep data (`big_rr.jsonl`, `big_struct.jsonl`, `rr.jsonl`, `tanner.jsonl`) | < 1 s | `summary.out`, `PASS_EXPRJ_SUMMARY` |

Included but not replayed: `xd/` `scan3.py`, `scan4.py`, `info.py`, `hill3.py`, `snd.py`, `sched.py` and their logs (twists,
non-square gates, hill-climb, soundness, schedule SA, `cut3.log`); `xd/gf.py` is a copy of expNC's field tables. `xm/` `chart.py`,
`inc.py`, `numopt.py`, `ref.py`, `jrank.py` and the sweep/continuation logs (`sweep_*.log`, `inc_*.log`). `rj/` `cc2.py`
(graphs, cuts, HiGHS max-concurrent flow), `sweep.py`, `big.py`, `girth.py` (girths; `girth.json`); these need networkx/scipy and
take minutes to hours.

Sources: expXD, expXM, expRJ (expRB, expRF via closures-3). Prepared with Claude Code (Claude Opus 5.5) assistance. Not
independently reviewed; the H(q) geometry facts are standard and were verified by computation for q ≤ 5.
