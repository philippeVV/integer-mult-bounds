# Full PG(2,16), non-reversible (expXX): the best schedule is |R| = 909, the completion fails, and the complex cap is below threshold

Status: exact checks plus heuristic search, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017. Inputs: `../contract-slack/` (expXW: non-reversible pricing, R4, the live lead), expXQ (Singer
cyclic method, min-cost-flow F, gate-by-gate verifier), expXL (split / cut-ascent identity), expNC (the F₁₆ Hermitian code,
dim 900), expRF (`markov.py`).

## Verdict: the full PG(2,16) host is closed for the Hermitian route

| Item | Value |
|---|---|
| Question | expXW's live lead: on the full PG(2,16) with non-reversible gates (points [17,8] 8→9, lines [17,9] 9→8), the Hermitian F₁₆ code has dim 900 and s = 900/4641 = 0.19392 > 0.19300. Is there a causal schedule good enough (briefed target \|R\| ≤ 906), and can the complex side follow? |
| **Best schedule** | **\|R\| = 909 = dim + 9**: Singer cyclic order m = 11, c = 91 (P_i at 2(11i mod 273), L_l at 2((11l − 91) mod 273) + 1, class orientation: points in-degree 8, lines 9). Verified gate by gate over F₁₆ (`sim16.py`): arcs are exactly the 4641 Singer incidences, Γ − R acyclic (Kahn), all 546 gate in-sets are information sets, **rank(Φ − I) = 9**, so dim W = 900 (an independent confirmation of expNC's dimension). Improves expXW's 1209 by 300. A second 909 schedule: m = 97, c = 182 |
| **The target \|R\| ≤ 906 was wrong** | R4's completion (c = \|R\| − dim dead sessions, zero-born at the head and erased at the tail, s_c = (dim − c)/(E − 2c)) is exact only if the dead set covers the support of Φ − I. Here **Φ − I has full support** (909 × 909, nnz 764,432): the minimum dead cover is **906** (CP-SAT, OPTIMAL), leaving **3** live sessions (s = 3/3735); the second 909 schedule gives 907 (2 live). Because the code is sound, the cut circuit's own code W′ ⊇ {w ∈ W : w_C = 0} gains nothing from splits, so **an exact schedule needs \|R\| = 900 exactly**, and for any local MDS codes on full PG(2,16) it needs F(σ) = dim ≤ Boot ≤ 901 |
| **Boot ∈ [900, 901]** | upper bound 901 (order SA `bsa.c`, 6·10⁹ moves from random, `b5.txt`), down from expXW's 907; lower bound 900 = dim C. The cut-ascent bound gives only 884 (below) |
| **Cost identity for unequal k (new)** | for any prefix S of σ, Σ_{v∈S}(k_v − j_v) = ½∂S + ½(L_S − P_S) (L_S, P_S = lines and points in S), hence **cost(σ) = ½(total ascent of the cut profile + #rising lines − #rising points)** and cost(σ) ≥ max_t ½(∂S_t + L_t − P_t). The cyclic order alternates P and L, so its rising points and lines balance (909 = ⌊1819/2⌋); the 901–902 Boot orders have the same max cut (≈ 1817–1819) and win by putting points in the rising half, which is exactly what breaks the back-arc transport (a 901 order has F = 1227). With the Laplacian spectrum {0, 13, 21, 34} and \|S\| = 273: min[∂S + L − P] ≥ 21·546t² + 546t + 1774.5 ≥ 1768, so Boot ≥ 884 |
| Search | Singer cyclic (all 144 units m × 273 c): 909, next 915; **exact DP over all merges** of the cyclic point and line sequences (all m, c, 17 cyclic out-windows; all 24,310 out-sets at (11, 91)): **909, provably optimal in this family**; two-level lex Singer orders (10⁷): 911; merge-path SA with exact F (1.7M F-evals): 909; joint order/orientation SA (3·10⁸ / 10⁹ moves): 909 from cyclic, 995 from random; soft-in-degree homotopy (2·10⁹): 909 / 1204; ILS with C primal-dual mcf: 909; CP-SAT LNS: contiguous 100-windows OPTIMAL at 909, seam blocks 909 (not proven), from the Boot-902 order 1287 → 1109; cost/F homotopy SA (2 × 2.2M exact F-evals): from the Boot-901 order F 1227 → **988**, from cyclic it leaves 909 at once and ends at 987. The unstructured landscape bottoms out about 80 above the Singer optimum |
| Excess trend | cyclic excess over dim: 0 (q = 9), **9 (q = 16, full)**, 31 (q = 25), 157 (q = 49) (`../exact-instances/`) |
| **Complex side** | expRF's Chebyshev–Markov cap with the free law of PQP at unequal traces α = 8/17, β = 9/17 (girth 6: moments j ≤ 2 are free, m₁ = αβ = 0.24913, m₂ = 0.18707; M = min{∫p dν : deg p ≤ 2, p ≥ 0 on [0,1], p(1) = 1}): **M = 0.18147 < 0.19300**, so a ≤ 0.0640, **κ ≤ 0.0602 (×35.4)** on this host even with a perfect schedule. The half-rate check reproduces expRF's 2/11 = 0.18182. Premise: flat-leverage local codes (expRF T2) |
| Multiplier | **×0** realised (no exact schedule; the 909 schedule completes to 3 sessions). Honest ceiling: exact schedule with F = 900 *and* a complex code at the Markov cap, κ ≤ min(0.103 bit, **0.0602 complex**) = 0.060 (×35), conditional on two open items, the first of which the search makes look out of reach. If R4's completion were valid (it is not): 909 would give s = 891/4623 = 0.19273, κ = 0.0991 (×58.3) |

## Reading

- The non-reversible lever of `../contract-slack/` does admit full PG(2,16) on the bit side by dimension, but a **dimension
  above threshold is not enough**: exactness needs F = dim, and the Singer family is exactly 9 short, with the whole of Φ − I
  carrying the defect. Dead-session completion cannot shave a sound code's excess unless the defect is localised.
- The complex side fails the Go condition **on its own**: no flat-leverage complex code on this host can reach 0.193, whatever
  the schedule. Reopening odd-degree full planes on the complex side needs non-flat local codes, which is a paper question first.

## Next steps and stop tests

1. **Close this host for the BH/Hermitian route.** Stop test met for Singer-structured orders (exact DP); open only for exotic
   orders.
2. If anything is tried, target an **exact pair (code, schedule)**, not \|R\| ≤ 906: local MDS codes (column scalings or twists
   per Singer class) for which Φ − I on the 909 schedule is supported on ≤ 9 rows or columns. **Stop** if random column
   scalings never bring the minimum cover below about 900.
3. Complex side: any reopening of odd-degree full planes must first beat M(8/17, 9/17) = 0.1815 with non-flat local codes.

## Files and reproduction

| Script (run from `xx/`) | Checks | Runtime | Output |
|---|---|---|---|
| `uv run -q --no-project --python 3.12 --with numpy --with scipy python -I -B markov17.py` | Markov caps: ½/½ = 2/11 (expRF check), 8/17, 9/17 = 0.18147 < 0.19300, κ_max = 0.0602; PG(2,16)−M = 2/11 | < 1 s | `markov17.out` (scouting `markov17.log` identical), `PASS_EXPXX_MARKOV` |
| `uv run -q --no-project --python 3.12 --with numpy --with ortools python -I -B sim16.py seed_11_91.txt sched909.json` | structure (incidences, in-degrees 8/9, \|R\| = 909, Kahn), 546 gates over F₁₆, rank(Φ − I) = 9, full 909 × 909 support, minimum dead cover 906 OPTIMAL (CP-SAT, one worker), (Φ − I)[S,S] = 0 on the 3 live sessions; order, arcs, R, rank and cover size equal to `sched909.json` | 3 s | `sim16.out`, `PASS_EXPXX_SIM909` |

Changes against the scouting copies: `sim16.py` and `code16.py` insert their own directory into `sys.path`; `code16.py`
caches `code16.pkl` next to itself (`run_all.sh` deletes it); `sim16.py` runs CP-SAT on one worker, compares with the stored
`sched909.json` instead of writing it, and prints a PASS marker; `markov17.py` appends a PASS check.

Included but not replayed (C via `gcc -O2 -lm`, Python via uv with numpy and ortools): `singer.py`, `common.py` (graph, cyclic
orders, Boot cost, exact F via ortools min-cost flow), `cyc16.py` (→ `cyc16.log`), `seed.py`, `seedo.py`, `pg.py`, `gf.py`
(expNC's code builder, used by `code16.py`), `bsa.c` (Boot SA; `b5.txt`, cost 901), `osx.c`, `osm.c`, `ils.c`, `msa.c`,
`mdp.c` (exact merge DP), `lex.c`, `hsa.c` (cost/F homotopy; `hsa_b5.txt`, `l_hsa_*.txt`), `lns.py` (CP-SAT LNS), `anal.py`;
inputs `D16.txt`, `graph16.txt`, `cyc_11_91.txt`; second schedule `sched909b.json` (m = 97, c = 182).

Source: expXX. Prepared with Claude Code (Claude Opus 5.5) assistance. The 909 schedule, its rank and the dead cover are exact;
the cost identity and the 884 bound are short author-checked arguments; the searches are heuristic except the merge DP (exact
within its family) and the CP-SAT window optimality; the Markov cap inherits expRF's premises. Not independently reviewed.
