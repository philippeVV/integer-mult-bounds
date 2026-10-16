# Laminar, helper-free decoders (idea 6): no laminar decoder compresses (realised ×1.00, side effect ≤ ×1.05)

Status: scoped closure, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017.

## Verdict

| Item | P120 | P112 | P104 | E7 |
|---|---|---|---|---|
| known decoder (G/2 type), C/V | 72/120 | 64/112 | 56/104 | 49/63 |
| its corrections (every free pair, both directions) | 6720 | 5824 | 4992 | 2016 |
| laminar? | **no** | **no** | **no** | **no** |
| edge-disjoint bowties (each alone forbids a laminar order) | 465 | 407 | 344 | 136 |
| corrections a greedy laminar schedule makes direct | 441 (6.6%) | 406 (7.0%) | 374 (7.5%) | 232 (11.5%) |
| best laminar G/2-gauge decoder | γ = 1 | – | γ = 1 | γ = 1 |
| span MILP over the compressing pieces (P120 pool), three laminar allowances | **C = 120 = V (identity), proven optimal in the family** | | | |

- **Realised: ×1.00.** No laminar decoder with γ < 1 exists in any searched family.
- **Side effect:** making the greedy 6.6–7.5% of corrections direct on today's word is worth **≤ ×1.04–1.06** at s = 9
  (uncharged envelope); s = 3 stays at ×0.42.
- Not yet a theorem for arbitrary rational decoders (see "Open").

## Changed hypothesis

Idea 6 (expL5): minimise C subject to the corrections (the free off-diagonal nonzeros of N) admitting a **laminar order**, in
which every correction is a direct gate t_a += c·x_b with no relay or helper. This generalises expNT's reader-chain condition
to arbitrary correction graphs. The decoder objective changes from "min C" to "min C subject to a chain order", and the helper
model changes from relays to direct gates.

## Model (`lam.py`, `dfs.py`)

- Source x_b starts at ⟨p_b⟩, target t_a at ∅. A correction (b → a) needs (a, b) free and joins frames (OR-reach).
- Legality: t_a stays inside p_a^⊥, so every label c that reaches t_a along a time-increasing path must have B_ca = 0.
- The relaxation favours the idea (corrections get their own phase; sources carry no centre or gather frames; targets no
  deliveries). Infeasibility in this model therefore holds for every real schedule. A subgraph of a laminar correction set is
  laminar, so any obstruction contained in E kills E.

## Results

**Bowtie obstruction (certified).** Centre c, two free triangles {c, a, b} and {c, d, e}, the four cross pairs required. With
all 12 corrections, no laminar order exists: `dfs.py` proves it exhaustively without a solver, and CP-SAT (`lam.py`) agrees.
The minimal directed sub-obstruction has 9 corrections (`minbow.py`). Every one-way orientation of the bowtie is laminar, as
are a single free triangle and bidirected stars.

**The compressing decoders are far from laminar (`maxlam.py`).** G/2-type decoders have N_ab = ±½ on every free pair, so E is
the whole free graph, with 52–56 corrections per port and hundreds of disjoint bowties. A greedy laminar schedule makes only
about 3.6–3.7 corrections per port direct (a lower bound on what can be made direct).

**Gauge family (scouting, not replayed).** Block-gauge decoders ⊕ G/2|_{L_i} are laminar iff every block's full free graph is
laminar. Randomised growth finds laminar-closed sets of at most 14 lines (P120, P104) and 11 (E7), with gauge gain −12 to −20.
The puncture ladder (delete the line in the most bowties) reaches gain 0 at V = 64 / 56 / 42 and the first laminar set at
V ≈ 10–11. An anneal with laminarity as a hard constraint found nothing better than gain −1.

**Span MILP (`span_milp.py`).** N = Σ_j λ_j G/2|_{L_j} + diag(δ) over a pool of P120 (72), 60 random P112_F (64) and all 120
E7_α (49), with laminar allowances none / the greedy schedule's 13 bidirected pairs / 40 disjoint free triangles. HiGHS proves
the optimum is the identity, C = 120 = V, in all three. Every line must be separated from each of its 56 free partners by
pieces containing it but not the partner; pieces with gain contain thousands of free pairs, and cancelling them costs more
than the gain.

## Rank lemma (proved, elementary) and what is open

For helper-free decoders (every term has ρ_k ≥ 1): **C ≥ #terms ≥ rank_Q N ≥ V − fvs(E)**, where fvs is the directed
feedback-vertex number of the correction graph (2-cycles count). If E[I] is acyclic, N[I, I] is unit triangular in a
topological order. So one-way correction graphs cannot compress, and γ needs fvs(E) ≥ (1 − γ)V. Laminarity forbids bowties but
allows matchings and disjoint free triangles (fvs = V/2, 2V/3), so rank alone gives only γ ≥ ⅓–½.

**Open:** an arbitrary rational N whose free support E is sparse (≤ about 4 per port) and bowtie-free while N is cheap. This is
expL2 §4's question ("does any compressing decoder have d̄ ≪ V/2?") with the sharper target "bowtie-free E".

## Pricing (`price.py`)

expL2 pair model: s = 3 is expNE's η₃ with helper mass ρ; s = 9 is expNB's cal9 and 'nohelp'; bit twin at h = 8, complex at
h = 9, x = γ on both sides. The control reproduces today (s = 9 cal9 ×1.009, s = 3 ×0.40).

| scenario | s = 3 | s = 9 | best |
|---|---|---|---|
| realised: no laminar compressing decoder | – | ×1.00 | **×1.00** |
| E8 today plus partial laminarity (6.6–7.5% direct; helper mass × (1 − f)) | ×0.42 | ×1.04–1.06 over today | **≤ ×1.05**, uncharged |
| hypothetical laminar E8-type, γ = 4/7, ρ = 1/16 | ×1.35 | ×2.40 (nohelp) | ×2.40 |
| γ = 0.65 / 0.70 / 0.80 | ×0.24 / ×0 / ×0 | ×1.57 / ×1.15 / ×0.54 | |

At s = 3 the deficit is D = 2 − 3γ, so the E8-type s = 3 rescue needs γ well below ⅔ and ρ ≈ 0. Direct gates absorb at most
about 7% of transfers, taking ρ from 2.31 to about 2.16, far from the ρ ≤ 1.5 the s = 3 route needs.

## Next step and stop test

1. **Close idea 6 as a decoder lever** (recommended).
2. (Optional) Hand "E bowtie-free, ≤ 4 per port" to the correction-density LP (expL2 §5). **Stop:** if it forces d̄ ≥ 5 at
   γ ≤ 0.7, laminarity is incompatible with compression and idea 6 is closed as a theorem.
3. (Optional, ≤ ×1.05) Partial laminarity on the shipped word. **Stop:** if the replayed helper mass falls by < 3%.

## Files and reproduction

Run each from this directory. Runtimes under a 4G cap, one worker.

| Script | Checks | Runtime | Output |
|---|---|---|---|
| `python3 -I -B dfs.py` | exhaustive bowtie certificate and controls (no solver) | 1.9 s | `dfs.out`, `PASS_DFS` |
| `python3 -I -B minbow.py` | minimal directed sub-obstructions (9 corrections) | 0.7 s | `minbow.out`, `PASS_MINBOW` |
| `uv run -q --no-project --python 3.12 --with ortools python -I -B lam.py` | CP-SAT laminar model, replay and abstract controls | 0.9 s | `lam.out`, `PASS_LAM` |
| `python3 -I -B maxlam.py` | greedy direct fraction and disjoint bowtie packing on P120/P112/P104/E7 (writes `maxlam.json`) | 57 s | `maxlam.out`, `PASS_MAXLAM` |
| `uv run -q --no-project --python 3.12 --with numpy --with scipy python -I -B span_milp.py` | joint MILP in the span of the compressing pieces (writes `span_milp.json`) | 21 s | `span_milp.out`, `PASS_SPAN_MILP` |
| `python3 -I -B price.py` | pricing, reads `maxlam.json` (writes `price.json`) | 0.2 s | `price.out`, `PASS_EXPQS_PRICE` |
| `e8.py`, `e8map.json`, `sets.py` | E8 line model (map copied from closures-2) and port systems; P104 is read from the repository's `research/biframe-decoder-ladder/RESULTS.json` (override the root with `IMB_REPO`) | — | — |
| `rho_min.py`, `kappa_pair.py`, `transfer.py` | pricing engines (expL2, expNB, expNE), the closures-2 copies | — | — |

Not included: the growth, ladder and gauge-anneal searches (`grow.py`, `ladder.py`, `lamgauge.py`) and the random census
(`probe.py`, `minobs.py`). Their results are quoted above.

Sources: expQS. Prepared with Claude Code (Claude Opus 5.5) assistance. The bowtie lemma (exhaustive) and the rank lemma are
proved; gauge-family and greedy numbers are search results; pricing values are planning numbers. Not independently reviewed.
