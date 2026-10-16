# Bootstrap (expXE, expXL): a field-free ceiling, schedulable ⇒ dim C = Boot, and the graph-exclusion table

Status: scoped theorem plus exclusion table, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017.

## Verdict

| Item | Value |
|---|---|
| Bootstrap theorem (any field) | for every graph, every field and MDS local codes of dimension k_v: **dim C ≤ Boot(Γ, k) := min over vertex orders σ of Σ_v (k_v − j_v^σ)^+**, j_v^σ = #edges to earlier vertices. Any order gives a valid bound |
| Schedulable ⇒ extremal | if C is schedulable (Problem X's sense), then **dim C = |R| = Boot(Γ)** and the topological order of Γ − R is an optimal bootstrap order. Problem X splits into a purely combinatorial condition Boot(Γ)/E > 1/(1+(g−1)^{1−a}) and codes that are exactly bootstrap-extremal |
| Tightness | complex: PG(2,3) 8 = Boot, biaffine(4) 10 = Boot (10/64 = 0.156 beats expRF's 1/8, still < 2/11). char 3: Braverman–He's PG(2,9) code has 166 = Boot (`../suppliers/xj`). char 2: PG(2,16) 900 vs ≤ 914, PG(2,8) 114 vs ≤ 117, PG(2,16)−M 772 vs ≤ 784 |
| Exclusions | X(·, 1/9) (κ ≥ 0.1) is impossible **in every characteristic** on PG(2,q) for q ≤ 13, PG(2,4)−M, PG(2,16)−M, W(q) for q ≤ 7, W(4)−M, and every degree-4 graph tested. Not excluded: PG(2,17) (0.198), PG(2,27), PG(2,32)−M, PG(2,64)−M (0.221), W(11) (0.168) |
| Correction in this replay | W(7): expXE's greedy value 483 gives Boot/E = 0.1509, which is *above* the exact threshold 0.15063 (the table rounded it to 0.151). With 300 restarts the greedy reaches **479 (0.1497)**, so the "excluded" verdict stands, with a margin of 3 sessions. Also W(5) ≤ 126 and W(4)−M ≤ 29 |
| Limits of the method | Boot ≥ cutwidth/2, so Boot/E stays bounded below on expanders; the theorem **cannot close high-girth expander families**, and it cannot see the characteristic (char-0 vs char-2 needs other tools) |
| Multiplier | **×0** realised. Field-free ceiling where not excluded: s ≤ Boot/E ≤ ~¼ at g = 6 (a ≤ 0.317, ×142); PG(2,64)−M Boot/E ≤ 0.221 (κ ≤ 0.178, ×105) |

## Changed hypothesis

Ceilings for Problem X were previously spectral (T1–T3: positivity, characteristic 0 only). The bootstrap bound is a rank
argument that uses MDS only through information sets, so it holds over every field, and it equals the boundary of the
gauge-invariant dual of `../complex-dual/` (expXG Theorem C).

## Theorem and proof

**Theorem.** Γ a (multi)graph, F any field, A_v ⊆ F^{E(v)} MDS of dimension k_v, C the Tanner code. For every order σ of V,
dim C ≤ Σ_v (k_v − j_v)^+.

*Proof.* Let w_e ∈ C* be the coordinate functional x ↦ x_e; these span C*. Put V_v = span{w_e : e ∋ v}, dim V_v = ρ ≤ k_v.
For J ⊆ E(v), |J| = j, MDS gives dim span{w_e : e ∈ J} = dim C|_J ≥ ρ − (k_v − min(j, k_v)). Processing vertices in order σ,
with U the span of earlier V_u (which contains the j_v back-edge functionals): dim(U + V_v) − dim U ≤ (k_v − j_v)^+. Sum. ∎

Remarks: for half-rate d-regular graphs the cost is ½Σ|j_v − k_v|. Without MDS the bound fails (coordinate codes on a
half-factor give s = ½).

**Corollary (schedulable ⇒ extremal).** Take an Eulerian orientation and R with Γ − R acyclic; let σ be a topological order of
Γ − R. Then |R| ≥ #back arcs = Σ_v(k_v − in_earlier(v)) ≥ Σ_v(k_v − min(j_v, k_v)) = cost(σ) ≥ Boot ≥ dim C = |R|. ∎
Conversely, an order with cost(σ) = dim C gives a schedule iff the back-arc transportation problem is feasible (each vertex
with j_v > k_v sends j_v − k_v back arcs to earlier neighbours with deficit, which receive exactly k_u − j_u); by expNC, R is
then automatically an information set. The greedy optimal orders on PG(2,3) and PG(2,9) are **not** flow-feasible (2/8 and
92/166), so a search over optimal orders is needed.

**Consequence for the bit side.** The Hermitian code on PG(2,16)−M (dim 772) is schedulable only if Boot(PG(2,16)−M) = 772;
today Boot ≤ 784 (expNI's schedule, through the corollary).

## Exclusion table (greedy upper bounds on Boot; `xe/check_boot.py`)

Thresholds for X(·, 1/9): 1/(1+(g−1)^{8/9}) = 0.1930 (g = 6), 0.1506 (8), 0.124 (10).

| graph | g | E | Boot (this replay) | Boot/E | verdict for κ ≥ 0.1 |
|---|---|---|---|---|---|
| PG(2,3) | 6 | 52 | 8 | 0.154 | excluded |
| biaffine(4) | 6 | 64 | 10 | 0.156 | excluded |
| PG(2,4)−M | 6 | 84 | 10 | 0.119 | excluded |
| PG(2,5) | 6 | 186 | 31 | 0.167 | excluded |
| PG(2,7) | 6 | 456 | 79 | 0.173 | excluded |
| PG(2,9) | 6 | 910 | 166 (= BH dim, so exact) | 0.182 | excluded |
| PG(2,11) | 6 | 1596 | 298 | 0.187 | excluded |
| PG(2,13) | 6 | 2562 | 489 | 0.191 | excluded |
| PG(2,16)−M | 6 | 4368 | ≤ 784 (expNI) | ≤ 0.1795 | excluded |
| PG(2,17) | 6 | 5526 | 1094 | 0.198 | not excluded |
| PG(2,27), PG(2,32)−M, PG(2,64)−M | 6 | — | 4408, 7046, 58847 (expXE) | 0.208, 0.208, 0.221 | not excluded |
| W(3) | 8 | 160 | 20 | 0.125 | excluded |
| W(4)−M | 8 | 340 | 29 | 0.085 | excluded |
| W(5) | 8 | 936 | 126 | 0.135 | excluded |
| W(7) | 8 | 3200 | 479 | 0.1497 | excluded (margin 3) |
| W(11) | 8 | 17568 | 2950 | 0.168 | not excluded |
| random bipartite d = 4 / 6 / 8 / 12 (n = 1000) | small | — | — | 0.082 / 0.114 / 0.131 / 0.153 | excluded unless g ≥ 18 / 12 / 10 / 8 |

## expXL: Boot vs 772 on PG(2,16)−M — undecided, Boot ∈ [772, 784]

- **Lower side.** 772 comes only from dim C of the Hermitian code. Linear-ordering/feedback LPs give 0, as does any
  relaxation that is consistent on each neighbourhood: x ≡ ½ satisfies every dicycle and triangle inequality. The
  spectral bound gives Boot ≥ 751 (λ₂ = 11).
- **Cutwidth cannot decide it.** An order with max prefix cut exactly 1544 = 2·772 exists (`xl/o3_52.txt`, cost 791).
- **Exact identity (new).** For any order σ and prefix S,
  cost(σ) = ½∂S + Σ_{v∈S}(j_v^S − 8)⁺ + Σ_{v∉S}(later-internal_v − 8)⁺, i.e. cost = ½ × the total ascent of the cut
  profile. So Boot = min_S [½∂S + P(G[S]) + P(G[V∖S])], where P(H) = 0 iff H is 8-degenerate, and expQD's cut + viol
  bound is an equality.
- **Upper side.** About 10¹⁰ annealing moves never beat the Singer cyclic order: 784, which passes the back-arc flow
  check (784/784) and is therefore a valid schedule. Other bests: 787 (subset search), 789 (Baer seed), 791 (cut ≤ 1544).
- **PG(2,64)−M₆₁₇.** Boot ≤ 57599 (Baer-interleaved order plus annealing), down from 57657 cyclic. So Boot/E ∈
  [0.2154, 0.2163], provided dim C is the same for this matching.
- **Verdict.** Open. A certificate needs min-bisection raised from about 1505 to about 1542 together with 9-core
  control on both halves (days of SDP branch-and-bound). 784 stays the operating value.

## Next step and stop test

1. Bit side (expXL): search vertex orders of PG(2,16)−M for cost exactly 772 with a feasible back-arc flow. **Stop** if a
   certified lower bound shows Boot > 772. Repeat at PG(2,64)−M against 0.2154·E.
2. Complex side: only the not-excluded graphs (PG(2,q), q ≥ 17; W(q), q ≥ 11) can host a solution; dense numerics cannot
   reach them, see `../seeds-negative/` (expXM) for the structured attempt on PG(2,17).

## Files and reproduction

| Script (run from `xe/`) | Checks | Runtime | Output |
|---|---|---|---|
| `python3 -I -B check_boot.py` (stdlib; builds PG(2,q), q = 7, 9, 11, 13, 17) | greedy Boot ≤ expXE value on 14 graphs, and the exclusion verdict of each against the exact threshold | 1 s | `check_boot.out`, `PASS_LIMITS_BOOT` |
| `env GF=graphs_small.json uv run -q --no-project --python 3.12 --with networkx python -I -B boot3.py PG3 100 8` | tie-broken greedy reaches 8 on PG(2,3); its back-arc flow is infeasible (2 of 8) | 1 s | `boot3_pg3.out`, `PASS_EXPXE_BOOT3` |
| `python3 -I -B bootrand.py d 1000 10`, d = 4, 12 | random bipartite Boot/E ratios | < 1 s | `PASS_EXPXE_BOOTRAND` |
| `python3 -I -B check_order.py graph16_59.txt ord784.txt` (in `xl/`, stdlib) | cost 784 = ½·ascent, max prefix cut 1568, back-arc flow 784/784 | < 10 s | `xl/check_784.out` |
| `python3 -I -B check_order.py graph16_59.txt o3_52.txt noflow` | cost 791 with max prefix cut 1544 = 2·772 | < 5 s | `xl/check_1544.out` |

Inputs: `graphs_small.json` (expXG), `graphs_gq.json` (W3, W5, W7) and `graphs_gq2.json` (W11) from expXE `gqp.py`
(included). `boot2.py` is the bucket-queue greedy used by all checks. Not included: expXE's complex Grassmannian optimiser
(`opt.py`) and its logs (the r = 8 PG(2,3) solution is replayed in `../complex-dual/xg/`), and the 11 MB Singer graph file
(PG(2,q) with q ≥ 16 and the −M variants are not rebuilt here; their values are quoted from expXE).

Sources: expXE, expXL. Prepared with Claude Code (Claude Opus 5.5) assistance. The theorem is elementary but has not
been independently reviewed.
