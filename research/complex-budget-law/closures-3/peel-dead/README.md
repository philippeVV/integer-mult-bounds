# The peel is dead: no deficit with ≤ W + 2 climbs, nor with any number of independent extra climbs (a = 0)

Status: scoped closure, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017. The contract is upstream `scripts/finite_bit_contract.py`: XOR circuits with gates of any
arity, one common rational frame per gate, endpoint identities T_{ρ(w)} − S_w = I, and s = Σ_e rank(frame(head) − frame(tail)).

## Verdict

| Item | Value |
|---|---|
| Peel with one climb per role route (u = 1, μ = 0) | **dead**: every single climb is exactly +I, so s = Wm and ν = 0 (Theorem B). Realised **×0** |
| Peel with μ ≥ μ₀ > 0 mass per family at rank ≤ r₀ | a ≤ ln(1 + ν/μ₀)/ln(m/r₀) (Lemma P′). Today's ν = 5/14, μ₀ = 1: ≤ 0.071 at m = 73, ≤ 0.022 at m = 10⁶. With the stage nesting (μ ≥ 2): ≤ 0.038 at m = 73. This sits inside the framework's realistic ceiling ≲ 0.03 (×20) |
| One extra climb in the whole circuit | s ≥ Wm (Proposition C) |
| k = 2 extra climbs (any W, m, arity, circuit) | **δ ≤ (5/3)·R**, R = rank E₁ + rank E₂ (Theorem E, expQX) |
| any fixed number k of extra climbs | **δ ≤ C(k)·R** (Theorem F), so ν ≤ C(k)·μ and a ≤ ln(1 + C(k))/ln(m/r₀) → 0: **the peel is dead for any fixed k** |
| **≤ W + 2 climbs (k ≤ 2), any frame, any arity** | **s ≥ Wm** (Theorem G, expRG): upgrades Theorem E to δ ≤ 0, and covers I ∈ span(E₁, E₂). **A strict deficit needs ≥ W + 3 climbs** |
| The δ = R hyperbola | **not realisable** by any circuit: all 200 rectangle-feasible class graphs violate a general entropy cut, and exact F2 CP-SAT says UNSAT for all 200 (positive controls: 892/892 SAT) |
| **Any number k of independent extra climbs** (H: I, E₁..E_k independent over Q) | **n_res ≤ k and δ ≤ 0 for every k** (expRH: Theorem L, locality with window 0; Theorem U, single use). In the level setting the peel has **a = 0 exactly** |
| **Open residue** | **dependent extras** (I ∉ span E, but the E's dependent; Theorem G covers dependent k ≤ 2), and **strict deficits outside the level setting** (every N − W climbs span I). Conjecture T is now proved in its count half under H |

## Changed hypothesis tested

The framework ceiling (`../../framework-ceiling/`) needs 14 hypotheses (table in `hyp_drop.py` §A and expQH). Dropping any one
of them alone gains ×1.0, except four. The top-ranked allowed violation was the **peel**: drop the stage template I3, the data
floor I4 and the halving guard I5 together, so that each family pays its m − u coordinates in **one** child of rank Θ(m).
Then the exact tilt root tends to **a → ν** (×210 fused at today's ν = 5/14, ×588 at ν = 1), and the root decays like
1/ln m in any per-family mass μ at O(1) rank (expQH Lemma P, `hyp_drop.py` §B/§B2). The question was therefore:

> **Q-peel.** Is there a balanced XOR network computing the swap with ν > 0 in which every role's frames take only two values
> along its route (one nonzero climb)?

## Proof sketch

**Theorem A (level cut; any characteristic, any arity, no routing).** Let f be any additive map on Mat_m(Q) into an abelian
group with f(I) ≠ 0. Then at least W physical edges have f(Δ_e) ≠ 0.
- Label each node by λ(v) = f(M_v). A gate is one node, so its ports share one label.
- The sinks at level t are a function of the level-t sources and the edges entering level t. By the endpoint identity, the
  sink of w is at level λ(source w) + f(I) ≠ λ(source w), so every level-t sink carries an input from outside level t.
- With independent uniform input bits, #sinks(t) ≤ H(entering_t) ≤ #entering(t). Summing over t gives the claim.

Equivalent form (A_Z): any (#climbs − W + 1) climbs span I over Z.

**Theorem B (the single-climb peel collapses).** If every route has at most one edge with Δ ≠ 0, that edge's Δ is D_w = I.
- With ≤ W climbs, (A_Z) gives I ∈ Z·D_w, so D_w = I/n_w with n_w a nonzero integer.
- Around a ρ-cycle C the endpoint identities telescope: Σ_{w∈C} D_w = |C|·I, so Σ 1/n_w = |C| with |1/n_w| ≤ 1. Hence n_w = 1.

So the family pays exactly m, ν = 0, and the tilt root is a = 0. This covers every W, m, gate arity and circuit.

**Lemma P′ (sharpened Lemma P).** For one family with paid mass m − ν, of which μ > 0 sits at rank ≤ r₀ and the rest anywhere,
a ≤ ln(1 + ν/μ)/ln(m/r₀). Proof: rank-r mass M contributes M r^{−a} m^{a−1} to Ψ(a), which decreases in r; at the root Ψ = 1.
There is no O(1/m) term, and no placement assumption.

**Proposition C (one extra climb).** If all routes are single-climb except one extra climb E ∉ Q·I, then by (A_Q) every
D_w = c_w(I + t_w E) with distinct nonzero t_w. The deficits are geometric multiplicities of distinct eigenvalues of E, so they
sum to ≤ rank E, and s ≥ Wm.

**Two classes, alternating (A…B…A).** Not excluded by Theorem A, but such a route has c ≥ 2 climbs of equal rank r, so
ℓ̄ ≥ ln c. That is the template regime the framework ceiling already covers.

## Few extra climbs: Theorems E and F (expQX)

Setting: W main climbs plus k designated extra climbs E₁..E_k with I ∉ span_Q(E); δ := Wm − s.

**Step 1 (I-coefficient 1).** Take f linear with f(I) = 1, f(E_i) = 0. (A_Z) on {D, E₁..E_k} gives f(D) = 1/x with x ∈ Z ∖ 0;
telescoping over a ρ-cycle C needs ≥ |C| main climbs on C, and there are exactly W, so every f(D) = 1. If I, E₁..E_k are
Q-independent, D = I + Σ P_i E_i with P ∈ Z^k.

**Step 2 (graded class graph).** Contract the Δ = 0 edges into classes. Main climbs go up exactly one integer level, E-climbs stay
on a level, and each demand goes up exactly one level.

**Lemma R (rectangle cut).** Theorem A's entropy cut works for any node set L. With L = (levels < t) ∪ A′ ∪ B′ and
X(a, b) = #demands a→b − #main climbs a→b, every rectangle satisfies X(A″, B′) ≤ c_A(A″) + c_B(B′) (E-arcs leaving A″, entering
B′). A main climb whose class pair also carries a demand is exactly I, so the number of non-I main climbs per level pair is
≤ ΣX⁺, and
- (i) ΣX⁺ ≤ k_A + k_B if k_A·k_B = 0; (ii) ΣX⁺ ≤ 2 if k_A = k_B = 1; (iii) ΣX⁺ ≤ k_A + k_B + 2k_A·k_B in general.
- Machine check (CP-SAT, exhaustive over E-digraphs with multi-arcs): max = k for k ≤ 2 (`rect_lemma.py`); k = 3 reaches 4
  (`rect_max.py`); k = 4 reaches 5 (scouting, not replayed).

**Kernel facts.** Every nonzero x lies in the kernels of at most k climbs, and ker D ⊆ U := Σ Im E_i (dim u ≤ R). So the
corank sum of any k + 1 main climbs is ≤ k·u.

**Theorem E (k = 2).** At most 4 main climbs are ≠ I (Lemma R), and averaging the triple bound gives Σd ≤ (4/3)·2u, so
**δ ≤ 8u/3 − R ≤ (5/3)R < 2R**. If E₂ ∥ E₁, δ ≤ 0. The case I ∈ span_Q(E₁, E₂) is not covered, but then R ≥ m already.

**Theorem F (bounded k).** n ≤ Σ_t f(k_t) ≤ 2k + k²/2 non-I main climbs, and **δ ≤ (n·k/(k+1) − 1)·R**, e.g. ≤ 6.5R at k = 3.
So ν ≤ C(k)·μ per family, and Lemma P′ gives a ≤ ln(1 + C(k))/ln(m/r₀) → 0 for every fixed k. With extras spread out, only O(k)
routes can be deficient at all; Θ(W) deficient families need k = Θ(W), or extras concentrated on a few levels with k ≳ √W.

**Search (`conic_search.py`).** A_Z-admissible deficient sets without the Lemma R cap: 12,000 random E-pairs (k = 2) and 2,300
E-triples (k = 3) give δ_max = 0 (scouting; a 500-sample k = 2 slice is replayed here). The structured divisor hyperbola
(αβ = c, c ≤ 720) gives at most 4 admissible points, i.e. δ = R exactly. Nothing exceeds 2R.

**Still open.**
- **Unbounded k.** C(k) grows because the kernel bound uses the global U = Σ Im E_i. The missing statement is a **locality
  lemma**: D − I is supported on E's within O(1) levels of the unmatched climb. With it, δ ≤ C·R with C independent of k, and the
  peel is dead outright. Stop test: a k = 4 class graph with an unmatched climb whose P uses a far E, under all cuts. **Now proved under
  independence (H) with window 0 by expRH (Theorems L and U below); open only for dependent extras.**
- **Is the δ = R hyperbola configuration realisable as a circuit?** No (expRG, next section).

## k ≤ 2 closed: Theorem G, Lemma S and the class-graph screen (expRG)

**Lemma S (level split; linear circuits).** In the level setting, level t's subcircuit is a bijection (src_t, in_t) →
(snk_t, out_t) with out_t ∈ span(src_t) and snk_t ∈ span(in_t), #snk_t = #in_t, #src_t = #out_t. (Theorem A is tight at every
level; induct from the lowest level.) So the circuit is a chain of independent two-layer networks.

**Class-graph screen.** Frames are constant on Δ = 0 components, so contract each class to one free gate: **s ≥ m·F_*(class
graph)**, and if the coarse circuit is Q-realisable then s ≥ Wm (upstream's char-0 argument verbatim). A deficit therefore needs a
class graph that is flow-deficient **and** characteristic-dependent. Also F_*(class graph) = W ⇒ n_res ≤ k.

**The hyperbola is not realisable.** It needs n_res = 4 unmatched main climbs with only 2 E-arcs, so each arc would have to serve
both adjacent level pairs. Every class graph that passes the rectangle cuts with n_res = 4 (200 of them, layouts 2,2,2 / 1,2,2,1 /
2,3,2; 0 on 1,2,1) violates a non-rectangle entropy cut (for example L = {0, 2, 4}: 2 sinks fed from outside, 1 entering edge, i.e.
one arc cannot carry both an in-bit and a src-bit), and exact F2 CP-SAT with one free gate per class gives **UNSAT for all 200**.
Positive controls: every all-cut-feasible class graph on the same layouts and arc placements is F2-SAT, **892/892**.

**Theorem G.** Any circuit in the contract with at most W + 2 climbs has s ≥ Wm, for every m, every frame and every arity.
- N ≤ W − 1 is impossible (Theorem A).
- If some designation of k = N − W ≤ 2 extras has I ∉ span_Q(extras): Lemma S; **arc purity** (with ≤ 2 arcs, each arc carries
  downstream information of one type only; case analysis over single, chain, fan-in, fan-out, parallel, antiparallel); a typed
  rectangle inequality and typed Lemma R give n_res ≤ 2; the residual is Eulerian and passes the cut condition, so
  Hu / Rothschild–Whinston (2 commodities) routes it integrally, and telescoping gives s ≥ Wm.
- If every pair of climbs spans I: the climbs commute (span(I, X)), are simultaneously triangular, and scalar Theorem A per
  diagonal coordinate gives s ≥ Wm.
- **The arc-purity case analysis is the delicate step. It is backed by an exhaustive entropy-cut check that does not use purity:**
  max n_res = k for k = 1, 2 on every layout tried (9 layouts for k = 2), and 429,640 cut-feasible configs with k ≤ 4 all route.

**Conjecture T (OPEN).** In the level setting every linear (F_q) circuit has n_res ≤ k; its residual then routes, so s ≥ Wm.
Evidence only: for k = 3, 4 the cut condition alone is not enough (about 65,000 cut-feasible class graphs with n_res > k), but all
64,680 tested are F2-UNSAT (192 with arc 3-cycles untested), and 1,598,148 cut-feasible configs with n_res ≤ k all route. If T
holds, a deficit needs every set of N − W climbs to Q-span I (the regime of upstream's signed-source templates).

**What a strict deficit must look like:** ≥ W + 3 climbs; non-commuting climbs; a class graph that is flow-deficient and not
realisable over Q or any F_p prime to the frame denominators; and (conditionally on T) no level structure. Nothing searched meets
all of these. κ: ×0.

## Every k under independence: Theorems L and U (expRH)

Hypothesis **H**: I, E₁, …, E_k are linearly independent over Q (levels, classes and the general node-set cut as above).

- **Lemma 1 (level balance).** dem(t → t+1) = D(t → t+1) for every t.
- **Lemma 2 (component balance).** X(K, K′) = 0 for every E-component K at level t and K′ at level t+1 (rectangle cuts over unions
  of components have bound 0, and the terms sum to 0 by Lemma 1).
- **Lemma 3 (forest + hub).** Under H each level's E-graph is a forest (a cycle gives Σ ±E_i = 0), and all demands from K to K′
  use one class pair (a₀, b₀), the hub: M_{b₁} − M_{b₂} = M_{a₁} − M_{a₂} lies in span E(K′) ∩ span E(K) = 0.
- **Theorem L (locality, window 0).** A main climb D: a → b with a ∈ K, b ∈ K′ has D − I = π_K(a → a₀) + π_{K′}(b₀ → b), the signed
  E-sums along the tree paths to the hub, with coefficients in {0, ±1}. D − I uses only E's in D's own two components (levels t
  and t + 1), never t − 1 or t + 2. (The support statement needs only I ∉ span E.)
- **Theorem U (single use).** Each E-arc lies on the path of at most one main climb in the whole circuit. Cases U1 and L1 follow
  from rectangle cuts. Cases U2 and L2 need a **three-level entropy cut** L = (levels ≤ t − 2) ∪ (level t − 1 ∖ 𝒦₁) ∪
  (level t ∖ X) ∪ 𝒦₂, into which the only entering E-arc is e; cancelling with Lemma 2 gives Σ D(K₁ → X) + Σ D(Y → K₂) ≤ 1. The
  three-level cut is necessary: with rectangle cuts alone, CP-SAT finds an E used twice.
  **Review this step first:** the sink/entering-climb accounting of the three-level cut is the least checked part of the proof.
- **Corollary.** ker D embeds into ⊕_{e ∈ path(D)} Im E_e and the paths are disjoint, so Σ_D (m − rank D) ≤ Σ rank E_i:
  **s ≥ Wm, δ ≤ 0, and n_res ≤ k, for every k, W, m and arity under H.** This proves the count half of Conjecture T under H, and
  explains expRG's k = 3, 4 "cut-feasible, n_res > k" graphs: they violate the hub rule, i.e. they are not frame-consistent.

**Coverage of the checks (honest).** Brute force over every node-set cut with explicit potentials: 3 levels × 2 classes, k ≤ 4,
multiplicities ≤ 2 (34,696 cut-feasible graphs; C0–C4 hold, max per-E use 1, random diagonal plantings give max δ/R = 0), and
2 levels × 3 classes at k ≤ 1. CP-SAT: up to 5 levels × 3 classes, k ≤ 5, **mostly with rectangle cuts only**; the stop test
("some deficient climb uses an E outside levels t − 1..t + 2", k = 4) is INFEASIBLE, and the negative control without cuts finds a
far use in 0.1 s. The single-use check with the three-level cut ran for T = 4, p = 3, k ≤ 4 (k = 4 in 570 s; k = 5 stopped).
Dependent labels (repeated equal matrices, T = 3, p = 2, k ≤ 3): hubs can split and a label can be used up to 3 times, but never
more often than it has arcs, and max δ = 0 (not replayed).

**Still open:** dependent extras (Lemma 3 and Theorem U lose their proofs; Theorem G covers k ≤ 2), and strict deficits outside the
level setting. The latter is the framework's real remaining question, not a peel question.

## Routability of small gadgets (expQY; scoped, no κ impact)

expQQ conjectured that every two-port GL(2,F2) circuit computing a permutation matrix is token-routable (ρ ∈ S(C)), which would
give s ≥ Wm for every frame assignment on two-port gadgets.

- **Lemma 1 (reduction to CNOT form).** Every g ∈ GL(2,F2) is uniquely τ^x·E with E ∈ {1, E_ab, E_ba}. Relabelling by any routing
  (the token frame) gives C = P_{π_n}·C̃_x and S(C) = π_n ∘ S(C̃_x); choosing x_k = type(g_k) leaves only transvections. So the
  conjecture is equivalent to its CNOT-only case.
- **Theorem 2 (Bruhat, any field).** If every gate is a GL(2,F) gate on an adjacent pair of one fixed wire order and the product is
  monomial D·P_ρ, then ρ ∈ S(C). Proof: GL(2)_{i,i+1} ⊂ B ∪ Bs_iB, and the BN-pair axiom tracks the reachable Bruhat cells exactly
  as token routing does; Bruhat cells are disjoint.
- **Corollary 3.** If some routing and token order make every coding gate meet adjacent tokens, the conjecture holds for that
  circuit.
- **The conjecture is still OPEN.** No counterexample: the C2-prefiltered exact MITM hunt found 0 realisable non-routable circuits
  among 34,768 candidates (W = 5, n = 8–12; W = 6, n = 10–12); 1.11 M matrix-preserving rewrite moves (W = 4–7, n ≤ 20) found 0
  non-routable; the GL(2,F3) hunt (W = 4, n = 6) found 0. Three-port GL(3,F2) scans (W = 4, n ≤ 4 and mixed n ≤ 5; W = 5, n ≤ 5;
  W = 6, n = 4) found **no C2-passing non-routable ρ at all**.
- **The two C2-passing non-routable instances** (W = 4, pairs (02)(12)(13)(03)(02)(01), ρ = (0,3,1,2) or (1,3,0,2)) show that the
  LGV/C2 condition is necessary but not sufficient; **both are F2-unrealisable**.
- κ consequence: none (×0). If the conjecture holds, deficits need gates with ≥ 3 ports.

## Searches (consistent with the theorems; not needed for them)

- Two-port GL(2,F2) gadgets, every realisable all-role permutation, exact union-find frame solving under the two-class
  constraint: s = Wm in every case (W = 3 exhaustively to n = 5 here; the scouting run also did W = 4 to n = 6).
- Upstream's §4 meet-in-the-middle on the balanced 10-role Fano graph: 279,936 / 279,936 halves, 0 words with
  ρ(0, 1, 2) = (7, 0, 3), and 14,064 realisable permutations, **all token-routable**. So upstream §1/§2 already give s ≥ Wm for
  every frame assignment on that graph.
- Not included (quoted from scouting): the W = 4 enumeration (`enum_W4.log`), random MITM graphs, three-port scans and the
  CP-SAT check on upstream compiled circuits (stopped without an optimality proof). None found a non-routable permutation.

## Multiplier

- u = 1 peel: **×0** (ν = 0 exactly).
- μ ≥ 1 per family: ≤ ×39 unfused before any other charge (a ≤ 0.071 at m = 73), and it falls like 1/ln m. That is already
  inside the framework ceiling.
- Any fixed number k of extra climbs (μ → 0): a ≤ ln(1 + C(k))/ln(m/r₀) → 0 (Theorem F); k ≤ 2: no deficit at all (Theorem G); any k under H: δ ≤ 0, a = 0 (Theorems L, U).
- Realised: ×1.00.

## Files and reproduction

Run each from this directory. The runtimes are under a 4G memory cap, one worker.

| Script | Checks | Runtime | Output |
|---|---|---|---|
| `python3 -I -B hyp_drop.py` | expQH: the 14-hypothesis table (§A), peel roots (§B), Lemma P on 284 cases (asserted, §B2), stage nesting (§C), width chain (§D), outer caps (§E) | 3.4 s | `hyp_drop.out`, `PASS_EXPQH_HYPDROP` |
| `uv run -q --no-project --python 3.12 --with mpmath python -I -B lemma_p.py` | Lemma P and P′ at 60 digits: 284 grid cases and 3,000 random placements, 0 violations | 51 s | `lemma_p.out`, `PASS_EXPQQ_LEMMAP` |
| `uv run -q --no-project --python 3.12 --with numpy --with scipy python -I -B qq_enum.py 3 5` | exhaustive two-port gadgets W ≤ 3, n ≤ 5: Theorem A at m = 1, monotone and arbitrary two-class bounds = 1 | 12 s | `qq_enum.out`, `PASS_EXPQQ_ENUM` |
| `python3 -I -B qq_mitm.py fano10` | upstream §4 Fano MITM: 14,064 realisable ρ, 0 non-routable, 0 hits on (7, 0, 3) | 79 s | `qq_mitm.out`, `PASS_EXPQQ_MITM` |
| `qq_core.py` | circuits, union-find with integer offsets, two-class gluing solver, LP and scalar bounds (library) | — | — |
| `uv run -q --no-project --python 3.12 --with ortools python -I -B rect_lemma.py` | Lemma R for k ≤ 2: max ΣX⁺ = k on every E-digraph, p, q ≤ 4 (CP-SAT, one worker) | 8.8 s | `rect_lemma.out`, `PASS_EXPQX_RECT` |
| `uv run -q --no-project --python 3.12 --with ortools python -I -B rect_max.py` | the maxima by k: 0, 1, 2, 4 for k = 0..3 | 30 s | `rect_max.out`, `PASS_EXPQX_RECTMAX` |
| `python3 -I -B conic_search.py 2 2 500 1 2 rank1 12 4` | A_Z-exact search, k = 2, m = 2, 500 rank-1 pairs, P-box ±12: δ_max = 0 | 2.0 s | `conic_k2.out`, `PASS_EXPQX_CONIC` |
| `uv run -q --no-project --python 3.12 --with ortools python -I -B hyper_cg.py L 1000` for L = `2,2,2`, `1,2,2,1`, `2,3,2` | expRG: every rectangle-feasible n_res = 4 class graph with 2 arcs is F2-UNSAT (4 + 4 + 192 = 200) | 0.9 s, 0.8 s, 1.2 s | `hyper_cg_{222,1221,232}.out`, `PASS_EXPRG_HYPER` |
| `uv run -q --no-project --python 3.12 --with ortools python -I -B ctl_realise.py 2,2,2 '2-3;3-2' 192` | positive control: 192/192 all-cut-feasible class graphs are F2-SAT (the other 700 controls run the same way) | 1.0 s | `ctl_realise.out`, `PASS_EXPRG_CTL` |
| `uv run -q --no-project --python 3.12 --with ortools python -I -B route_lek.py 4 2,2,2` | Conjecture T evidence: 19,776 cut-feasible configs with n_res ≤ 4 all route | 3.7 s | `route_lek.out`, `PASS_EXPRG_ROUTE` |
| `rg_core.py`, `max_nres.py`, `route.py` | class graphs, entropy cuts, exact F2/F_p realisability (CP-SAT, one worker); max n_res model; residual routing (libraries) | — | — |
| `uv run -q --no-project --python 3.12 --with ortools python -I -B classgraph_sat.py --claim far` | expRH stop test: no deficient climb uses a far E (T = 5, p = 3, k = 4, rectangle cuts, one worker): INFEASIBLE | 10 s | `rh_stop_far.out`, `PASS_EXPRH_FAR_INFEASIBLE` |
| `uv run -q --no-project --python 3.12 --with ortools python -I -B classgraph_sat.py --claim far --cuts none --expect FEASIBLE` | negative control: without cuts a far use exists | 0.6 s | `rh_negctl.out`, `PASS_EXPRH_FAR_FEASIBLE` |
| `uv run -q --no-project --python 3.12 --with numpy python -I -B rh_brute_check.py 3 2 4 2 2 3` | independent brute force (`brute_check.py` of expRH), every node-set cut: C0–C4, single use, δ ≤ 0 | 30 s | `rh_brute.out`, `PASS_EXPRH_BRUTE` |
| `python3 -I -B qy_check.py all` | expQY: Lemma 1 on 300 random circuits, Theorem 2 over F2 and F3 (exact MITM), the C2 gap example, a small hunt | 15 s | `qy_check.out`, `PASS_EXPQY_CHECK` (also prints `ALL OK`) |

Sources: expQQ, expQH, expQX, expQY, expRG and expRH, copied with a PASS marker (and, where needed, a here-relative import line) added.
`rect_max.py`: a loop-order slip in the scouting copy (the size filter read `ktot` before the `ktot` loop) is fixed here.

Prepared with Claude Code (Claude Opus 5.5) assistance. Theorems A, B, E, F, G, L, U, Lemmas R, S, 2 and 3, Proposition C and the expQY Lemma 1 /
Theorem 2 are elementary but have not been
independently reviewed.
