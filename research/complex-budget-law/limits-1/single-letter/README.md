# Single-letter theory (expXB): κ is a log-Cesàro mean of seed exponents; one seed per side is enough; Problem X

Status: theory lane, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`. Multipliers are over
κ = 0.0017.

## Verdict

| Item | Value |
|---|---|
| Question | can level-dependent ("growing") seeds beat a fixed seed, and how much of the gap to κ → 1 is the recursion's fault? |
| Varying-seed law (Theorem 1) | **κ is set by the log-length-weighted (Cesàro) mean of the seed exponents a_j actually traversed, not by their liminf.** Alternating seeds a = 0.05 / 0.40 give a_eff = 0.224 (Cesàro 0.225) |
| Overheads | additive glue h(y) ≤ y^{1−a−ε} per node is free; a constant multiplicative loss K per level costs ln K / ln b_j (fatal at fixed b, amortised to 0 when b_j → ∞); Moore-bound seed growth forces ln b_j ≲ √(C ln y), which still gives a_eff → 1 |
| Single seed suffices | **one finite seed with a ≥ 1/9 on each side (bit and complex) gives κ ≈ a/(1+a) ≥ 0.1 by itself** (×59). Level-scaling only matters for approaching a_∞ and κ → 1 |
| Amplifier | on the Braverman–He seed the repo recursion gives a = ln(157/147)/ln 5 = **0.0409**; BGS gives 0.00298 printed (×13.7), 0.00412 α-optimised (×9.9), 0.00198 all-N (×20.7). **The repo recursion is a 10–20× better coding-gap amplifier** (denominator ×6.97, numerator ×1.97) |
| Bit side, dimension | already met: Hermitian Singer PG(2,64)−M has s = 0.2154 > 0.193 (girth 6). What is open is the schedule (now equivalent to Boot = dim C, `../bootstrap/`) and frames |
| Complex side | the requirement s ≥ 1/(1+(g−1)^{8/9}) exceeds every proved cap for flat codes, and no non-flat construction is known (`../complex-dual/`) |
| Open problem | **Problem X(F, a)** below, needed for F = F_{2^e} (bit) and F = ℂ (complex) |
| Multiplier | **×0** realised. Conditional ceilings in the table below |

## Changed hypothesis

The earlier framework fixed one seed for all levels. Here the seed may depend on the scale x (number of digits). This
separates two questions that were mixed before: *how good can one seed be* (a finite problem) and *how much does the
recursion lose* (nothing beyond the Cesàro averaging and amortised overheads).

## Theorem 1 (log-scale Cesàro law)

Model: T(x) = Σ_r (n_r/W)·T(r x/m) + h(x), T = O(1) for x ≤ 1, seed S(x) with histogram (n_r, m, W). p(x) = 1 − a(x) is
the root of W⁻¹Σ n_r (r/m)^p = 1 (the tilt identity of the budget law). Put Φ(x) = exp ∫₀^{ln x} p(e^ℓ) dℓ. Then

  ln T(x) = ln Φ(x) + E(x),  |E(x)| ≤ Σ_{levels j on the worst chain} L_max,j·ω_j + ln(1 + Σ_j h(x_j)/Φ(x_j)),

where ω_j is the oscillation of p over the level window [ln x_j − L_max,j, ln x_j].

*Proof sketch.* With G = T/Φ and q_r(x) = n_r(r/m)^{p(x)}/W (a probability vector by the tilt identity),
Φ(rx/m)/Φ(x) = (r/m)^{p(x)}·e^{±L_r ω}. Hence G(x) ≤ e^{Lω}Σ_r q_r G(rx/m) + h(x)/Φ(x) and the mirror inequality from below;
iterate down the deepest chain. ∎

**Corollary.** a_eff(u) = Ā(u) − E/ln u with Ā(u) = (ln u)⁻¹∫₀^{ln u} a(ℓ) dℓ. So κ is governed by liminf Ā, which is
≥ liminf a_j and equals lim a_j when it exists. T(u) = O(u^{1−a_∞}) exactly iff Σ_j L_j(a_∞ − a_j) < ∞ (for Hermitian,
a_∞ − a(q) ≈ 0.85 q^{−1/2}, so q_j = 4^j suffices).

**Single-letterisation.** Ψ_{S⊗S′} = Ψ_SΨ_{S′} (expQG), and `../complex-dual/` Theorem E gives the dual counterpart
Θ̂(S⊗S′) ≤ Θ̂(S)Θ̂(S′). Together with Theorem 1: a seedwise bound a(S) ≤ a₀ closes every tensor product and every Cesàro
mixture of such seeds at a₀. **The κ problem is single-letter: it is decided by the best finite seed on each side.**

## Problem X (the remaining open problem, exact statement)

- Γ: finite graph, every degree even (bounded), girth g. F: a field.
- At each vertex v, A_v ⊆ F^{E(v)} of dimension deg(v)/2, **MDS** (every deg(v)/2 coordinates are an information set).
- Tanner code C = {x ∈ F^{E(Γ)} : x|_{E(v)} ∈ A_v ∀v}; s := dim C/|E| (pure excess: naive count gives 0).
- C is **schedulable** if there is an Eulerian orientation and an edge set R, |R| = dim C, with R an information set of C
  and Γ − R acyclic.

> **Problem X(F, a).** Find (Γ, F, A_v) with C schedulable and s ≥ 1/(1 + (g−1)^{1−a}).

For κ ≥ 0.1 (a = 1/9) the thresholds are 0.193 (g = 6), 0.151 (8), 0.124 (10), 0.106 (12). For a = ½ they are
1/(1+√(g−1)). The two instances needed are **X(F_{2^e}, a)** (bit side; the gap must not lift to char 0) and **X(ℂ, a)**
(complex side). `../bootstrap/` shows that schedulability forces dim C = Boot(Γ), and `../complex-dual/` maps what is
proved about X(ℂ).

## Conditional ceilings (all unrealised, ×0 today)

| family | a_∞ | κ_∞ = a/(1+a) | × over 0.0017 |
|---|---|---|---|
| Hermitian Singer PG(2,q)−M, char 2 (bit only) | 0.317 | 0.241 | ×142 |
| same family, odd q (`../suppliers/`, arithmetic supplier in char p) | 0.317 | 0.241 | ×142 |
| growing-girth Tanner at constant s = ¼ | → 1 | outer-capped at ¼ | ×147 |
| generalised polygons (g ≤ 16), s = ¼ | ≤ 0.594 | outer-capped at ¼ | ×147 |
| far-pair expanders (`../seeds-negative/` expRJ) | 0.458 at 10⁵ roles | 0.31 | ×185 |
| growing lattices (`../lattices/`) | → 0 under the template | — | ×0 |

## Next step

None for the theory. The open items are Problem X on each side (see the index README). Stop test for any growing-family
proposal: compute a_j along the family; if it decreases, the family is worth its best member, never more.

## Files and reproduction

| Script | Checks | Runtime | Output |
|---|---|---|---|
| `python3 -I -B xb.py` (stdlib) | A: Cesàro law on 4000-level simulations; B: additive overhead threshold; C: multiplicative overhead; D: amplifier comparison on BH; E: Hermitian a(q) and the Tanner threshold table s_req(a, g) | 2 s | `xb.out`, `PASS_EXPXB` |

Scope: Theorem 1 is an elementary Akra–Bazzi-type argument, not reviewed. The overhead h ≈ ln N_j for on-tape seed
application is an assumption. The amplifier comparison compares guarantees, using expRL's c₁ = 67.94. Problem X's
schedulability clause follows expNI's exact-schedule model; frames realising uniform rank m/b are a further requirement.

Source: expXB. Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.
