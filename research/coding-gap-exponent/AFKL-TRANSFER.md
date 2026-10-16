# κ for oblivious multipliers is capped by the undirected coding-gap exponent

This is local only: nothing was pushed or posted, and there is no κ claim. The base is main `0fdd7d82`. Multipliers are over κ = 0.0017.

The papers (not redistributed here) were read as data. AFKL is arXiv 1902.10935v1, BGS is arXiv 1608.06545v1, and the two October 2026 refutations are arXiv 2610.10108 and 2610.09367. The only numerics are in `exponents.py` → `exponents.out` (standard library, under 1 s: `python3 exponents.py`).

## 0. Verdict

1. **Neither stop test fires.**
   - NCC enters AFKL's multiplication theorem *linearly*, through a single "flow rate ≥ coding rate" step.
   - No upstream step is inherently data-dependent.
2. **The quantitative transfer holds (Theorem A).** Write γ(N) for the undirected coding/flow gap on graphs with ≤ N edges. Then every bounded-degree circuit for the n-bit shift (and hence for multiplication) with s wires satisfies s ≥ n·log_{2c} n / (4·γ(s)).
3. **Corollary B (oblivious TMs).** κ_obl ≤ η_* := liminf_N log γ(N) / log log N.
   - Only "⇒" holds: κ → 1 *requires* gaps of (lg N)^{1−o(1)}.
   - The converse claimed in an earlier local scouting note (κ → 1 "⇔" the gap is saturated) is not implied, because a large coding gap does not yield a multiplier.
   - Since γ(N) = O(log N) (flow–cut), η_* ≤ 1. **No known theorem constrains κ below 1.**
4. **The supervisor term is not in the multiplication theorem.**
   - It belongs only to AFKL Thm 3 (depth-3 circuits / Valiant). There it survives with a 1/g-linear degradation (§1.3).
   - The far-pair term is purely graph-theoretic and does not depend on g.
5. **Obliviousness audit (§2).** Upstream's multiplier (§§02–08) moves its heads as a function of n and the fixed parameters only, with two exceptions, both removable at O(1) cost:
   - (a) sign-dependent "complement if negative" passes;
   - (b) variable-length signed-magnitude records at the interfaces to the cited Harvey–van der Hoeven routines.

   §09 (division/√) and §10 (transposition) are data-dependent, but they are corollaries and lie outside the multiplier. **As written the multiplier is not literally oblivious. A padded variant is oblivious at ×O(1) cost,** modulo an unaudited assumption that HvdH's internal routines are oblivious on fixed-width inputs.
6. **Refutation check (§3): confirmed.** Both papers exist, are dated 7 Oct 2026, and refute Li–Li.
   - Braverman–He: an F₉-linear code on 182 vertices, coding rate 1 against flow ≤ 147/157, giving gap Ω(log^ε n) via BGS.
   - Zhang–Li–Li: subcubic binary-linear networks with routing ≤ K_β/(log n)^c, formalized in Lean.
   - **Neither states a numerical exponent.** BGS amplification of the BH seed guarantees c₂ = ln(1 + ε/2)/ln(16c₁²) ≤ **0.0080** (c₁ ≥ 2, ε = 10/147). That is only 1.6–4.7× today's κ = 0.0017.
   - ZLL's c ≤ 0.035·ln(1/ρ₀) with ρ₀ < 1 unspecified, and presumably tiny.
7. **Consequence.** If upstream's padded multiplier is oblivious, it already *is* an undirected coding-gap family of exponent ≥ κ = 0.0017. That is the same order as the best amplified exponents known. κ = 0.017 (×10) would give a gap exponent above every one that is currently explicit.

**Multiplier: ×1.** This is a necessary condition, not a construction.

## 1. Quantitative AFKL

### 1.1 Definitions

- An instance is an undirected graph Ḡ with unit capacities and k source–sink pairs.
- NC(Ḡ) is the coding rate. In AFKL's entropic definition, an orientation of Ḡ as a DAG carrying H(A_e) ≤ c_e bits is enough. The undirected coding rate is at least that.
- MCF(Ḡ) is the max concurrent flow rate.
- γ(N) := sup{NC/MCF : |E(Ḡ)| ≤ N}. It is nondecreasing in N, since isolated edges can be appended.
- Li–Li (NCC) is the statement γ ≡ 1.
- Write γ_c(N) for the same supremum over max-degree-≤ 2c graphs arising from bounded fan-in/fan-out circuits. Everything below holds with γ_c ≤ γ.

### 1.2 Theorem A (AFKL Thm 2, made quantitative)

**Statement.** Let C be a boolean circuit with arbitrary gates, in-degree and out-degree ≤ c, and s wires, computing the shift s(x, ℓ) on n bits. Then

  s ≥ n·log_{2c} n / (4·γ_c(s)).

Because shifting x by ℓ is multiplication by 2^ℓ (AFKL p. 3), the same holds for n-bit multiplication. Here C is the multiplier with its second operand hardwired, as AFKL do.

**Proof.** This is AFKL §4 verbatim, with the conjecture step replaced.

1. *Far pairs.* In Ĉ (max degree ≤ 2c), at most √n nodes lie within distance ½·log_{2c} n of x_j. Averaging over ℓ gives some ℓ₀ for which at least n − 2√n ≥ n/2 pairs (x_j, y_{j+ℓ₀−1}) are far. No g enters here.
2. *Coding.* Fix ℓ₀. C with unit capacities is a coding solution of rate 1 for the n pairs (s_j = x_j, t_j = y_{j+ℓ₀−1}).
3. *Replaced step.* Instead of "MCF = NC = 1" we have MCF ≥ NC/γ_c(s) ≥ 1/γ. So a flow sends 1/γ units per commodity with Σ_j f^j(e) ≤ 1 per edge (normalized flows f^j of value 1, scaled by r = 1/γ).
4. *Length accounting.* Each far commodity's unit flow has Σ_e f^j(e) ≥ ½·log_{2c} n. Hence
   |E| ≥ r·Σ_{far j} Σ_e f^j(e) ≥ (1/γ)·(n/2)·(½·log_{2c} n). ∎

**Shape of the dependence.** NCC enters exactly once, as a linear factor 1/γ. The proof also yields the *certificate* form, which is the one ZLL use: MCF ≤ |E| / (#far · dist). So a circuit with s ≪ n·log n wires is itself an explicit coding-gap instance, with

  gap ≥ (n/2)(½·log_{2c} n)/s.

### 1.3 The supervisor term (AFKL Thm 3, depth 3 only)

**Where it lives.** The supervisor u appears only in the depth-3 / Valiant theorem (§5). It is not in Thm 1/2.

**What changes.** Rerun AFKL (1)–(4) at flow rate k/g:
- the flow leaving u through the u·a_ℓ and u·t_ℓ edges is ≤ (g/k)·2Σc_ℓ;
- Lemma 4 gives Σc_ℓ ≤ 3n/k + O((n/k)·√(ε lg k)) + (ε/8)·n·lg(2/ε);
- the total demand is n/k.

So u carries a bounded fraction of the flow provided **k ≳ 20g and ε·lg(1/ε) ≲ 1/k ≈ 1/g**. That is, ε = Θ(1/(g·lg g)), against a constant in AFKL. Equation (4) then gives c·log c ≳ lg n / g.

**Conclusion for Thm 3.** It also degrades gracefully: the hypothesis tightens by ~1/(g lg g) and the conclusion loses a factor 1/g. None of this is needed for multiplication or TMs.

### 1.4 Corollary B (oblivious TMs)

**Statement.** Let M be a fixed-tape, fixed-alphabet *oblivious* TM, meaning its head positions at each step depend only on n. Suppose M multiplies n-bit numbers in time T(n). Then

  T(n) ≥ n·lg n / (C_M·γ_{c_M}(C_M·T(n))),

with constants C_M and c_M depending on M only.

**Exponent form.** Hence if T = O(n·lg^{1−κ} n), then γ(N) = Ω((lg N)^κ) for all large N. This uses N_n = Θ(T(n)), N_{n+1}/N_n → 1, and monotone γ. So

  **κ_obl ≤ η_* := liminf_{N→∞} log γ(N) / log log N ≤ 1.**

**Why it holds.**
- The standard simulation turns an oblivious TM into a circuit:
  - one O(1)-size gate block per step, reading the state and the k scanned cells;
  - each cell value is wired to the next visit of that cell;
  - so fan-in and fan-out are ≤ O(k·log|Σ|·|Q|) and there are O(T) wires.
- The upper end is the flow–cut gap: NC ≤ sparsity in undirected graphs, and MCF ≥ sparsity/O(log k) (LLR / Aumann–Rabani). This is the "at most O(log n)" quoted in BH §1.

**Scope.**
- For a *general* TM, Pippenger–Fischer gives circuits of size O(T lg T). The bound then reads T·lg T ≳ n·lg n/γ, which is satisfied by T = n. **So the transfer says nothing about non-oblivious machines.**
- Even under NCC, the transfer caps an AFKL-type bound at Ω(n) for oblivious TMs once γ = Θ(log), which is consistent with κ = 1.

## 2. Obliviousness audit of upstream

Upstream is at `upstream/build/sections`. Only the steps where head motion could depend on data values are listed. Everything else (scans, digit moves, controlled shifts, splits and merges, counters) is driven by addresses and lengths that are functions of n and fixed parameters. `02-streams.tex:61-112`: "controls" are address fields.

| step | where | data-dependent? | cost to make oblivious |
|---|---|---|---|
| input parse x#y, counter for n | 08:52-58 | no (lengths only; promise \|x\| = \|y\|) | — |
| fixed-radix scans, controlled shifts, rational unit scalings, counters | 02 Lemmas 2.1–2.2 | no: offsets depend on address control fields | — |
| fixed-width signed add/sub with carries "using finite control" | 02:198-208 | no: carry is state, the head sweep is fixed | — |
| XOR interchange, parking and descriptor stacks, padding | 04:205-240, 04:356-385 | no: stack contents have n-determined lengths | — |
| packed selected-bit changes and exceptional repair | 05:181-188 ("correction depends only on addresses"), 05:395-410, 08:960-980 (bad fraction "independently of the payload"; stable bit sort on address keys) | no | — |
| row padding, role splits, recursion schedule | 05:555-590 | no | — |
| sign/magnitude ↔ two's complement, "complement and increment if negative, otherwise retain" | 06:522-540, 08:757-763 | **yes**: whether a pass is made depends on the sign | always make the pass, writing either the complemented or the unchanged bits: ×O(1) |
| removal of leading sign-extension bits *to a fixed width* w_c | 08:419-424, 08:485-488 | no: the output width is fixed | — |
| adapters to cited HvdH Gaussian/exponential routines: "ordinary binary magnitude… separator", "emits the signed magnitude after removing leading zeroes" | 07:202-240 | **yes**: field lengths depend on the values, and the resampled vectors are data | pad every field to p + 1 digits plus sign: O(p) per field, which is what the current bound already pays. **Assumption (unaudited):** HvdH's internal routines and their O(n log n) multiplier move obliviously on fixed-width zero-padded inputs. Their algorithms are length-driven (label-and-sort keys are indices), but upstream cites them, it does not reprove them |
| calls to the established multiplier on fixed-width magnitude blocks | 06:490-500, 08:757 | widths are fixed by upstream. §10:20-24 says the cited interface uses "standard binary", i.e. stripped | feed zero-padded fixed width: ×O(1). Same assumption |
| §09 division/√: leading-bit normalization (h, e depend on b, a), zero-input branch, ±1 corrections | 09:12-17, 09:144-150 | **yes** | outside the multiplier. An oblivious normalization shift is a shift problem, so it costs O(M(n)). Not needed for E1 |
| §10 transposition: strip leading zeros of the product | 10:26-31 | **yes** | outside the multiplier (output formatting). Theorem 1's machine itself outputs exactly 2n bits (08:1003-1006) |

**Conclusion.**
- No inherently non-oblivious step exists.
- Upstream's multiplier *as written* is not oblivious, because of the two "yes" rows inside §§06–08. A padded variant M_obl with the same O(n·lg^{1−κ} n) bound is oblivious under the HvdH assumption.
- This matches `community-reconciliation-fourth100-review.md:47,127`, where #634's oblivious-compilation claim is "not established". The audit gives the concrete list that an establishment would need.
- `find_overlap.py 'oblivious coding gap AFKL'` gives 0 hits.

## 3. The October 2026 refutations

| paper | what it shows | exponent |
|---|---|---|
| Braverman–He, *The Multiple Unicast Conjecture is False*, arXiv 2610.10108 (7 Oct 2026) | A deterministic F₉-linear code on a 182-vertex bipartite subgraph of the PG(2,9) incidence graph. It supports 157 unicast sessions at rate ≥ 1, while every fractional MCF has rate ≤ 147/157 (Thm 1.1). Cor 1.2: gap Ω(log^ε n), ε > 0, via BGS. It notes an O(log n) upper bound by flow–cut. | Not stated. BGS guarantees ε ≤ ln(1 + 5/147)/ln(16c₁²) ≤ **0.0080**: 0.0080 / 0.0067 / 0.0056 / 0.0045 for c₁ = 2 / 3 / 5 / 10 (`exponents.out`). Here c₁ = l/w + 1 comes from the seed's dual (BGS App. F). |
| Zhang–Li–Li, *Network Coding Can Beat Routing in Undirected Multiple-Unicast Networks*, arXiv 2610.09367 (7 Oct 2026) | A finite connected subcubic network with unit capacities and leaf terminals, where a binary linear block code beats fractional routing (Thm III.1). The seed is a reversible integer circuit, and the routing bound is a *length certificate* ρ = m/(kb), the AFKL accounting in reverse. Thm VI.2 amplifies it by a degree-preserving tensor plus subdivision: for every β, girth ≥ n^β, coding → 1, routing ≤ K_β/(log n)^c. Works over every field and finite abelian group, and is Lean-formalized. | c = ln(1/ρ₀)/ln(4Γ²) (eq. 31). Γ ≥ b₀ > 7.8·10⁵, so c ≤ 0.035·ln(1/ρ₀). ρ₀ is not given numerically and is presumably near 1. |

**Also noted.** BH §1 records that OpenAI's *Finite tensor savings and exact Fourier circuits* (preprint, 25 Sep 2026, listed in `upstream/formalization.yaml:377`) "also implies a refutation": a sub-(n log n) circuit plus the AFKL cyclic-shift argument. That is exactly Theorem A's certificate form, used in the contrapositive.

**Assessment of the repo's record.** The repo's description (PR288 notes `:260-303`, `stronger-rank-screens.md:164-180`) is accurate. It should say that the gap exponent is unspecified and ≤ ~0.008 from the BH seed, not "(lg|G|)^c, c < 1" with c anywhere in that range.

## 4. What this means for κ

**Where today's κ sits.** A padded oblivious upstream multiplier at κ = 0.0017 is a coding-gap family with exponent ≥ 0.0017. That is consistent with, and comparable to, the BH and BGS ε ≲ 0.008.

**What higher κ would require.**

| target κ | implied gap exponent | how that compares with known families |
|---|---|---|
| 0.017 (×10) | ≥ 0.017 | larger than any exponent currently guaranteed by the literature |
| 0.317 (the paused Tanner value) | ≥ 0.317 | — |
| → 1 | (lg N)^{1−o(1)} | saturates flow–cut, an extremal network-coding result |

These are necessary conditions on the network-coding side. They are not barriers: nothing *forbids* η_* = 1.

**Non-oblivious escape.** The only escape from this cap is data-dependent head motion. Corollary B is vacuous for general TMs, and no multiplication mechanism of that kind is known.

## 5. Next step and stop test

1. **(Light, ~1 h) Settle the HvdH assumption.** Grep HvdH Annals 2021 §§2–5 for value-dependent lengths (normalization, leading-zero stripping, data-dependent branches) in the exponential/Gaussian routines and the multiplier.
   - *Stop:* if any step needs value-dependent length that cannot be padded within O(1), restrict Corollary B to "upstream with an oblivious substitute for that routine".
2. **(Optional) Compute c₁ for the BH seed.** Use the dual of the 147/157 LP and the published deleted-edge list (`checks/bh_sign_gauges.py` data). This gives the exact guaranteed ε to compare with κ.
   - *Stop:* if ε < 0.0017, an oblivious upstream already beats every published gap exponent. That would be worth a discussion post on the network-coding side, with no κ change.

## 6. Citations

✓ = read in full (local copies, not redistributed); R = as recorded in this repository.

- ✓ P. Afshani, C. B. Freksen, L. Kamma, K. G. Larsen, *Lower Bounds for Multiplication via Network Coding*, ICALP 2019, arXiv 1902.10935. §2 (Conj 1, p. 6), §4 (Thm 2 proof, p. 8), §5 (supervisor, eqs. 1–4, pp. 10–11).
- ✓ M. Braverman, S. Garg, A. Schvartzman, *Coding in Undirected Graphs Is Either Very Helpful or Not Helpful at All*, ITCS 2017, arXiv 1608.06545. Thm 1.1, Thm 3.10, App. F (c₂ = log(1 + ε/2)/log(16c₁²)).
- ✓ M. Braverman, Z. He, *The Multiple Unicast Conjecture is False*, arXiv 2610.10108 (2026-10-07).
- ✓ X. Zhang, B. Li, Z. Li, *Network Coding Can Beat Routing in Undirected Multiple-Unicast Networks*, arXiv 2610.09367 (2026-10-07). Lean: github.com/iqua/li-li-conjecture (not built).
- ✓ M. Braverman, Z. He, *Undirected Multicast Network Coding Gaps via Locally Decodable Codes*, arXiv 2510.18737 (abstract only).
- R OpenAI, *Finite tensor savings and exact Fourier circuits*, preprint, 25 Sep 2026 (cited by BH and in `upstream/formalization.yaml`; not read).
- ○ N. Linial, E. London, Y. Rabinovich (1995); Y. Aumann, Y. Rabani (1998): MCF within O(log k) of sparsity.
- ○ N. Pippenger, M. J. Fischer, JACM 26 (1979).
- ○ D. Harvey, J. van der Hoeven, *Integer multiplication in time O(n log n)*, Annals 2021.

## Scope

- Theorem A and Corollary B are elementary rewrites and have not been reviewed. The §1.3 constants are order-of-magnitude.
- The audit is a reading of §§02, 04–09 and 10 at the cited lines. §03 (motifs) is algebra with no tape motion beyond §02's primitives, and was skimmed only.
- No supplier, outer or κ number is changed.

Prepared with Claude Code (Claude Opus 5.5) assistance; not independently reviewed.
