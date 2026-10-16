# HvdH obliviousness audit, Corollary B made unconditional, and the BH/BGS exponent

No κ claim. It follows up `AFKL-TRANSFER.md` in this directory. Base: main `0fdd7d82`; guidance refreshed against `fbdbcb53`.

κ values used:
- κ_sel = 0.0016727, the selected result on main;
- κ_ours = 0.0017149, our PR #2679;
- κ_best = 0.0017037, the best claimed result, PR #2474.

Multipliers are given over 0.0017.

Sources (local copies, not redistributed):
- the HvdH text is HAL hal-02070778v2, which is the Annals 2021 text. It was converted with `pdftotext -layout` → `hvdh.txt`. **"l.N" below means line N of that file.**
- the #634 sources are #634's notes 1 and 1a, fetched read-only at head `189b5771`.

Scripts:
- `bh_c1.py` → `bh_c1.out`;
- `bgs_dense.py` → `bgs_dense.out`.

Both use only the standard library, take under 1 s, (`python3 bh_c1.py`, `python3 bgs_dense.py`).

## 0. Answers

1. **Every value-dependent step in HvdH §§2–5 pads to an oblivious step at a cost of ×O(1).** There are exactly three kinds of data influence:
   - (D2) variable-length "standard binary" integers;
   - (D3) the data branch in Lemma 2.2 (Scaling);
   - (D4) sign-conditional handling in the rounding ρ₀, the Kronecker packing of Lemma 2.5, and signed calls to M.

   Every other branch, loop count, window, sort comparison or jump is a function of indices, and so of n. The table is in §1.

   One routine needs a small explicit substitute rather than plain padding:
   - **Lemma 2.2's square root and division.** Clamp a to a′ = max(a, 2^{2p}+1), always run the "else" branch, then select. This keeps the dummy branch well defined (b′ ≥ 2^p+1, never 0) and its operands at a fixed width.
   - **The unspecified auxiliary O(p^{1+δ}) multiplier/divider of §2.1.** Instantiate it with the padded HvdH multiplier itself. Every call has size < n, so obliviousness follows by strong induction on n.
2. **Corollary B therefore holds unconditionally for upstream's padded multiplier.** "Unconditionally" here means with no remaining assumption on HvdH. It still rests on upstream's Theorem 1 and on our unreviewed Theorem A. The precise statement is in §2.
3. **The BH seed has c₁ = 67.94** (computed exactly; see §3). The guaranteed exponent is **ε = 0.00298** with BGS Thm 1.1 as printed, and 0.00412 with α optimised inside the same proof. **It is not below 0.0017**, so the claim "an oblivious upstream beats every published guarantee" **does not fire**:
   - κ_ours/ε = 0.58 along BGS's sequence.
   - For the all-N exponent (the liminf in Corollary B), BGS's own App. F recursion still guarantees 0.00198 (or 0.00260 with α optimised). Our κ is 0.87× (0.66×) of that.
   - Only reading the printed c₂ loosely as an all-N bound (c₂/2 = 0.00149) would put κ above it. That is an artefact of slack, not a result.
4. **Overlap is heavy: #634 (ShikeChen01, closed, retained in the 4th-100 review) already has most of this.**
   - Its Lemma B′ gives an oblivious machine at the same κ, its Cor. B″ gives the circuit form, and its §1.5 states "κ ≤ γ" as an open short claim.
   - Its note 1a audited HvdH Lemmas 2.1–2.5, 2.12–2.14 and 4.8–4.12 and found the same D2/D3/D4.
   - It left HvdH Thm 1.1 internals as **"Not established; replaceable"** (row 13), proposing Schönhage–Strassen plus the `08-assembly.tex:996` p^κ margin.

   **New here:**
   - (a) the audit of HvdH §3, Prop. 4.7's sort, Thm 4.1 and §5, which closes row 13 *without* the SS replacement or the margin;
   - (b) the clamp substitute for Lemma 2.2;
   - (c) the quantitative Theorem A transfer, which settles #634 §1.5;
   - (d) the first explicit exponent for the BH⊗BGS family, including an independent check of BH's flow certificate.

   **No discussion thread exists on obliviousness or network coding.** The searches covered discussions for "oblivious", "network coding" and "634", and PRs/issues for oblivious, coding gap, AFKL, Li–Li and Braverman. The hits were only #634, #288 and #349. The canonical place is therefore **a new single discussion, "κ and undirected coding gaps (oblivious padded multiplier)"**, with a short pointer on #1167 (larger-κ structure hub) and a one-line credit/pointer comment on #634. Do not post the "beats every guarantee" sentence.

**Multiplier: ×1.** This is a necessary condition on the network-coding side, not a construction.

## 1. Audit of HvdH 2021 §§2–5

**What upstream cites.** Upstream cites:
- Thm 1.1 as the "established multiplier" (`06-transforms.tex:454-460`, `07-resampling.tex:42`);
- Lemma 2.12 (`07:44-47`);
- Lemmas 4.8, 4.9 and 4.12 (`07:150-152`), which pull in Lemmas 2.13–2.14, Cor. 2.10, Lemma 2.4 and Lemma 4.11;
- Thm 4.2, Lemma 4.6 and Lemma 5.1, which are mathematics only (`07:107`, `08:123`).

So the whole of §§2–5 is in scope through Thm 1.1.

**Model.** HvdH use a multitape TM with a fixed number of tapes and "integers encoded in the usual binary representation" (l.36-37, l.49-50, §2.1 l.422).

**Obliviousness criterion.** A machine is oblivious when its head positions are a function of n. A step passes if:
- **(N)** its inputs are functions of n;
- **(A)** its extents are fixed by indices or counters; or
- **(F)** it is fixed-width arithmetic with carries or signs held in finite control.

A key observation is that **every exponential and Gaussian weight in HvdH is computed from indices (j, k, s, t, α, σ), never from the data.** Whatever their internal branching, those routines are (N).

| # | HvdH routine (section, line) | Value-dependent length, loop or head motion? | Fixed-width padding oblivious at O(1)? |
|---|---|---|---|
| 1 | §2.1 integer format, "standard binary" (l.422); §2.2 C̃∘ element as an integer pair (x, y), \|x\|,\|y\| ≤ 2^p (l.436) | **Yes (D2):** lengths depend on values | Yes. Use sign plus magnitude padded to p+1 bits, and n+1 bits for M(n) inputs. HvdH already charge O(p) per coordinate. |
| 2 | §2.1 auxiliary M(n) = O(n^{1+δ}), quotients ⌊x/y⌋ and ⌊(x/y)^{a/b}⌋ (l.427-431); no algorithm specified | Depends on the instantiation | Yes, given a choice. Data-touching uses are only Lemma 2.4/Cor 2.10 products and Lemma 2.2. Instantiate with padded HvdH recursively (calls have size p ≪ n, so induction) or with padded SS. Divisions and roots occur on data only in Lemma 2.2; see row 4. |
| 3 | ρ₀ round-toward-zero, floor or ceil by sign (l.472); Lemma 2.1 ½(u±v) (l.511) | **Yes (D4):** the sign picks floor or ceil | Yes. In sign-magnitude, round-to-zero is a plain truncation of the magnitude, with no branch. |
| 4 | **Lemma 2.2 Scaling (l.527-560)**: "If a ≤ 2^{2p} … take ṽ := cũ" (l.534), "Suppose instead a > 2^{2p} … b := ⌈a^{1/2}⌉ … ρ(2^p cũ/b)" (l.540) | **Yes (D3):** a data branch, and a √/÷ of data | Yes, with a substitute. Compute a′ = max(a, 2^{2p}+1) by a fixed-width compare and mux, then always run the else-branch on a′. Then 2^p < b′ < 2^{2p+1}, since \|x\|,\|y\| ≤ c·2^p ≤ 2^{2p}, so the operands fit fixed widths of 2p+2 and 4p+2 bits. Newton √/÷ uses a fixed iteration count. The leading-bit normalisation varies over ≤ p+1 positions: an oblivious barrel shift costs O(p log p) ⊂ O(p^{1+δ}). Finally, mux between cũ and ṽ. Total ×O(1) on O(mp^{1+δ}). Used in Prop 3.4 (l.1086), Prop 5.2 (l.2192) and Prop 5.3 (l.2230). |
| 5 | Lemma 2.4 / Cor 2.10, product in C̃∘ (l.711-724, l.819) | Only via D2/D4 | Yes. Fixed-width p-bit Gaussian-integer product, then truncation. |
| 6 | Lemma 2.5, Kronecker: "zero-padding (or one-padding … negative) plus carry and sign handling" (l.738); 4·M(3rp) on signed integers; unpack "read off last 3p bits … stripping" (l.744); scale and round (l.751) | **Yes (D4):** sign-dependent passes | Yes. Slots are fixed at 3p bits, and carries and borrows are single bits in finite control. For the signed M call: always complement-and-increment and select (as in upstream 06:522), multiply the padded magnitudes, then always negate and select. O(rp) extra, already charged. |
| 7 | §2.3 slice rearrangement / fast transposition [4, App.] (l.584); Lemma 2.11 tensor products (l.828) | No: recursion on fixed shapes | Yes, once records are fixed width (A). |
| 8 | §2.8 π, log 2, e^z (l.870-879); Lemmas 2.12–2.14 (l.881-933), incl. the branch "check whether j > kp" (l.902) and τ (l.924) | Branches only on indices | (N): already oblivious. |
| 9 | §3 Lemma 3.2 FFT over R: butterflies (l.996), multiplication by y^{−2rk/t} as a cyclic shift with sign flips, interleave | Shift amount is index-driven | Yes (A+F). Negation of a fixed-width word is a fixed pass. |
| 10 | Props 3.3–3.5 (l.1019-1155): multivariate FFT, pointwise Lemma 2.5, Lemma 2.2 scaling, S̃/T̃/Ũ with ζ_j from Lemma 2.12 | Only through rows 4 and 6 | Yes. |
| 11 | Thm 3.1 Bluestein: η_{j₁…j_d} by incremental (j+1)² over index tuples (l.1175-1182); pre- and post-multiply | No | Yes (N, then row 5). |
| 12 | §4 Lemma 4.8 D̃′ (l.1641): weights via Lemma 2.14, times data | No | Yes (N, then row 5). |
| 13 | Lemma 4.9 S̃′ and Remark 4.10 (l.1662-1756): window \|j−sk/t\| < m, the head "jumping over gaps" at exceptional rows, the sum Σa_j (l.1721); \|T̃_k\| ≤ 1 is proved, not tested | Jumps depend on k, s, t, m only | Yes (A). The accumulator has fixed width p+O(log m). |
| 14 | Lemma 4.11 Ẽ (l.1757): exactly 2m terms | No | Yes. |
| 15 | Lemma 4.12 J̃′ (l.1820-1863): Neumann iteration count n = ⌈p/α²θ⌉ (l.1826), no convergence test | No | Yes. |
| 16 | Prop 4.7(ii) (l.1864-1890): C̃ "copies its input in order, skipping t−s" entries (l.1879); P̃_s⁻¹ and P̃_t by "label-and-sort … merge sort" (l.1881-1883) | Merge comparisons are on labels tj mod s, which are index functions | Yes. The comparison outcomes, and hence head motion, are n-determined. Records must be fixed width. (Upstream §07 removes these sorts anyway.) |
| 17 | Thm 4.1 proof (l.1893-1910): tensor via Lemma 2.11 | No | Yes. |
| 18 | §5.1 parameters b, p, α, γ, T, r, t_i (l.2079); primes by Eratosthenes (l.2135) | No | (N). |
| 19 | Prop 5.4 (1): split f, g into N = ⌈n/b⌉ b-bit chunks (l.2243); overlap-add of 3b-bit H_j (l.2258) | Only if inputs are stripped | Yes. Require exactly n input bits, padded. Do the overlap-add as a fixed number of fixed-width passes with the carry in finite control (F). |
| 20 | Prop 5.4 (2): CRT layout "by attaching labels and sorting" (l.2273) | No (index labels) | Yes, as row 16. |
| 21 | Prop 5.4 (3): multiply by 2^{2b}S and round to the nearest integer (l.2297) | Rounding carry | Yes. Run the rounding increment over the full fixed width, not stopping early. |
| 22 | Base case n < n₀, "any convenient … classical O(n²)" (l.2041) | Schoolbook "add if bit = 1" | Yes. Use AND-gated additions. n₀ is a constant, so cost is O(1) either way. |
| 23 | Recursion M(n) → (12T/r)·M(3rp) (Prop 5.4, Cor 5.5) | Sizes are functions of n | Yes. The whole machine is oblivious by strong induction on n. |

We looked for these and found none in §§2–5: data-keyed sorting, convergence tests or early exits, run-time overflow or clipping tests (every disk bound is proved), and data-dependent recursion. §4.4 and §5.4 are remarks and are not used.

**Conclusion of (1).** Padded HvdH (rows 1–6 fixed) is an oblivious multitape multiplier with M_obl(N) = O(N log N). No Schönhage–Strassen replacement is needed, so #634's reliance on the `08-assembly.tex:996` margin "Any fixed power of log p is at most p^κ" is unnecessary. (That route is still valid as an alternative.)

## 2. Corollary B, unconditional form

**Definitions.**
- Let γ(N) be the supremum of NC/MCF over undirected unit-capacity multiple-unicast instances with ≤ N edges (`AFKL-TRANSFER.md` §1.1).
- Let M_obl be upstream's multiplier with `AFKL-TRANSFER.md` §2's two fixes:
  - sign passes always made (06:522-540, 08:757-763);
  - padded HvdH adapters (07:202-240);
  - and with every cited HvdH routine replaced by its padded form from §1 (rows 1–6).

**Corollary B (padded upstream).**
- M_obl is an oblivious fixed-tape TM.
- If upstream's Theorem 1 holds with exponent κ, then M_obl multiplies n-bit integers in time O(n·lg^{1−κ} n).
- By Theorem A of `AFKL-TRANSFER.md` applied through the standard oblivious-TM-to-circuit simulation,
    γ(N) = Ω((lg N)^κ) for all large N, i.e. η_* := liminf log γ(N)/log log N ≥ κ.
- For main's selected supplier, η_* ≥ κ_sel = 0.0016727.

The ingredients still used are only upstream Theorem 1 itself and Theorem A of `AFKL-TRANSFER.md` / Cor B, all of which are unreviewed.

**Alternative proof, cited from #634 §1.2.** For this corollary, a *circuit* for HvdH (they note their results hold in the Boolean circuit model) spliced in at n-determined positions also suffices. That gives a second proof of the corollary that does not use §1.

## 3. c₁ for the Braverman–He seed and the guaranteed exponent

**Seed.** BH §4 and App. A: H = PG(2,9) incidence (182 vertices, 910 edges); R = 175 deleted edges; sessions B = R∖(C∪J), |B| = 157; G has 735 edges.

`bh_c1.py` rebuilds G from App. A (F₉ = F₃[θ]/(θ²+1), the ordered triples, and the R list). It checks the following:
- |B| = 157 and |E(G)| = 735.
- **All 157 session distances are exactly 5.** This is an independent check of BH's flow certificate 157·5 = 785 > 735.

**Dual.** The dual is uniform, w_e = 1/785. It is feasible because Σ_i dist = 157·5·w = 1, and its objective is z(D) = 735/785 = 147/157. So 1+ε = 157/147, with every weight positive.

**α-extension (BGS App. C).** The terminals share vertices (a vertex meets up to 5 deleted edges), so the α-extension is needed:
- add K = 2·157 = 314 leaves, each of capacity z(D)(1+ε) = 1;
- give the leaves weight α/(K(1+ε)). The script asserts that the leaves add exactly α·z(D).

**c₁.** By BGS App. F, c = l/w, where l is the largest session distance in D* and w is the smallest weight in D*.
- The smallest weight is a leaf weight: 0.0770 of an edge weight at BGS's α = 5/152.
- So c = (5/785 + 2w_leaf)/w_leaf = **66.94**, and **c₁ = c + 1 = 67.94**.

| exponent (BH ⊗ BGS) | BGS α = 5/152 | α optimised (same proof) | κ_ours / it |
|---|---|---|---|
| Thm 1.1 c₂ = ln(1+ε/2)/ln(16c₁²), along the sequence | **0.00298** | 0.00412 (c₁ = 262) | 0.58 / 0.42 |
| App. F recursion, exact inequality, along the sequence | 0.00396 | 0.00519 | 0.43 / 0.33 |
| **all-N (liminf; the quantity in Cor. B)** from the App. F recursion, = sequence/2 | **0.00198** | 0.00260 | 0.87 / 0.66 |
| printed c₂ read as all-N, c₂/2 (loose) | 0.00149 | 0.00206 | 1.15 / 0.83 |

**Why the all-N exponent is half.** BGS give gaps only on a doubly exponential sequence. Each step squares the gap but doubles log log(size). For N just below m_{i+1}, the best guaranteed gap is still gap_i.

**Sensitivity.** If App. C's "k" counted pairs rather than leaves (157 rather than 314), the leaf weight doubles, giving c₁ ≈ 35.5 and c₂ ≈ 0.00337: larger again.

**Other families.**
- A non-uniform dual cannot change the picture much: the leaf weight dominates w.
- ZLL's c = ln(1/ρ₀)/ln(4Γ²) ≤ 0.035·ln(1/ρ₀), with ρ₀ not stated numerically. It would beat 0.002 only if ρ₀ < 0.94; not computed.
- The 33-vertex seed of dysfunctor/#634 (149/150) is far weaker.

**Verdict on (3).** ε_guaranteed ≥ 0.00298 > 0.0017, so the claim **does not fire**.

**Precise statement that does hold.** If upstream's Theorem 1 holds, the padded oblivious multiplier gives an explicit, *every-N* family with η ≥ κ_sel = 0.0016727. For comparison:
- That is within a factor 0.85 of the best all-N exponent guaranteed by BGS's recursion for the BH seed (0.00198, or 0.00260 with α optimised).
- It beats that guarantee only once κ > 0.00198 (×1.16 over 0.0017), or κ > 0.00260 (×1.53) against the optimised-α version.
- It beats the along-sequence guarantee only once κ > 0.00412 (×2.4).

## 4. Overlap and where to post

- **#634** (ShikeChen01, closed; retained in `docs/research/community-reconciliation-fourth100-review.md:47,127` as "oblivious compilation … not established").
  - Its note 1 §1.2 Lemma B′ and Cor B″, and its note 1a rows 1–27, are the prior art for §1 and for Corollary B's premise.
  - Its §1.5 lists "κ ≤ γ [P, short]" as open: that is Corollary B, which Theorem A of `AFKL-TRANSFER.md` supplies with constants.
  - It says the "implied gap (log N)^κ is consistent with the published (log n)^c families" and "no useful value of γ is known". §3 gives the first numbers.
- **No discussion covers obliviousness or network coding.** The other hits are only #288 and #349 (PR288 sign-gauge, which cite BH qualitatively).
- **Canonical thread.** Open one new discussion: "κ and undirected coding gaps: the padded upstream multiplier is oblivious (HvdH §§2–5 audit), η_* ≥ κ, BH⊗BGS exponent 0.0020–0.0041". Give it the standard 100–150-word status block. Then add:
  - a pointer on #1167 (it bounds what a framework change must achieve on the NC side);
  - a one-line comment on #634 crediting note 1a and saying that row 13 is now closed without SS.

  #1166 is not needed, since no independent check is being requested. A PR is not warranted: this is a scoped research note with a ×1 multiplier, unless the maintainers want the audit table in `docs/research/`.

## 5. Next step and stop test

1. **(Optional, ~1 h) Make "η_* ≥ κ_sel" rigorous on the BGS side.** Check whether tensoring A_i with A_{i−1} or with the seed fills the doubly exponential gaps. That would raise the all-N guarantee toward the sequence value, 0.0040–0.0052.
   - *Stop:* if mixed tensoring is covered by BGS Lemma 3.8 without extra girth cost, the all-N guarantee is ≈ the sequence value. Upstream's κ is then 0.33–0.43× of it, and the "comparable to known families" remark becomes "below known families by ×2–3".
2. **(Optional) Compute ZLL's ρ₀ from their Lean seed parameters**, to close "every published guarantee".
   - *Stop:* if ρ₀ ≥ 0.95, ZLL's exponent is < 0.002 and BH⊗BGS is the best published guarantee.

## Scope

- §1 is a reading audit of HAL v2 at the cited lines. It is not a formal TM construction.
- Theorem A and Corollary B (`AFKL-TRANSFER.md`) are unreviewed.
- §3's numbers are bounds read off the BGS proofs, not new constructions.
- No supplier, outer stage or κ number is changed. Nothing was posted.

Prepared with Claude Code (Claude Opus 5.5) assistance; not independently reviewed.
