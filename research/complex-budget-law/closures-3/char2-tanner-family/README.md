# Family record (rev-5 expansion gate): char-2 Hermitian Tanner codes on PG(2,q) − M. Decision: PAUSE

Status: **family record, not a closure.** LOCAL-ONLY, author-checked, not reviewed, no PR, no canonical thread yet. No κ claim.
Base: main `0fdd7d82`; guidance revision 5 (`docs/research/research-family-status.md` at `fbdbcb53`). Multipliers are over
κ = 0.0017. This record follows the five points that revision 5 asks for before a family is expanded.

## Read this first: the family is BIT-ONLY

κ = min(complex, bit), and the complex supplier binds today (it consumes bit calls through its change-of-basis adapters, so
κ ≤ min(a_bit, a_complex) holds structurally). Every result in this family is a **bit-side** primitive. **On its own it is worth
×1.0 in κ**, however large its bit-side gain. With the complex side at its realistic framework cap (a ≈ 0.03), a bit-only route
is worth ×17 at most (expL6 §E4, `l6/estimates.py`). The literal complex twin of the layout is dead: the char-0 Teichmüller lift
has S ≤ 141 (expNK), which gives a ≤ −1.11 by expQK's formula.

## Decision: **PAUSE** (upgraded from expL6's draft "NARROW")

Reasons:
1. **Bit-only** (above): ×1.0 in κ until some complex-side mechanism exists, and none is known.
2. **Every tested frame family is certified ×0.** GRS gives nullity 0 on all 142 fixed-scale runs; the ray test (all
   scalings, any denominator) has a trivial gcd for every natural char-2 label family; F1 (cut-aligned frames) has rank drop 0
   on every cut (expQR). The smallest certified rigid + obstructed gadget (q = 8) routes integrally, so it is ×0 for every
   frame (expQI).
3. **No exact schedule.** The best |R| is 784 against S = 772 at q = 16 (57657 against 57360 at q = 64). Every algebraic block
   layout is killed by a divisibility theorem (expQU); the cut + viol bound is the binding quantity (expQD).
4. **No code exceeds 772** on `schedule784`: dim W of every Frobenius-twisted GRS pair equals a Hilbert function, maximal at
   772, and any single-column replacement stays ≤ 781 < 784 (expQE).

**Reopen triggers** (any one):
- an **exact char-2 schedule** (|R| = S, a realised causal circuit) **and** Lemma 3's first-order 2-adic GRS test passing on it
  (expQR §2.4; the test can run only once a realised char-2 circuit exists);
- a **frame family in B** = {N : det Z_A(N) = 0 for all realisations A} (expQR Lemma 3), i.e. one that survives `qr/ray.py`;
- a **complex-side counterpart**: a char-0-only coding seed with S/E > 1/6 and no char-2 realisation (expL6 C1). **C1 is now
  closed** (`../complex-coding-seeds/`): flat-leverage (ζ and Walsh–Hadamard) codes have S0/E ≤ 1/4 on any layout and a ≤ 0 at
  every girth ≡ 0 mod 4, every char-0 gap reappears mod 2, and only an unreached sliver at girth 6/10 remains. **The complex side
  has no coding-seed route**, so this trigger has no known candidate mechanism.

Without an outer contract that lifts the ¼ cap (expL6 O1/O2), even full success is capped at κ < ¼.

**Inputs.**
- **expQR (GRS): filled in below (§5).**
- **expRA (MC_4): filled in (§5b). Not certified;** the decision stays PAUSE.

## 1. Sources, target and evidence status

**Declared target.** A **bit primitive** (local module, not a complete supplier) in the `finite_bit_contract` / 04-swap
contract: a helper-free, garbage-free causal circuit of square [16,8] MDS gates over F16 (F2-linear 32-bit pointwise gates) on
PG(2,16) − M, realising its sessions with rational frames, with s < W·m. Milestone: local module; no integrated κ is targeted.

| source | content | status | replayed here |
|---|---|---|---|
| expNC | the seed: Hermitian x⁴ code, dim W = 772 on PG(2,16) − M, sound split, δ_nd ≥ 0.0684; odd-degree 900 / 0.1687 (non-circuit) | author-checked, exact | via `qk/` (772) |
| expNK | char-0 screen: q = 8 rigid and obstructed; generic char-0 dim 0 at q = 16; Teichmüller lift ≤ 141 | author-checked | no |
| expQB | Singer-block model: q = 16 rigid (39311 = gauge) and first-order obstructed; q = 4 not rigid | author-checked, exact | `qi/qb.py` used as a library |
| expNI | cyclic schedule \|R\| = 784 (`data/schedule784.json`, verified); cut lemma \|R\| ≥ bisection/2; spectral bound ≥ 751 | author-checked | `qe/verify_sched.py` |
| expNJ | small MDS gadgets: no realisable deficit (unrealisable in every field, or full flow F_* = W ⇒ s ≥ Wm for every frame). **No local gadget exists; the frame question is global** | author-checked, exact census | no |
| expNY | the generic-realisation screen (GRS): s ≥ rank_Q Z_G for any char-0 gate assignment; only terminal frames enter | author-checked | `ny/check_telescope.py` |
| expQK | S(q) = q(q² − q^{3/2} + 1)/4 (exact at q = 4, 16, 64); a → 0.317; q = 64 cyclic layout 57657 vs S = 57360 | author-checked | `qk/qk_singer.py` (q = 4, 16, 64) |
| expQD | Singer-graph bisection with **cut 1542** exists (so the bisection obstruction is false); joint (S, M) search to 1534; **cut + viol lemma** | author-checked | `qd/verify_cut.py` |
| expQE | max dim W = 772 = h(7,7) (**Hilbert function match**); bit-exact `schedule784` simulation rank(Φ − I) = 12; single-column replacement ≤ 781 | author-checked, exact | `qe/hilb.py`, `qe/verify_sched.py`, `qe/colrep_bound.py` |
| expQU | **divisibility theorem**: block layouts have \|R\| divisible by \|H\|; gcd(S, n) = gcd(√q − 2, 3); Baer nesting obstruction | author-checked | `qu/divcheck.py` |
| expQI | q = 8 (smallest rigid + obstructed gadget) **routes integrally**, F_* = W, so ×0 for every frame; q = 16 SA finds nothing for m ≤ 4 | author-checked | `qi/flow.py` |
| expQV | GQ T2*(O) seed: every code misses the cut requirement; charged a ≤ 0.086 (optimistic) | author-checked | `qv/table.py` |
| expQR | GRS on `schedule784` + 12 C/J lanes: 142 runs nullity 0; ray test trivial gcd; F1 dead; Lemmas 1–3 | author-checked | `qr/` |
| expL6 | draft five-point record, five changed contracts, planning arithmetic | synthesis | `l6/estimates.py` |
| expRB, expRF | C1 (complex-side coding seeds) closed: no coding-seed route on the complex side | author-checked | `../complex-coding-seeds/` |
| expRA | the frame-independent quantised-metric screen MC_m ≥ Wm (m ≤ 4) on q16 and `schedule784`: **not certified**; exact LP = m·E_int/5; MC_5 exact with room | author-checked; rational certificates | `ra/check_flow.py` |

## 2. Boundary specification

| item | specification |
|---|---|
| inputs | the S session symbols (F16, i.e. 4S bits per address) on the source ends of the R back arcs of a causal Eulerian 8-in/8-out orientation |
| outputs | the same symbols at the sink ends (sound ⇔ the 273-step monodromy Φ = I ⇔ dim W = \|R\|; `qe/verify_sched.py`) |
| role permutation ρ | the session own-pair map. **Unknown** whether it composes into the bank exchange the 04-swap endpoint needs, or how session roles are packed into the X/Y banks. Fixed roles would climb the full m (Lemma W′) |
| frames | rational m × m on every gate and terminal. Must have 2-adic denominators (T1), be non-diagonal (T2), have no odd-characteristic realisation (C3), and lie in B (GRS, expQR Lemma 3) |
| owned / borrowed state | none if \|R\| = S. If \|R\| > S: \|R\| − S zero-ancilla (C) and garbage (J) lanes |
| required reads | each gate reads its 8 in-arcs (32 bits) pointwise |
| dirty compensation | none when exact. Otherwise the 12-dimensional defect (784 vs 772) needs C/J lanes; expN1 scenario B prices their completion at δ ≈ 0 |
| restoration | helper-free ⇒ nothing. With C/J: uncompute, priced as above |
| field | F16 gates ⇒ char-2-only (rigid and obstructed). By construction **no complex (char-0) twin** of the same layout |

## 3. Obligation ledger (proved / constructively upper-bounded / unknown; unknowns are not zero)

| obligation | status | value / evidence |
|---|---|---|
| coding gap dim W = 772 (q = 16) | **proved** (exact, three independent models) | NC, NI, QB; `qk/` reproduces 9 / 772 / 57360 |
| why 772: global-form count | **proved** for Frobenius-twisted GRS pairs: dim W = h(7,7) of the twisted incidence ring, 512 / 708 / 772 / 708 | QE, `qe/hilb.py` |
| soundness (information-set split) | **proved** at q = 16 | NC, QB |
| char-2-only at this code (q = 8, 16) | **proved** (rigid and first-order obstructed) | NK, QB |
| no char-0 realisation on *any* component | **unknown** | NK §4, QB §3 |
| routing ceiling δ ≤ 0.0684 | **proved upper bound** (distance-5 routing) | NC |
| exact causal schedule \|R\| = S | **unknown; every structured attempt fails.** Constructive 784; spectral lower bound 751. Bisection is not the obstruction (cut 1542 exists, `qd/`), but **cut + viol** is: the 1542 bisection has internal degree ≥ 9 on every vertex, so cut + viol ≥ 1546 > 1544; the best known min_S[cut + viol] is 1568 (→ 784). Block layouts: \|R\| ≡ 0 mod \|H\| and gcd(S, n) = 1 at q = 16, so only \|H\| = 1 remains (QU) | NI, QD, QU |
| a code making 784 exact | **upper-bounded at 772** for global-form GRS codes; any single-column replacement ≤ 781 (`qe/colrep_bound.py`); general **unknown** | QE |
| C/J completion of the 784 schedule | **upper-bounded at ×0** (completion charge) | N1 scenario B |
| frames with s < Wm | **unknown, and every tested family is certified dead.** q = 8: integral routing ⇒ ×0 for every frame (QI). q = 16: GRS nullity 0 on 142 runs, trivial ray gcd (all scalings), F1 drop 0 (QR). The quantised metric needs m ≥ 5 (QI, evidence) | QI, QR, NJ |
| frame-independent GRS certificate | **neither proved nor refuted.** Positive and switch certificates are provably impossible (all 297,606 session pairs mutually reachable; switch cycles ≤ 728 < 772). A single realisation kills all frames iff L_A is a space of nilpotent matrices (Lemma 2) | QR |
| frame efficiency η (share of the flow cap realised) | **unknown**. The ceilings assume η ≈ 0.3–0.73 | — |
| complex counterpart (κ = min) | **upper-bounded dead for coding seeds**: the literal twin (S ≤ 141 ⇒ a ≤ −1.11), and every flat-leverage char-0 code (C1 closed: S0/E ≤ 1/4, a ≤ 0 at girth ≡ 0 mod 4, girth-6/10 sliver unreached). Without a complex mechanism the κ value is **≤ ×1.0** | NK; expL6 §E3/§E4; expRB, expRF (`../complex-coding-seeds/`) |
| outer conversion | **proved constraints of the current recipe**: κ ≤ min(a/(1+a), ¼); checker guard κ < 1/33. q = 64's κ ≈ 0.164 already exceeds the guard | assembly-domain-audit |
| larger q: S(q) closed form | **proved at q = 4, 16, 64** (computed); general form conjectural. Schedule slack to the spectral bound shrinks like ~2q^{−3/2}; cyclic layouts miss by more | QK |
| other geometries (GQ T2*(O)) | **upper-bounded**: S/E stays below the cut requirement at q = 8 and 16 (deficit 4.6% → 6.5%); spectrally charged a ≤ 0.086 at q = 16 (optimistic) | QV, `qv/table.py` |
| hosting, precision, bootstrap stock, all-size uniformity | **unknown, not priced** | — |
| amplify-then-compile, odd-degree full PG | **closed** (Lemma W′; completion charge) | QN, N1 |

## 4. Why a next subtask could change the decision

- **What would disappear:** helpers and completion lanes (with an exact \|R\| = S schedule), and the complex template (with a
  complex-side counterpart).
- **What replaces them:** coding gates, and 2-adic non-diagonal frames on the whole network.
- **Why it would not come back:** an exact sound schedule has no C/J lanes to uncompute, and frames are internal.
- **Allowance for the missing modules.** Frames: to beat the bit record a_bit = 0.0017334 at L = ln 8, δ ≥ 0.0036 suffices
  (η ≥ 0.05 of the 0.0684 cap). For κ the allowance is irrelevant until the complex side moves. **That fact dominates.**
- **What can still flip it:** only the reopen triggers above. The decisive tests that were running when expL6 drafted
  "NARROW" have now reported:
  - **expQR (GRS): fired against every frame family**, including F1 and all rescalings (§5). Not a frame-independent kill.
  - **expQD/expQU (schedule):** the bisection threshold did not fire (1542 exists), but cut + viol and the divisibility theorem
    leave no structured candidate.
  - **expQE (codes):** 772 is the maximum in the global-form family.
  - **expRA (MC_4): not certified** (§5b). m ≤ 4 frames are neither excluded nor found; it does not produce a reopen trigger.

## 5. expQR in detail (GRS on `schedule784`)

**Circuit.** `schedule784` (PG(2,16) − Singer matching d0 = 59, 546 gates, 4368 arcs, 784 back arcs) with 12 C + 12 J lanes,
W = 772, in two lane modes (`same`, `pair`). T_G for random integer 8×8 gates mod p = 33,554,393 is fully dense, rank 772, and
every source reaches every sink.

**Algebra.** With source frames N_k and sink frames I + N_k (the endpoint identity), Z_A(N) = (I − K_A(N))Ã with
K_A(N) = ÃNÃ⁻¹ − N (`qr/check_lemma.py`).
- **Lemma 1 (frames scale freely).** tN satisfies the contract for every t, so Z_A(tN) is singular ⇔ 1/t ∈ spec K_A(N); algebraic
  t of degree d is realised by rational frames of size md.
- **Ray test.** For a family direction N0, a trivial gcd of charpoly(K_A(N0)) over 3 independent realisations proves that **no
  rescaling at any denominator** survives GRS.
- **Lemma 2.** One realisation A kills every frame assignment iff L_A = {K_A(N)} is a space of nilpotent matrices. Necessary:
  A_ii(A⁻¹)_ii = 1 and B + Bᵀ = 2I (B = A∘A⁻ᵀ), 298,378 equations on a variety of dimension ≤ 31,360. **Positive and switch
  certificates are impossible** here (no cancellation-free triangular order exists since F_* < W; switch cycles ≤ 4368/6 = 728 <
  772). The Λ + rank-X route is vacuous (needs rank X ≲ 26; the dimension count gives ≳ 593).
- **Lemma 3.** The route survives GRS iff some frame set lies in B = {N : det Z_A(N) = 0 ∀A}. Every genuine saving lies in B
  (telescoping), so GRS can never refute a real saving; a frame-independent kill needs a Nullstellensatz certificate. B is decided
  2-adically, and for denominator-2 frames the first-order test is det(I + [N̄, Ē]) ≡ 0 (mod 2) for every Ē in an affine F16
  space. **It can run only once a realised char-2 circuit exists.**

**Results.**

| test | result |
|---|---|
| nullity of Z_G, ten label families (D-roots, gate-orthogonal, rank m/2, E8 roots and rank-4 frames, F16/Hamming, realified Hermitian, two-class, common control), m = 2..8, 3 G seeds + `pair` | **142 runs, max nullity 0** (`qr/fam_summary.txt`) |
| ray test, 3 independent G, m = 2, 4, 8 | **gcd degree 0 for every family** |
| F1 cut-aligned frames: rank drop of T[ρA,A] + T[ρĀ,Ā] on singsa 1542, the order prefix, kah1544 and pg1534 layouts, 3 G, 3 side rules | **0 everywhere** (threshold 0.02·W ≈ 15–21). F1 is dead by its own stop rule |

**Reading.** GRS is decisive against **frame families**, not against the route. With generic gates, rank T[ρA,A] equals the
directed arc-disjoint max-flow, and on these 16-regular expanders it is never below |A|; nothing forces low-rank blocks.

## 5b. expRA: the quantised rank-metric screen MC_m (frame-independent, m ≤ 4)

MC_m := min over integer edge ranks r ≥ 0 of Σ_e r_e + Σ_k (m − d_r(c_k, p_k))⁺, on the internal edges. Any frames with s < Wm give
MC_m < Wm, so MC_m ≥ Wm would kill every frame assignment of size m, independently of the frame family. Instances: the split
layout `q16` (772 sessions, 3596 internal edges) and the causal schedule `s784` (784 sessions, 3584 internal edges).

| instance | W | exact LP / m (= E_int/5) | m = 4 | m = 5 |
|---|---|---|---|---|
| q16 | 772 | 719.2 | **2877 ≤ MC_4 ≤ 3088** | **MC_5 = 3596 exactly** (< Wm = 3860) |
| s784 | 784 | 716.8 | **2868 ≤ MC_4 ≤ 3136** | **MC_5 = 3584 exactly** (< Wm = 3920) |

- **Not certified: MC_m ≥ Wm for m ≤ 4 is still open.** Lower bounds are ⌈m·F⌉ from rational multiflows on 5-paths (`ra/`);
  upper bounds are r ≡ 0 (= Wm). Every heuristic (greedy, SA, isolation families, CP-SAT incumbents) gives MC_m = Wm, but that is
  evidence, not a proof.
- **Exact LP = m·E_int/5** on both instances: r ≡ m/5 is feasible, the dual matches, and pricing over all paths finds no violated
  row. The uniform 1/5 metric is LP-optimal and saturates every internal edge.
- **MC_5 is exact:** ⌈5F⌉ equals the r ≡ 1 cost E_int (every session path has length ≥ 5). So at m = 5 the metric relaxation leaves
  264 (q16) / 336 (s784) units of room: **m = 5 frames are not excluded by this screen**; there the open question is realisability.
- **Subadditivity** (MC_{a+b} ≤ MC_a + MC_b, with MC_m ≤ Wm) means MC_4 = 4W implies MC_1..3 = mW: **m = 4 is the whole claim.**
- **m = 1 reduces to a signed max-cut:** MC_1 = W − max_A [δ_R(A) − δ_int(A)]. Any MC_m ≥ Wm proof contains a max-cut optimality
  certificate on 546 vertices; the plain SDP provably cannot give it, and the cycle LP and CP-SAT bounds are far from 0.
- **Effect on the decision: none.** The screen neither excludes nor exhibits m ≤ 4 frames, and it is not a reopen trigger.
  Decision stays **PAUSE**. A future m = 1 max-cut certificate followed by an m = 4 branch-and-cut would turn it into a frame-
  independent kill for m ≤ 4.

## 6. Results replayed in this package (selected)

- **Cut and cut + viol (`qd/`).** The Singer graph has a balanced bisection with cut 1542 < 1544 = 2S, so expNI's cut obstruction
  does not rule out 772. But for any schedule and gap, |R| ≥ ½(c + viol), viol = Σ_v[(e_part(v, earlier) − 8)⁺ +
  (e_part(v, later) − 8)⁺]; the 1542 bisection has internal degree ≥ 9 everywhere, so it cannot be the middle of an exact schedule.
- **Hilbert function (`qe/hilb.py`).** h(7,7) of F2[P, L]/(twisted incidence) is 512 / 708 / 772 / 708 for the four twists, equal to
  the measured dim W. x⁴ is the complete intersection (4,1) + (1,4).
- **Schedule simulation (`qe/verify_sched.py`).** Bit-exact over F4096: all 17 gate types bijective, rank(Φ − I) = 12, dim W = 772,
  matching the Fourier count. So `schedule784` is not exact for x⁴.
- **Column replacement (`qe/colrep_bound.py`).** dim W(v) = base_d + #{η : v ∈ V_η} ≤ base_d + 273 ≤ 781 for every column.
- **Divisibility (`qu/divcheck.py`).** Block layouts (contiguous Singer-subgroup orbits) have |R| ≡ 0 mod |H| (exact min-cost flow,
  16 random orders). With gcd(S, n) = gcd(√q − 2, 3) = 1 at q = 16, only |H| = 1 can be exact.
- **dim W formula (`qk/qk_singer.py`).** S(q) = q(q² − q^{3/2} + 1)/4: 9, 772, 57360 at q = 4, 16, 64.
- **q = 8 routes integrally (`qi/flow.py`).** 82 edge-disjoint internal paths for every tested information set, so F_* = W and
  s ≥ Wm for every rational frame (expNJ C1).
- **GQ T2*(O) (`qv/table.py`).** Measured S/E 0.117 (q = 4), 0.131–0.133 (q = 8), 0.154–0.155 (q = 16), all below the spectral cut
  requirement, so no exact helper-free schedule exists for any of these codes. Spectrally charged (optimistic) a ≤ 0.086 at
  q = 16; at q = 8 the RS code charges to −0.003 and the non-MDS duplicate-column code (S = 680) to +0.015. Soundness was checked
  only for the RS code.
- **GRS telescoping (`ny/check_telescope.py`).** The identity Σ_e (b_e a_eᵀ) ⊗ (M_head − M_tail) = Z_G holds exactly mod p on
  Braverman–He scenario A with random non-commuting frames. This is the tool expQR applies.

## Files and reproduction

Run each from its subdirectory. Runtimes under a 4G cap, one worker. Python with numpy etc. runs through
`uv run -q --no-project --python 3.12 --with …`, written `UV --with …` below.

| Script | Checks | Runtime | Output |
|---|---|---|---|
| `qd/`: `python3 -I -B verify_cut.py` | balanced bisection, cut 1542, internal degree ≥ 9 (cut + viol ≥ 1546) | 0.1 s | `verify_cut.out`, `PASS_EXPQD_CUT` |
| `qe/`: `UV --with numpy python -I -B hilb.py` | Hilbert functions 512 / 708 / 772 / 708 | 0.3 s | `hilb.out`, `PASS_EXPQE_HILB` |
| `qe/`: `UV --with numpy python -I -B verify_sched.py` | bit-exact `schedule784` simulation: dim W = 772 (Fourier and monodromy) | 0.7 s | `verify_sched.out`, `PASS_EXPQE_SCHED` |
| `qe/`: `UV --with numpy python -I -B colrep_bound.py` | single-column replacement ≤ 781 for all 32 columns | 2.0 s | `colrep_bound.out`, `PASS_EXPQE_COLREP` |
| `qu/`: `UV --with numpy --with ortools python -I -B divcheck.py` | \|R\| divisible by \|H\| on random block orders | 0.5 s | `divcheck.out`, `PASS_EXPQU_DIV` |
| `qk/`: `UV --with numpy python -I -B qk_singer.py 2 1 0` / `4 2 0` / `6 3 0` | dim W = 9 / 772 / 57360 = S(q) | 0.1 s / 0.5 s / 34 s | `qk_singer_q{4,16,64}.out`, `PASS_EXPQK_SINGER` |
| `qi/`: `UV --with numpy --with scipy python -I -B flow.py 3 1 0 8 0` | q = 8: 8 information sets route integrally | 1.7 s | `flow_q8.out`, `PASS_EXPQI_FLOW` |
| `qv/`: `python3 -I -B table.py` | GQ T2*(O) table with the spectral cut charge | 0.1 s | `table.out`, `PASS_EXPQV_TABLE` |
| `ny/`: `UV --with numpy python -I -B check_telescope.py` | GRS telescoping identity on BH scenario A (imports `grs.py`, `bh_build.py`, `bh_data.py`; writes `bh_graph.json`) | 105 s | `check_telescope.out`, `PASS_EXPNY_TELESCOPE` |
| `qr/`: `UV --with numpy python -I -B check_lemma.py` | §5 algebra: Z = (I − K)Ã, E_ii spectrum, tr K², triangular determinant | 0.1 s | `check_lemma.out`, `PASS_EXPQR_LEMMA` |
| `qr/`: `UV --with numpy --with python-flint python -I -B families.py 0 same 2` | nullity 0 for every family at m = 2, seed 0 (one slice of the 142 runs in `fam_summary.txt`) | 6.1 s | `families_g0_m2.out`, `PASS_EXPQR_FAMILIES` |
| `qr/`: `UV --with numpy --with python-flint python -I -B ray.py same 2 F1_Droot_r1` | ray test, F1 family, m = 2: trivial gcd | 21 s | `ray_F1_m2.out`, `PASS_EXPQR_RAY` |
| `qr/`: `UV --with numpy --with python-flint python -I -B cutrank.py` | F1 cut-aligned rank drops, all 0 (writes `cutrank.json`) | 14 s | `cutrank.out`, `PASS_EXPQR_CUTRANK` |
| `ra/`: `python3 -I -B check_flow.py q16.json cert_flow_q16.json.gz` (and `s784.json cert_flow_s784.json.gz`) | stdlib exact-Fraction check of the rational multiflow certificates: MC_m ≥ ⌈mF⌉ for m = 1..5, MC_5 = E_int exactly | 0.7 s each | `check_flow_{q16,s784}.out`, `PASS_EXPRA_FLOW` |
| `l6/`: `python3 -I -B estimates.py` | expL6 planning arithmetic (imports `../../peel-dead/hyp_drop.py`) | 5.2 s | `estimates.out`, `PASS_EXPL6_ESTIMATES` |

Data (`ra/`): `q16.json`, `s784.json` (expRA instances) and the gzipped certificates `cert_flow_*.json.gz` (about 0.65 MB each).
Data (`data/`): `schedule784.json` (expNI, unchanged); `singsa_59_1.json`, `kah1544_ordsa1.json`, `pg1534_ordsa1.json` (expQD
layouts, converted from pickle to JSON). The full GRS census (all m, all seeds, both modes) and the remaining ray runs use the same
scripts with other arguments (`families.py gseed mode m-list`, `ray.py mode m family`; about 40 min per ray run at m = 8).

Not included: the SA, SDP and KaHIP searches of expQD, the code searches of expQE (`t2.py`, `s1.py`, `s4.py`, `s5*.py`,
`boot_sa.py`), the layout searches of expQU, the ½Z frame anneals and metric searches of expQI, the GQ dilation-block solver of
expQV, and expNJ's gadget census. Their results are quoted above.

Prepared with Claude Code (Claude Opus 5.5) assistance, using many sub-agents. It records other scouts' unreviewed local results;
nothing here has been independently reviewed.
