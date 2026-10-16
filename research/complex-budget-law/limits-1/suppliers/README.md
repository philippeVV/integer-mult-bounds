# Suppliers (expXF, expXJ, expXN, expXS): movement vs arithmetic, no exponent-2 arithmetic for integers, and odd-characteristic Hermitian codes

Status: theory plus small exact checks, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017.

## Verdict

| Item | Value |
|---|---|
| Question (moonshot "bit side only") | can integer multiplication run on finite-field (char 2) suppliers alone, so that the complex side stops binding? |
| What the two suppliers are | the "bit" supplier is the **movement** supplier (address interchange, XOR on data bits): char 2 in **every** recipe. The "complex" supplier is the **arithmetic** supplier (H0 butterflies on the coefficient ring R): char R. **Every candidate recipe keeps both**; changing the outer only changes R |
| Lemma A (proved) | **no coefficient group of exponent 2 can compute integer products**: any transform scheme computes a Z-bilinear map, and c(2a, b) = 2c(a, b) maps to 0 when every group has exponent 2, so only c mod 2 is computed; lcm of the exponents must exceed max\|c\| ≈ N·2^{2b}. Additive char-2 FFTs compute carry-less products |
| Lemma B′ (proved) | the radix-ℓ DFT kernel has \|det\|² = ℓ^ℓ, so a radix-ℓ layer needs gcd(ℓ, char R) = 1: in the arithmetic supplier, frame denominators and gate characteristic are always coprime |
| Lemma B (Greenberg) / C (Ax–Lefschetz) | Z/ℓ^K with K unbounded, and every outer whose characteristic grows (single large NTT prime, Fermat rings, Fürer, DKSS, HvdH modular), are in the **char-0 class**: a bounded seed with a gap at infinitely many characteristics has the same gap in char 0 |
| The one escape: MSP | the multi-small-prime outer (≈ 3 lg n/ln lg n odd primes p ≈ 4.2 lg n after Bluestein) costs only constants (CRT volume ×1.00) and puts the arithmetic supplier in char p, with seeds that must grow with p |
| Law | **κ ≤ min(a_mov, a_arith)**. With today's bit side, MSP alone is worth **≤ ×1.01** (the bit-endpoint cap) |
| Conjecture M (expXF) | "p-integral frames make char-p-only gaps unusable" — **refuted as a mechanism** by expXJ: phase-frame cost lives in F₂, the field of the labels, not in the gate field (Lemma Y) |
| Phase frames (expXJ) | modulo free corrections a phase frame **is its polar form β ∈ Sym_m(F₂)**; edge cost c(Δβ) = rank_F₂ (+1 if alternating), endpoint Δβ = I |
| Small-m law (Theorem Ψ) | for m ≤ 3, **s ≥ Wm for every phase frame on every circuit realisable over any field** (char 0, char p, liftable or not). BH has s_min = 157m at m = 1, 2, 3. A frame family needs m ≥ 4 and non-commuting β |
| Odd-q Hermitian codes | S(q) = [(q+1)(q²+q+1) − √q(q²+1)]/4 on the full PG(2,q), q = p²: **166, 3449, 26434** at q = 9, 25, 49 (S/E = 0.182, 0.204, 0.216 → ¼). q = 9 is exactly Braverman–He's code |
| Multiplier | **×0** realised (no frames). Conditional on frames, schedule and soundness on **both** suppliers: ×97 at q = 49 unfused, ×142–187 as q → ∞; with today's bit side ≤ ×1.01 |

## Changed hypothesis

"Make the arithmetic run in characteristic 2 so the bit-side Tanner gaps apply to both suppliers." Lemma A kills it for
integers. The surviving variant moves the arithmetic supplier to odd characteristic p (MSP) and asks whether char-p-only
coding gaps can be used there.

## Lemmas (expXF)

**Lemma A.** Let c = 𝒞(𝒜a ∘ ℬb) with additive 𝒜, ℬ, 𝒞 into/out of groups G_i, ∘ bilinear, c ∈ Z^{2N−1} the convolution.
Then c(2a, b) = 2c(a, b); if every G_i has exponent 2 the left side is 𝒞(0∘ℬb) = 0, so the scheme computes c mod 2 only.
No carry pass can recover the lost information. ∎ (Side remark: for F₂[x] multiplication an all-char-2 recipe exists:
Schönhage's triadic transform.)

**Lemma B′.** |det DFT_ℓ|² = ℓ^ℓ, so the layer is invertible over R iff ℓ ∈ R^×. The movement supplier escapes only because its
frames are address permutations. ∎

**Lemma B.** A seed is a finite polynomial system (gates, rank conditions dim W ≥ k, soundness, schedule). Solvable mod ℓ^K for
all K ⇒ an ℓ-adic solution (Greenberg) ⇒ a char-0 gap. Char-2-only Hermitian codes are first-order obstructed (expQB), so
they already die over Z/4. ∎

**Lemma C.** "Gates exist on this incidence with dim W ≥ k and a sound schedule" is an existential sentence of bounded size;
true over F_P for infinitely many P ⇒ true in a char-0 ultraproduct ⇒ true over Q̄. Gaps specific to char P need seeds that
grow with P. ∎

**Exponent losses of outers** (`xf.py` D): Schönhage–Strassen/nested-Nussbaumer depth lg lg n with ×2 expansion drives
κ_eff → 0 for every a; Fürer's log* depth and single-level outers keep a. Genuine twiddle rows cap κ at a/(1+a) unfused.

## Phase frames reduce to Sym_m(F₂) (expXJ)

Upstream frame C_Q = H^{⊗m} diag(i^{Q(x)}) H^{⊗m}, Q a Z/4 quadratic form. Constant phases, 2[l·x] terms (translations X_l)
and Fourier shifts are free; each changes Q by 2·(linear) + const and preserves the polar form β(Q) = B mod 2. A kernel
(aI + bX_v)^{±1} adds vvᵀ to β, since [v·x] ≡ (ṽ·x̃)² mod 4. So **edge cost = c(β_head − β_tail)**, c(A) = rank_F₂A, or
rank + 1 if A ≠ 0 is alternating (brute-forced at m ≤ 3), and the endpoint condition is β_out(ρw) − β_in(w) = I.

- **Lemma X (exponent mismatch).** An additive ψ from phase data into a char-p space (p odd) is zero: 2β = 0 forces ψ = 0. The
  GRS row of T1/expNY (Z ≡ I mod p) has nothing to act on. Operator-level reduction works but gives only the concave screen
  rank(C_VC_U⁻¹ − I) = 2^m(1 − 2^{−ρ}), too weak.
- **Lemma Y (char-blindness).** s_min(incidence, ρ, m) = min over β of Σ c(·) after terminal elimination; it never mentions
  gate coefficients. So achievable savings depend on the incidence only; liftability does not enter.
- **GRS-F₂ (char-2 twin of T1).** For any char-2 gate assignment G, ψ = id is 1-Lipschitz, so s ≥ rank_{F_{2^k}}[T_G[j,k](I + β_j
  + β_k)]. If the incidence realises ρ in char 2 then Z = I and **s ≥ Wm for every phase frame**: the arithmetic supplier needs
  incidences that are *not* char-2 realisable (upstream's ½-coefficient network is one; a char-p-only incidence is another).
- **T2 twin.** Realisable over some field + terminal β_in simultaneously triangularisable over F₂ ⇒ s ≥ Wm (block-triangular Z_G,
  diagonal blocks of total generic rank W). Savings need non-commuting β.
- **Theorem Ψ.** For m ≤ 3, ψ(β) = diag(λ₁(β), …, λ_m(β)) with F₂-linear λ_i lifted to {0, 1} ⊂ K, λ_i(vvᵀ) = 1_{S_i}(v)
  (m = 2: S = {e₁}, {e₂}; m = 3: S = {e₁, 111}, {e₂, 110}, {e₃, 011}) is 1-Lipschitz for c on all pairs and has ψ(β + I) −
  ψ(β) = diag(±1), for every field K. Telescoping with the circuit's own K-realisation gives s ≥ Wm. ∎ (exhaustive check;
  at m = 4 no diagonal-linear ψ exists: four disjoint RM(2,4) supports of size ≥ 4 do not fit in 15 points.)
- **Control.** On expNJ's two unrealisable W = 6 incidences the model does express savings: s = 4m < 6m at m = 1, 2, 3.

## Odd-characteristic Hermitian Tanner codes (expXJ, `xj/odd_singer.py`)

Full PG(2,q), q = p², degree q + 1, square [q+1, (q+1)/2] gates: point code = values of degree-(q−1)/2 forms on the pencil
(doubly extended RS, MDS); line code = the same with x ↦ x^√q twisted coordinates (PΓL-canonical). Computed by the Singer/trace
model in F_{q³} with a character-block decomposition (wrap cocycle ρ = g^{(√q−1)(k−1)}), validated at q = 9 against dense
elimination (166, = BH's nondeterministic code).

| q (p) | E | S = dim W | S/E | a (b = 5) | κ = a/(1+a) (×) | ½·spectral bisection bound | slack |
|---|---|---|---|---|---|---|---|
| 9 (3) | 910 | 166 | 0.1824 | 0.068 | 0.064 (×37) | 159.25 | 6.75 |
| 25 (5) | 16926 | 3449 | 0.2038 | 0.153 | 0.133 (×78) | 3417.75 | 31.25 |
| 49 (7) | 122550 | 26434 (predicted blind) | 0.2157 | 0.198 | 0.165 (×97) | 26348.25 | 85.75 |
| ∞ | | | → ¼ | → 0.317 | 0.24 (×142); fused ≈ ×187 | | q^{3/2}/4 |

The closed form was fitted on q = 9, 25 and predicted q = 49 before that run. Compare char 2 on PG(2,q)−M (expQK):
q(q² − q^{3/2} + 1)/4. The schedule obstruction is the same as in char 2: an exact helper-free schedule needs a
near-spectrally-tight bisection. MSP fit: primes p ≡ 3 mod 4 make Z[i]/(p) = F_{p²} exactly this code's field; seed size
~p⁶ = polylog n (admissible by `../single-letter/` §overheads, and as Lemma C requires).

## Rev-5 obligation ledger for MSP

| obligation | status |
|---|---|
| char 2 impossible for the arithmetic layer | **proved** (Lemma A) |
| radix coprime to the gate characteristic | **proved** (B′) |
| Z/ℓ^K and growing-characteristic rings are char-0 class | **proved** modulo Greenberg and Ax |
| MSP outer: exact residue convolution after Bluestein, CRT volume ×1.00, conversion once | **proved by reading** (upstream skeleton; not written as rows) |
| phase-frame identities over F_p[i] | **proved** (polynomial in ½, i) |
| char-p seeds with a_p ≥ a uniformly in the prime window | **constructively upper-bounded** in dimension only (S(q) above); soundness, rigidity and schedule **unknown** |
| phase frames realising a saving for a char-p-only incidence | **unknown**; impossible at m ≤ 3 (Ψ), needs non-commuting β at m ≥ 4; expXN: nullity 0 for all 23 families at m = 4–8 (evidence) and a two-class theorem |
| hosting hundreds of characteristics on one fixed TM | **unknown, not priced** |
| movement side a_bit | unchanged; Problem X(F_{2^e}) schedule open (`../bootstrap/`) |

## expXN: phase frames on BH at m = 4–8, closed at ×1.0 (evidence + a two-class theorem)

- **Screen.** 23 terminal families × m = 4..8 × 3 random char-2 gate draws over GF(2^16). That is 345 exact ranks; every
  one equals Wm, so every GRS-F₂ nullity is 0. The families cover:
  - generic symmetric and alternating forms;
  - the expNY signed-polar twin;
  - symplectic and orthogonal orbits;
  - F₉ labels through fixed F₂-embeddings that break Ψ;
  - an m = 4 adversarial walk.
  12–21 families per m generate all of M_m(F₂), so the non-commuting case is genuinely tested.
- **Two-class theorem (proved).** For β = X on a session set S and Y off it:
  det Z = Π_μ det T_SS · det T_S′S′ · det(I − μ²M_S), with M_S = T_SS⁻¹T_SS′T_S′S′⁻¹T_S′S.
  If one draw has no nonzero root of χ(M_S) in F_{2^d} for d ≤ 15, the nullity is 0 for every m ≤ 15 and every X, Y.
  All 35 partitions tested pass (`xn/twoclass.log`).
- **Char-2 realisability of BH** is undecided: SAT on 157 sessions, about 2.6M clauses. Either answer gives the same
  verdict.
- **Verdict.** The finite-field arithmetic route (odd-p codes + H0 phase frames) closes at ×1.0 on BH. The q = 25
  Hermitian lift is not pursued, because its precondition (a family with positive nullity) failed.
- **Scope.** Evidence over the families listed, plus a proof for the two-class case. GL(m, F₂) address maps in frames
  are not covered.

## expXS: bit-side frames at m = 5 on PG(2,16)−M — the exact MC₅ pattern is infeasible

- **Lemma T (any field).** If D₁ + … + D₅ = I with every rank D_i ≤ 1, each D_i is an idempotent of trace 1. So along
  every 5-path of an exact session the frame trace rises by ±1 per step, and each exact session forces its own flag-level
  trace function.
- **Certificate** (`xs/check_cert.py`, stdlib): two sessions disagree by 2 (s784) or by 4 (q16), which no sign choice
  reconciles. At most 8 sessions (a star) can be exact simultaneously.
- **Robust bound (any frames, any m).** s ≥ W − |X| + MC_{m−1}. For every rank-≤1 frame family at m = 5, the saving is
  at most 8 + (4W − MC₄). So the metric "room 264 / 336" at m = 5 is a relaxation artefact.
- **Still open.** Frames using some rank-2 segments.

## Files and reproduction

`xs/check_cert.py inst.json trace_inst.json` (stdlib, < 1 s) prints CERT_OK for s784 and q16 (`xs/check_cert.out`). The `xn/` scripts are kept as scouting sources with logs; they are not replayed (they depend on the expNY/expN1 BH build).


| Script | Checks | Runtime | Output |
|---|---|---|---|
| `xf/`: `python3 -I -B xf.py` (stdlib) | A exponent obstruction; B kernel determinants; C MSP prime counts and CRT volume; D recursion-depth law; E rank mod p of a fixed integral seed | 7 s | `xf.out`, `PASS_EXPXF` |
| `xj/`: `python3 -I -B psi_check.py` | Theorem Ψ at m = 2, 3 (exhaustive) and the m = 4 non-existence of diagonal-linear ψ | 1 s | `psi_check.out`, `PASS_PSI` |
| `xj/`: `python3 -I -B control.py` | exhaustive min s = 4m < 6m on the two unrealisable W = 6 instances | < 1 s | `control.out`, `PASS_EXPXJ_CONTROL` |
| `xj/`: `uv run -q --no-project --python 3.12 --with numpy python -I -B odd_singer.py 3 dense` | S(9) = 166 by blocks and by dense elimination; closed form | 3 s | `sq9.out`, `PASS_EXPXJ_SQ9` |
| `xj/`: `uv run … --with numpy python -I -B odd_singer.py 5` | S(25) = 3449 = closed form | 2 s | `sq25.out`, `PASS_EXPXJ_SQ25` |
| `xj/`: `uv run … --with numpy python -I -B odd_singer.py 7` (slow) | S(49) = 26434 = closed form | 23 s | `sq49.out`, `PASS_EXPXJ_SQ49` |

Included but not replayed: `xj/beta.py` (polar-form tables; library), `anneal.py` with `an_m2_rank.log`, `an_m3_rank.log` (BH
upper bounds 314, 471), `solve.py` (CP-SAT, one worker), and `grs_families.log` (GRS-F₂ nullity 0 on six families at m = 2, 3;
its scripts depend on the expNY/expN1 BH builders, which are not copied). `xj/bh_graph.json` is a copy of expN1's BH scenario-A
incidence (735 segments, 157 sessions).

Sources: expXF, expXJ, expXN, expXS. Prepared with Claude Code (Claude Opus 5.5) assistance. Lemmas A, B′, X, Y, GRS-F₂,
the T2 twin and Theorem Ψ are proved here; B and C use standard theorems on bounded seed systems; S(q) is exact at three q and
the closed form is a fit with one blind prediction. Not independently reviewed.
