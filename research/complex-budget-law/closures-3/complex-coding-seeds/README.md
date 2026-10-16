# Complex-side coding seeds (expL6 C1): closed, ×0

Status: scoped closure, author-checked, LOCAL-ONLY until published. No κ claim. Base: main `0fdd7d82`.
Multipliers are over κ = 0.0017.

## Verdict

| Item | Value |
|---|---|
| Question (expL6 C1) | is there a char-0 (Q(ζ)) Tanner coding network whose gap exists in characteristic 0 and dies in characteristic 2, with S/E above the a > 0 threshold? Only such an object could lift the binding complex side (κ = min(bit, complex)) |
| ζ / frequency codes on PG(2,q), q = 3…16 (expRB) | best char-0 S/E = 0.077 (q = 3) with MDS gates, falling like 1/(4q); a < 0 everywhere. **Chebotarev rigidity:** for prime N every DFT minor is nonzero in char 0, so S = k² exactly. **The char-2 gap is a mod-2 Chebotarev failure** (the same frequency set gives 82 at q = 8 and 708 at q = 16 mod 2, against 16 and 64 in char 0) |
| Walsh–Hadamard codes (expRF) | best char-0 S/E = 1/8 (biaffine planes, q = 4, 8, 16), a = −0.21; the excess reappears mod 2 in every one of ~830 rows (S2 ≥ S0) |
| Ceiling theorems (expRF) | **T1:** S0 ≤ Σ_e a_e b_e, so **S0/E ≤ 1/4 for flat-leverage codes on any layout** (girth-4 layouts, including every E8 incidence, can never pass). **T2/T3:** free-projection moments plus a Chebyshev–Markov atom bound give S0/E ≤ M_g with **M_g = 1/(2g) for even g**, exactly the threshold, so **a ≤ 0 at girth 4, 8, 12, 16** |
| What is left | only a sliver at **girth 6 (S/E ≤ 2/11, a ≤ 0.065) and girth 10 (a ≤ 0.017)**; no construction enters it (best 1/8) |
| Semicontinuity | rank only drops on reduction, so **a char-0 gap always reappears mod 2**; "char-0-only" can only mean "the reduced gates admit no valid orientation", never extra sessions |
| Multiplier | **×0** realised, and ×0 as a ceiling outside the girth-6/10 sliver |

**Consequence for the family record (`../char2-tanner-family/`): the complex side has no coding-seed route.** The bit-only
Tanner family therefore stays worth ×1.0 in κ, and the "complex-side counterpart" reopen trigger has no known candidate
mechanism (only non-flat-leverage codes with positively correlated leverages, close to coordinate codes, could beat T1).

## Changed hypothesis

The complex word is an E8 subspace template, capped at a ≈ 0.03 realistic. C1 replaces it by a coding network over Q(ζ) whose
gap lives in char 0 and dies in char 2 (the field-swapped twin of the char-2 Tanner route). Design target: by semicontinuity the
local codes must be 2-adically singular, and S/E > 1/6 is needed for a > 0 on girth-6 layouts (a = 1 − ln((E − S)/S)/ln b,
b = girth − 1).

## ζ / frequency codes (expRB)

**Model.** Singer model of PG(2,q): N = q² + q + 1, points and lines indexed by Z_N, incidence d ∈ D; for even q the Singer
matching is deleted (n = q), for odd q the full plane (n = q + 1). Gates are square [n, n/2] codes, and the Tanner dimension is
exact by Fourier: dim W = Σ_x dim(A ∩ diag(ξ^{xd}) B). Char-0 ranks are taken mod two split primes (≈ 1e9 and 1.9e9), which agree.

**Families.** Teichmüller-lifted homogeneous GRS (all line twists, including Frobenius lifts and complex conjugation; the ramified,
2-singular lifts are exactly C1's target class): S ≤ 16. Generic random codes: S = 0. Singer-frequency codes: S = k² exactly at
every prime N, and at composite N whenever the codes are MDS; excess only with non-MDS codes, and then at least as large mod 2.
Frequency hill-climb: best char-0 10/27/28/36/64 at q = 4/7/9/11/16, char-2 values of the same sets ≥ those (10/35/28/36/273).
Extended QR codes over Q(√p*): S ≤ 7.

**Mechanism.** On an abelian-group-equivariant layout over Q(ζ_N), Tanner excess equals the vanishing of DFT minors. In char 0 they
vanish only through subgroup structure (composite N), which on Sidon difference sets produced only non-MDS gates; every vanishing
persists mod 2. **ζ-structured codes are maximally rigid in char 0.**

**q = 3 bound (Tao's uncertainty theorem).** For prime N and any equivariant MDS A, B: f(c) = det[A; Λ_c B] has Fourier support in
the k-subset sums of D′; if those are distinct mod N, f ≢ 0, and Tao's theorem bounds #{c : A ∩ Λ_c B ≠ 0} ≤ C(n, k) − 1; dim = k at
most once. At q = 3 (D Sidon) this gives S ≤ 7 < E/6 = 8.7 for **every** equivariant MDS realisation. For larger q the bound is
vacuous and the table is the evidence.

| q | N | E | best char-0 S (MDS) | S/E | a (b = 5) | best char-2 S (MDS) | S/E |
|---|---|---|---|---|---|---|---|
| 3 | 13 | 52 | 4 | 0.077 | −0.54 | 6 | 0.115 |
| 4 | 21 | 84 | 4 | 0.048 | −0.86 | 9 (Hermitian) | 0.107 |
| 8 | 73 | 584 | 16 | 0.027 | −1.22 | **82** (same set mod 2) | 0.140 |
| 16 | 273 | 4368 | 64 | 0.015 | −1.62 | **708** (same set mod 2) | 0.162 |

(Full table, q = 3…16: `rb/table.out`.)

## Walsh–Hadamard codes and the ceiling theorems (expRF)

WH codes A = span{χ_s ∘ λ_v : s ∈ T}, |T| = n/2, are singular mod 2 (every ±1 character reduces to 1) and uncertainty-tight on
subgroup cosets: the opposite of ζ codes.

**T1 (leverage bound).** S0 = dim(U_P ∩ U_L) ≤ tr(Π_P Π_L) = Σ_e a_e b_e (in a simple graph the two projections share only diagonal
entries). Codes spanned by orthogonal unimodular characters have flat leverage a_e = 1/2, so **S0 ≤ E/4 on any layout**. It needs a
positive-definite inner product, so it fails mod 2, which is exactly where the char-2 gaps live.

**T2 (moments).** S0 ≤ tr((Π_PΠ_L)^j); closed walks shorter than the girth are tree-like, and with coordinate-transitive local codes
the moments equal those of two free projections of trace ½: μ_j = C(2j, j)/(2·4^j) (checked numerically on W(3)).

**T3 (Chebyshev–Markov atom bound, `rf/markov.py`).** S0/E = ν({1}) ≤ M_g = min{∫p dν_free : deg p ≤ g − 1, p ≥ 0 on [0,1], p(1) = 1}.

| girth | b | M_g | threshold 1/(2g) | a_max |
|---|---|---|---|---|
| 4 | 3 | 0.2500 | 0.2500 | 0 |
| 6 | 5 | 0.1818 (= 2/11) | 0.1667 | +0.065 |
| 8 | 7 | 0.1250 | 0.1250 | 0 |
| 10 | 9 | 0.1034 | 0.1000 | +0.017 |
| 12 | 11 | 0.0833 | 0.0833 | 0 |
| 16 | 15 | 0.0625 | 0.0625 | 0 |

The coding conversion's break-even point is exactly the Christoffel mass of the free law: flat half-rate codes can at best match
tree behaviour, which gives a = 0.

**Computations (scouting; `rf/table.out`).** Layouts PG(2,3), PG(2,4) − M, PG(2,7), PG(2,8) − M, PG(2,16) − M, biaffine(4, 8, 16),
W(3), W(4) − M, W(7), W(8) − M; natural, affine and random labellings; subgroup, weight and RM-like frequency sets. Best S0/E = 1/8
(biaffine, all q); S2 ≥ S0 in every row. n = 4 is closed by elementary means (every [4,2] WH code is a coset-pair code, so
S ≤ E/girth). The translation-equivariant biaffine family has an exact Fourier formula, and its hill-climb never beats q³/8
(`rf/equiv.py`). **E8 incidences** (120 lines vs 135 frames, either pairing; 240 roots vs 135 frames) **all have girth 4**, so T1
closes them without computing.

## Next step

None for C1. Optional paper checks: T2's free-product claim for non-transitive flat codes at j ≥ 3, and whether hexagon terms lower
the girth-6 constant 2/11 below 1/6 (if not, the sliver stays a documented, unattained ceiling).

## Files and reproduction

Runtimes under a 4G cap, one worker.

| Script (run from its subdirectory) | Checks | Runtime | Output |
|---|---|---|---|
| `rb/`: `uv run -q --no-project --python 3.12 --with numpy --with galois python -I -B c1.py 3` | the q = 3 family sweep (char 0 mod P and char-2 reductions): best char-0 S = 4 < E/6 (writes `c1_q3_full.json`); other q run the same way, about 15 min in all | 4.5 s | `c1_q3.out`, `PASS_EXPRB_C1` |
| `rb/`: `python3 -I -B table.py` | per-q digest of the full scouting sweep (`c1_all.log`): every char-0 S/E < 1/6, char 2 never below char 0 | 0.1 s | `table.out`, `PASS_EXPRB_TABLE` |
| `rf/`: `uv run -q --no-project --python 3.12 --with numpy --with scipy python -I -B markov.py` | the M_g table: M_g = 1/(2g) for even g, 2/11 at girth 6 | 0.6 s | `markov.out`, `PASS_EXPRF_MARKOV` |
| `rf/`: `uv run -q --no-project --python 3.12 --with numpy --with galois python -I -B equiv.py 8 3000 1` | exact equivariant biaffine formula and hill-climb at q = 8: best S0/E = 1/8 < 1/6 | 5.2 s | `equiv_q8.out`, `PASS_EXPRF_EQUIV` |
| `rf/`: `python3 -I -B table.py` | digest of the scouting sweep logs (`*.log`, copied): every S0/E < 1/6, S2 ≥ S0 in every row | 0.1 s | `table.out`, `PASS_EXPRF_TABLE` |

Not included: expRB's frequency hill-climb (`freqopt.py`, about 20 min) and expRF's sweep driver and libraries (`run.py`, `wh.py`,
`geom.py`, `moments.py`, `e8inc.py`); their logs are digested by the two `table.py` scripts.

Sources: expRB, expRF. Prepared with Claude Code (Claude Opus 5.5) assistance. T1–T3 are elementary but have not been independently
reviewed.
