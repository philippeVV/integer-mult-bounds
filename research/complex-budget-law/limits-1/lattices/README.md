# Lattices (expXA, expXH): Barnes–Wall fails by valuation, the Lorentzian E_n series has γ_n = 6/(n−1) → 0, and γ·(h−1)/u is the real figure of merit

Status: scoped closure of the growing-decoder lane, author-checked, LOCAL-ONLY until published. No κ claim. Base: main
`0fdd7d82`. Multipliers are over κ = 0.0017; pricing ratios are against the E8 P112 control in the same formulas.

## Verdict

| Item | Value |
|---|---|
| Question | is there a family of masked decoders whose γ (cost per port) → 0 as the window grows, and is it a step change? |
| Barnes–Wall (valuation lemma) | under parity labels the Gram decoder W/w is valid **iff every nonzero cosine has one 2-adic valuation**. E8 (cosines ±½) passes; **BW16 mixes ¼ and ½ and is invalid** (302,400 required pairs with N = ±½); BW32 and Leech inherit the obstruction. With labels compatible with G/4, h ≥ 53 and γ = 1.21. Deepening the 2-adic denominator is the wrong direction |
| Lorentzian E_n (LZD_n) | height-≤ 1 roots of I_{n,1} in k^⊥, k = (3; 1ⁿ): V = C(n+1, 3), decoder W/2 exact; n = 9 is P120, n ≥ 10 hyperbolic (signature (n−1, 1)). **C = n(n+1), γ_n = 6/(n−1) → 0** (exact for n ≤ 14, span method to n = 61). Below 4/7 from n = 12 |
| Frankl–Wilson families | k-subsets with the polynomial decoder P(W) = Π_{j<k/2}(W − 2j)/(k − 2j): γ exponentially small in h at k ≈ n/2 (FW_{48,8}: h = 47, γ = 0.024), with 10⁷–10⁸ ports |
| Why it is not a step change | the template pays per window: the label fraction is k = u/(h−1), so **the figure of merit is γ/k = γ·(h−1)/u**, not γ. E8 P112: 4.57. LZD: exactly 6 for every odd n (worse). FW beats 4.57 only at h ≥ 31, where the window term eats the gain |
| u = 1 lane ceiling | the free-decoder row (γ = 0) at h = 9: **×2.80 selected-word, ×5.02 data floor** vs E8; ×0.56 selected at h = 31. In every template regime the falling γ is eaten by m ∝ h |
| Fixed window (expXH step 1) | E_n roots ∩ (C + 2Z^{n,1}), n ≤ 16, heights ≤ 3–5, about 400 codes, plus E8 ⊕ ⟨−4⟩^ρ: **the largest valid set at h ≤ 10 is always 120 lines, an E8 copy**. Stop test fires: the u = 1 decoder lane closes at ×2.8 |
| Hermitian ports (tensor law, proved) | with the one-valuation condition over any O with u-dimensional residue labels (Z[i], Z[ω], unramified degree u): **B_C = D·(B_R ⊗ T)·Dᵀ**, so h = u·h_R, every term rank multiplies by u, γ_C = γ of the real pattern, and k = u/(u·h_R − 1) ≈ 1/h_R independent of u. Max valid sets over Z[i], Z[ω] equal the real counts (4, 10, 20, 36, 63 for n = 3…7); Gaussian / Eisenstein E8 keep only 12 / 4 of 120 lines |
| u-pricing | κ(u, h_R, γ) ≤ κ(1, h_R, γ) in both regimes: at h_R = 9, γ = 0, selected ×2.80 (u = 1), ×2.20 (u = 2), ×1.60 (u = 8) |
| Growing family in the Cesàro law | a_j = ν_j/Λ_j decreases, so a lattice family is worth its best member, never more (`../single-letter/`) |
| Multiplier | **×1.0** realised (decoders only, no unit or word built). Ceilings: ×2.80 selected / ×5.02 data floor (u = 1), ≤ ×2.20 / ×4.75 for any Hermitian u = 2 system; outside the template (peel, μ = 0) ×2.81 in the tilt root, regime unrealised |

## Changed hypotheses

1. (expXA) Growth without deepening denominators: keep the single denominator 2 (so W/2 stays exact under parity labels) and
   grow the **signature**. Root sets with inner products in {0, ±1} are bounded by ADE only in positive-definite signature; in
   I_{n,1} the E_n series continues past E8. This is a Frankl–Wilson-type modular-intersection construction: F₂ sees
   |T∩T′| mod 2 (rank ≈ n), Q sees a low-degree polynomial (rank n), the Q/F₂ min-rank separation that expQL Thm A requires.
2. (expXH) Freeze the window instead: multiplicity at fixed h (step 1), or u-dimensional labels from a Hermitian structure
   (step 2), hoping for γ → 0 at fixed k.

## Proof sketches

**Valuation lemma.** Integral Gram W = cG, labels B = (W + J) mod 2. B_ab = 1 exactly when W_ab is even, and such pairs need
N_ab = 0. So W/w is valid ⇔ every nonzero off-diagonal W_ab is odd ⇔ all nonzero cosines share one 2-adic valuation v and
v₂(norm) = v + 1. ∎ BW_{2n} ⊇ {(u, u)} inherits BW16's cosines; Leech has {0, ±¼, ±½}.

**γ_n = 6/(n−1).** Labels (x mod 2, 1) ∈ F₂^{n+2} under the dot form have form-rank n (odd n) or n+1 (even n); on the support of
coordinate term c (roots with x_c ≠ 0) the label span has form-rank exactly n (computed for n ≤ 30, 41, 61; **a general proof
is an open bookkeeping item**). So C = (n+1)·n, V = (n+1)n(n−1)/6. ∎ Lower bounds do not force γ ≥ c: the char-2 trace bound
gives rank ≥ V/h ≈ n²/6 for 2-integral decoders (LZD's rank is n, using the ±½ entries on free odd cycles); the PSD local ϑ
bound gives 11.95 at n = 11, above the indefinite W/2's rank 11, so the Lorentzian sign is used.

**Tensor law.** Absorb the unit c_a into v_a (diagonal unitary gauge: cost and validity unchanged). The label pairing is
Tr(λμ̄z_ab), z_ab = H_ab + 1 mod 2; validity forces z_ab ≡ 0 when H_ab ≠ 0, and z_ab ≡ 1 when H_ab = 0 or a = b. So the
uV × uV Gram is B_R ⊗ T with T = (Tr λ_iμ̄_j) invertible (Z[i]: I₂; Z[ω]: [[0,1],[1,0]]); Kronecker ranks multiply blockwise. ∎

**u-pricing.** With m = u + s(u·h_R − 1) and ν = u(1 − s/(s−1)·γ), ν/m → (1 − cγ)/(s·h_R): the u = 1 value with window
h_R − 1 replaced by h_R − 1/u, which is larger. ∎

## Data

| system | h | γ | γ·(h−1)/u | lp, data floor (vs E8) | selected word (vs E8) |
|---|---|---|---|---|---|
| E8 P112 (control) | 9 | 4/7 | 4.57 | ×1.00 (absolute ×4.05) | ×1.00 |
| LZD_13 | 13 | 0.500 | 6.0 | ×0.90 | ×0.74 |
| LZD_19 | 19 | 0.333 | 6.0 | ×1.18 | ×0.65 |
| FW_{32,8} | 31 | 0.107 | 3.22 | ×1.20 | ×0.49 |
| FW_{48,8} | 47 | 0.024 | 1.12 | ×0.91 | ×0.33 |
| γ = 0, h = 9 (free decoder) | 9 | 0 | 0 | ×5.02 | **×2.80** |
| Gaussian E8 (u = 2) | 6 | 1.00 | 2.50 | — (V = 12) | — |

## Next step

Bookkeeping only: retire "G/2^j on BW_{2^k} / Leech" by the valuation lemma; record LZD_n as the first proven γ → 0 family with
merit law γ·(h−1) ≡ 6; prove the per-term form-rank claim. A step change needs labels that are not tensor products (PR ∈
{1, …, u−1}) with a decoder other than Gram/2; no candidate exists. **Stop:** no structured system at h ≤ 17 with u ≥ 2 has a
2-adic decoder of rank < uV/h. Open micro-question (capped at the same ×2.8): can a complex decoder beat 4/7 on E8 labels at
u = 1 (an LM search for rank-7 masked decoders on P120)?

## Files and reproduction

Runtimes under a 4G cap, one worker.

| Script (run from its subdirectory) | Checks | Runtime | Output |
|---|---|---|---|
| `xa/`: `uv run -q --no-project --python 3.12 --with numpy python -I -B lor.py` | LZ/LZD and FW, exact decoder validity (Fractions), C = n(n+1) for LZD, n ≤ 14 | 3 s | `lor.out`, `PASS_EXPXA_LOR` |
| `xa/`: `uv run … --with numpy python -I -B scale.py` | LZD exact to n = 61 (γ = 6/(n−1)), FW_{n,k} to n = 128 | 1 s | `scale.out`, `PASS_EXPXA_SCALE` |
| `xa/`: `uv run … --with numpy python -I -B rankW.py` | rank W = n, signature (n−1, 1), all roots in k^⊥ | < 1 s | `rankW.out`, `PASS_EXPXA_RANKW` |
| `xa/`: `uv run … --with numpy python -I -B bw16.py` | BW16: 302,400 invalid pairs under parity labels; B₀ labels h₀ = 191, C = 2608; clique ≥ 53 | 1 s | `bw16.out`, `PASS_EXPXA_BW16` |
| `xa/`: `uv run … --with numpy python -I -B price.py` | pricing table; free-decoder row ×5.02 / ×2.80 vs E8 | < 1 s | `price.out`, `PASS_EXPXA_PRICE` |
| `xa/`: `uv run … --with numpy --with cvxpy --with clarabel python -I -B theta_lzd.py 9 11` | PSD local ϑ bound (8.00 at n = 9, 11.95 at n = 11) | 14 s | `theta_lzd.out`, `PASS_EXPXA_THETA` |
| `xh/`: `uv run … --with numpy --with ortools python -I -B part1.py 12 3 e8pad 60` | E8-padded code at n = 12: max valid set 120 (CP-SAT optimal) | < 1 s | `e8pad12_h3.out`, `PASS_EXPXH_PART1` |
| `xh/`: `uv run … --with numpy --with ortools python -I -B part1.py 12 3 rand 10 5` | five random [13, 8] codes: no valid set above 120 at h ≤ 10 | 16 s | `rand12_h3_quick.out`, `PASS_EXPXH_PART1` |
| `xh/`: `uv run … --with numpy --with ortools python -I -B part2.py` | Gaussian / Eisenstein E8: max valid 12 / 4 (optimal), tensor law h_C = 2h_R, term ranks double | < 1 s | `part2.out`, `PASS_EXPXH_PART2` |
| `xh/`: `uv run … --with numpy --with ortools python -I -B part2b.py gauss n 4 k`, n = 3…7; `eis` n = 5, 7 | Hermitian Lorentzian E_n: every valid port is a unit multiple of a real root; max valid = 4, 10, 20, 36, 63 | 0–15 s each | `part2b_*.out`, `PASS_EXPXH_PART2B` |
| `xh/`: `uv run … --with numpy python -I -B price2.py` | u-pricing; κ(u) ≤ κ(1) at every (h_R, γ) | < 1 s | `price2.out`, `PASS_EXPXH_PRICE2` |

Dependencies inside this branch: `xa/price.py` imports `kappa.py` from `../../../closures-3/large-ports-trace-bound/` (expNM),
`framework_ceiling.py` from `../../../framework-ceiling/` (expM1) and `hyp_drop.py` from `../../../closures-3/peel-dead/` (expQH);
`xh/price2.py` imports `../xa/price.py`. Not replayed (scouting outputs kept): `xa/fixedh.py` (weak random greedy at fixed h,
`fixedh.out`), and the larger expXH searches `blocks12/16.out`, `rand12/16.out`, `h5*_*.out` (heights ≤ 5, local search,
not certified) together with `e8pad.py`, `e8pad_ls.py`.

Sources: expXA, expXH. Prepared with Claude Code (Claude Opus 5.5) assistance. Pricing figures are planning ceilings, not κ
claims. Not independently reviewed.
