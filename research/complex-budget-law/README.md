# The complex budget as a deficit/log-mass ratio

*Analysis note; no new κ and no supplier. Base: origin/main `0fdd7d82`. The selected record is #1094: κ = 0.001672660197, with the complex side binding at a = 0.0016754627. The task is self-selected synthesis work, closest card `complex-masked-decoder`, scope changed: a structural budget law plus a lever map rather than a decoder.*

## 0. The answer in one line

$$
a^\* \;\approx\; \frac{\nu}{\Lambda},\qquad
\nu=\frac{D}{W}=u\Bigl(1-\tfrac{s}{s-1}\,\gamma\Bigr),\qquad
\Lambda=\frac1W\sum_r n_r\,r\ln\frac mr .
$$

- **ν** is the net number of unpaid coordinates per family. Here u is the port dimension, s = 9 is the number of calls per replica, and γ = C/(uV) is the decoder's copy charge per unpaid port coordinate. For the selected record, ν = 1 − (9/8)(4/7) = **5/14** exactly.
- **Λ** is the paid mass weighted by log-price, per family. Each paid coordinate costs ln(m/r), where r is the rank of the child that advances it. For the selected record, Λ = **212.51**.

The resulting estimate ν/Λ = 0.0016806 is 0.3% above the true root. The two-sided bounds below enclose the true root within ±0.06%.

In the ideal relaxation, κ can be at most about 23× higher. That gap factors as **×2.8 from the decoder (ν/u)** and **×8.2 from child size (Λ)**. The local word supplies 97% of Λ, and the helper registers alone supply 65%.

## 1. Exact identity and two-sided bounds

Fix a histogram n_r with 0 < r < m. Define the mass M = Σ r n_r, the deficit D = mW − M > 0, and the mass distribution π_r = r n_r/M. Let L_r = ln(m/r), ℓ̄ = E_π L, and Λ = Mℓ̄/W.

**Theorem 1 (tilt identity).** Ψ(a) = W⁻¹Σ n_r (r/m)^{1−a} = 1 if and only if
$$\mathbb E_\pi e^{aL}=1+\delta,\qquad \delta=D/M .$$
Equivalently, D = Σ_r r n_r (e^{aL_r} − 1). The achievable saving a\* is the unique root.

*Proof.* We have (r/m)^{1−a} = (r/m)e^{aL_r}, so Σ n_r (r/m)^{1−a} = (M/m)E_π e^{aL}. This equals W exactly when E_π e^{aL} = mW/M = 1 + δ. The function φ(a) = Σ r n_r (e^{aL_r} − 1) − D is convex and increasing, and φ(0) = −D < 0. ∎

**Corollary 2 (sandwich).** With L_max = ln m,
$$
\frac{\ln(1+\delta L_{\max}/\bar\ell)}{L_{\max}}\;\le\;a^\*\;\le\;\frac{\ln(1+\delta)}{\bar\ell}\;\le\;\frac{D}{W\Lambda}=\frac{\nu}{\Lambda},
\qquad a^\*\ \ge\ \frac{\nu}{\Lambda}\,m^{-\nu/\Lambda}.
$$

*Proof.*
- Jensen gives E e^{aL} ≥ e^{aℓ̄}, which is the upper bound.
- e^{aL} lies below its chord on [0, L_max], which gives the left bound.
- e^x − 1 ≥ x gives the linear upper bound.
- e^x − 1 ≤ xe^x gives the last bound.

Second order: a\* = (ν/Λ)(1 − (ν/Λ)Λ₂/(2Λ) + O(a²)), with Λ₂ = W⁻¹Σ r n_r L_r². ∎

**Check on the selected global profile.** The profile is {1:349604, 2:60192, 3:62585, 4:20448, 5:11808, 7:7, 8:1152, 48:4866} with W = 14336. All values below are exact interval arithmetic.

| quantity | value |
|---|---|
| root a\* | 0.001675462676808 (the certified endpoint lies in this bracket) |
| upper bound, Jensen / linear | 0.00167647 / 0.00168059 |
| lower bound, chord / exponential | 0.00167456 / 0.00166852 |
| second-order estimate | 0.00167544 (relative error −1.2·10⁻⁵) |
| ℓ̄; effective rank r_g = m e^{−ℓ̄} | 2.9254; **3.92** |

**Exchange rates per call.** Write d = D/N = 35.56 and Λ_call = 21,157.
- One deficit unit is worth **2.81%** of a.
- One log-unit is worth 0.0047% of a. One deficit unit therefore equals about 595 log-units.
- An extra rank-1 child that displaces idle mass costs ln(R/1) = 3.9 log-units, or 0.019%.
- The same child charged on top of coverage costs 2.81%.

This 150× split is why the overhead convention (idle-displacing versus on-top) decides most "batching" screens.

## 2. The numerator: coverage turns into decoder compression

PR1012's coverage identity says paid mass = (m − u)W + C·N_calls. ctoth's #1587 count fixes families per call at w = (s − 1)V/s. Therefore
$$\nu=\frac{uW-CN}{W}=u-\frac{s}{s-1}\frac{C}{V}=u\Bigl(1-\tfrac98\gamma\Bigr).$$

- **Selected complex side:** C = 64, V = 112, γ = 4/7, ν = 5/14.
- **Bit side, same law:** γ = 63/111 gives ν = 321/888, matching the bit deficit T(2Kv − L·cst).

Two consequences follow.
- **Break-even is γ < 8/9.** A naive private query per port coordinate (γ = 1) always loses exactly u/9 per family. Only a compressing decoder, like E8's quarter-term Gram decoder, makes ν positive.
- **Decoder floors translate directly into ν.** The eight-term all-pairs floor (γ ≥ 7/13, PR968/1077) gives ν ≤ 41/104, which is at most +10% on a.

**Heuristic: changing the lattice does not help.** Assume C = n² (n terms of rank n), m − 1 = n(n+1) and s = n + 1. Then a·ℓ̄ ≈ 1/(n(n+1)) − 1/V. Values:

| lattice | value |
|---|---|
| E8 punctured | 0.00496 |
| E7 | 0.00198 |
| 10-dim kissing configuration | 0.00517 |
| Barnes–Wall (16-dim) | 0.00321 |
| Leech | 0.00166 |

E8 is essentially optimal. This is a heuristic scaling, not a theorem.

## 3. The denominator is pinned by frame diversity

**Per-call decomposition of Λ** (exact, from the shipped `LOCAL-PROGRAMS.json.gz`; forward and backward agree):

| class | mass | log-mass | share of Λ |
|---|---:|---:|---:|
| helper climbs (480 registers; widths 9×240, 8×122, 7×30, 5×32, 4×56) | 3730 | 13802 | **65.2%** |
| data climbs (224 registers × 8) | 1792 | 6439 | 30.4% |
| centre (8 × rank 8) | 64 | 142 | 0.7% |
| idle remainder | 1646 | 774 | 3.7% |

The 240 width-9 helpers alone carry 8028 log-units, which is **38% of Λ**.

**Lemma 3 (subadditivity).** For σ = 1 − a ∈ (0, 1], Σ_i (r_i/m)^σ ≥ (Σ_i r_i/m)^σ. Two consequences:
- *(One-shot.)* If every register climbs its whole width w in one child, the moment only falls. So the "one-shot" profile gives a ceiling for every frame schedule and packing of the same register inventory.
- *(Fixed climb counts.)* Fix c climbs summing to w. The sum Σ (r_j/m)^σ is concave, so it is minimised at an extreme point: (c − 1) unit steps plus one step of w − c + 1. That profile gives a ceiling for every re-ranking with unchanged climb counts.

**Measured facts.**
- **All 3,486 local climbs are touch-separated:** a scalar touches the register between any two of its consecutive climbs.
- Every helper's entry dimension equals the dimension of its first-touch gate.
- The width-9 helpers are touched mostly at dimension 0 (the precondition, 8727 touches) and dimension 9 (the cleanup, 7963 touches). Even their middle-phase touches span dimensions 0/1 through 8.
- A late-entry probe asked whether the zero-frame prefix commutes to the first middle touch. It does for only 4 of 240 helpers, and then only up to dimension 1.

| denominator ceilings | κ multiplier |
|---|---:|
| re-ranking with the same climb counts | **×1.022** |
| one-shot data registers / helper registers / both | ×1.13 / ×1.36 / **×1.62** |
| all local children at the local dimension h = 9 | ×1.71 |
| remove the 240 width-9 helpers / all helper climbs / plus data one-shot | ×1.53 / **×2.46** / ×3.47 |

So Λ is fixed by *how many distinct frames each register is touched at*, to within about 2%. Breaking this needs fewer or narrower helpers, not better packing or ranking.

**A question on cross-register "batching" (#1779 and the diagonal-batch-commutation package).** A child is C^{⊗rf} acting on r slots of *one* register. Rank-1 transforms on k different registers form a direct sum I_k ⊗ C^{⊗f}: still rank 1 at k times the volume, not C^{⊗kf}. If that reading is right, the +39.5% figure needs a tensor embedding that should be made explicit. I have not audited those packages in full.

## 4. Lever map

The multipliers are ideal relaxation ceilings on a; each changes one lever unless stated. "×2" means κ = 0.003345 (a = 0.0033565).

| lever | invariant | ideal | what κ×2 needs | proven caps and obstructions | status |
|---|---|---:|---|---|---|
| idle packing | idle Λ (3.7%) | +0.87% | impossible | complex-local-word-ceiling | **capped** |
| frame re-ranking, fixed climbs | Λ_local | ×1.022 | impossible | Lemma 3 plus touch-separation (here) | **capped** |
| frame coarsening, fixed inventory | climb counts | ×1.62 | impossible alone; with it, needs C ≤ 56 (u=1) or ≤ 155 (u=2) | Lemma 3 one-shot; singleton 4→6 cap; coalescing-subset | capped at ×1.62; mostly unavailable |
| cross-register batching | — | 0 in transform semantics, if ⊕ ≠ ⊗ applies | — | question above | unclear |
| lattice / template (n, s) | ν/(n(n+1)) | E8 ≈ optimal | — | ctoth U ≤ v(N−d); heuristic in §2 | **capped** |
| line decoder, u = 1 | γ | ×2.80 at C = 0 | **C ≤ 29** (γ ≤ 0.26) | 8-term all-pairs γ ≥ 7/13; E-graph rank-8 rigidity, so ≥ 64; fixed-reader 9th term ≥ 64; singleton rank-9 ≥ 64; A8 nine-term ≥ 72; masked co-star k ≤ 6 impossible | essentially capped (only masked ≥9-term decoders with new readers are open) |
| plane ports, u = 2 | γ₂ | ×5.57 at C = 0 | **C ≤ 128** (γ₂ ≤ 4/7, the incumbent's own ratio) | 2-integral ≥ 224; half-integral ≥ 55; ρ≤2 ≥ 150; ρ≤3 ≥ 96; all-ρ≥7 ≥ 175 (a23); rank ≥ α+2; scalar-mask obstruction for parity lifts | **open, narrow** (window 96–128 for low-ρ, mixed-denominator terms) |
| helper inventory | helper Λ (65%) | ×2.46 | **helper climb mass ≤ 16% of today** | none, but touch-forced in today's flag schedule | **open, needs a new local algebra** |
| 0.01 | both | — | u=2 and C ≤ 111 with helper climbs removed; or C ≤ 136 also with data one-shot; or u ≥ 2.15 at C = 0 | no single lever reaches it: u=2 alone caps at 0.0094 | needs ν **and** Λ |

The bit side follows the same law and sits only 3.5% above the complex side. Every step change on the complex side must be mirrored there.

**Pricing conventions differ between plane models.** #2126 prices its fit at 112 over 56 *shared* matched planes (PR1825). The plane-port bounds on discussion #1587 (2-integral ≥ 224, etc.) assume 112 *distinct* per-port planes. These are different models, so the two numbers are not directly comparable. The u=2 thresholds below use the per-port-plane model.

## 5. Recommendation

**Most promising open lever: the coverage numerator through plane ports (u = 2).** Three reasons:
- It is the only single lever whose ideal headroom (×5.6) far exceeds 2.
- It is necessary for 0.01.
- Its κ×2 target is precise: γ₂ ≤ 4/7, which is E8's line compression ratio carried over to planes.

The Λ levers are either capped (×1.02 to ×1.62) or need a qualitatively different local algebra with about 84% less helper climbing.

The exchange rate gives a partial-upgrade form of the question. Upgrading a set S of ports earns (8/9)|S| deficit units per call. So κ doubles iff |S| − (9/8)ΔC ≥ 40. At ΔC ≈ 0, 40 ports (u ≈ 1.354) already suffice.

**First experiment (small and exact): a single orthogonal-frame upgrade.**
- Take an orthogonal frame S of 8 live ports and one common r with r·q_a = 0 and r·r = 1. Upgrade those 8 ports to planes ⟨q_a, r⟩ and keep the other 104 as lines.
- Find the minimum copy charge of a decoder for the resulting masked target.
- **Pass iff ΔC ≤ 7**, because 8 × 8/9 = 7.11.

**Step 0 is done here** (`frame_upgrade_probe.py`, self-built E8 model):
- P112 contains 1260 such frames.
- For each frame and r: |T| is 48, 52 or 56; there are 192–224 newly forced zeros; the correction block on each of S×T and T×S has rank_Q 8 and F2 pairing rank 4.
- So additive rank-1 repair costs 16 > 7 and fails.
- Under a matrix-coefficient term convention, the repair costs 4 + 4 = 8 > 7.11. That narrowly fails too.

The live question is whether re-choosing the incumbent eight terms can absorb the two rank-8 corrections for at most 7 extra units. Those eight terms are rigid at rank 8, so at least 9 terms are needed. This is a fixed-size exact problem: 112 ports, one frame, a support dictionary. It should be run on the selected labels from `DATA-PHASES.json` under the pricing convention settled above. A pass makes κ×2 a matter of repeating it on about 5 frames with shared corrections. A fail with a proof turns the u=2 route into a whole-assignment question with a sharper floor.

## Reproduce

Python 3.10+, standard library only. Inputs are sha256-pinned in the script.

```sh
python3 -I -B research/complex-budget-law/verify_theory.py          # ~2 min, PASS_DEFICIT_LOGMASS_LAW (output in full.out)
python3 -I -B research/complex-budget-law/verify_theory.py --fast   # ~40 s
python3 -I -B research/complex-budget-law/frame_upgrade_probe.py    # ~4 s
```

Controls reproduced exactly:
- ceiling 0.001690054109;
- centre-free ×2.80;
- plane-free 0.009410;
- the coverage identity 72W + 64N;
- the selected endpoint inside the root bracket.


Prepared by Philippe Van Velzen with substantial Claude Code (Claude Opus 5.5) assistance; no independent human review is claimed. Apache-2.0.
