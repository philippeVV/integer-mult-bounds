# The bit budget as a deficit/log-mass ratio (selected PR589 on PR555)

*Analysis note, companion to the complex-side note (#2278). No new κ, no supplier. Base origin/main `0fdd7d82`.
Selected record #1094: κ = 0.001672660197, complex a = 0.0016754627 binds, bit a = 0.001733427039589844234749.*

## 0. Answer

The bit primitive obeys the same law as the complex one, with almost the same numbers:

| | ν (deficit per family) | Λ (log-mass per family) | ν/Λ | true a |
|---|---|---:|---:|---:|
| bit (PR589) | 321/888 = 1 − (9/8)(63/111) | 207.92 | 0.0017386 | 0.0017334 |
| complex (#1094) | 5/14 = 1 − (9/8)(64/112) | 212.51 | 0.0016806 | 0.0016755 |

So any κ step must move **both** ν/Λ ratios by the same factor. On the bit side 54% of Λ is *hosting fragmentation*
(helpers carved out of idle windows), the line decoder is capped at ×1.12 by the proven 9/17 floor, and every u = 1 lever
that is still open is at most ×1.2–1.5. The cheapest route to a bit budget ≥ 0.0034 is the **same plane-port lever** as
on the complex side, and its threshold is the same exchange law: **|S| − (9/8)ΔC ≥ 38.5 per call** (≥ 17.6 for
a = 0.0025), i.e. C₂ ≤ 127 at full u = 2 (complex: C ≤ 128).

## 1. Exact identity and reproduction

Template (`global_k.py`, `assemble9.py`): m = 1 + 7L = 64, L = 9 stages, K = 4 pairs, T = 32 replicas, 288 local
invocations per call, W = 2KvT = 28,416 live families (all 126,432 helper residuals hosted, no banks). Global profile =
288 × local + REM.

**Theorem 1 (bit tilt identity, with the bad-residue surcharge).** With L_r = ln(64/r), N = number of children and
ε = 32·64²·N/10¹⁶ (pinned `moment.py`, `bad=True`),
Ψ(a) = 1 ⇔ Σ_r r n_r (e^{aL_r} − 1) + ε·64^a = D = 10,272. ε = 1.19·10⁻⁵ moves the root by 2·10⁻¹², so the complex
sandwich (Jensen / chord / linear) applies verbatim.

| quantity | value |
|---|---|
| record endpoint | own interval engine accepts c = 1733427039589844234749/10²⁴ and rejects c + 10⁻²⁴; pinned `moment.certify` proposes the same c |
| bounds: exp / chord ≤ a* ≤ Jensen / linear ν/Λ | 0.0017260 / 0.0017323 ≤ a* ≤ 0.0017336 / 0.0017386 (+0.30%) |
| second-order estimate | 0.00173341 (rel. error −1.15·10⁻⁵) |
| ℓ̄, effective rank 64·e^{−ℓ̄} | 3.267, **2.44** (complex 3.92) |
| control: eleven-stage price (m = 78, W = 8,880) | reproduced exactly; ν = 139/370 (+3.9%) but Λ = 218.9 (+5.3%), hence lower |

Exchange rates (first order): one copy unit per invocation = **2.80%** of a (matches bit-orthogonal-pair-decoders'
+2.8017%); one plane-upgraded port per invocation = 2.49%; one helper index (all 288 invocations) = 0.124%.
One copy unit is worth 22.7 helper indices.

## 2. Numerator

D = mW − mass = T(2Kv − L·cst) = W − 288·cst. Per family: ν = 1 − (L/2K)·γ with γ = cst/v = 63/111, exactly the complex
law with s = L = 9. The helper inventory does **not** enter D: helpers are idle-displacing (§3).

- Proven decoder floor (bit-orthogonal-pair-decoders, label decoders, ≤ 25 punctures): γ ≥ 9/17 ⇒ ν ≤ 55/136.
  Priced on the selected word: **×1.119** at most. Attained only at P102 (cst 54), where the hosting wall 4v = 408 < 439.
- cst needed at u = 1 with the shipped word: ≤ 47.2 for a = 0.0025, ≤ 28.8 for 0.0034 — both below the floor (58.8 at v = 111).

## 3. Denominator: hosting fragmentation

**Slot identity (derived, not read off).** Per pair the template leaves idle windows of 7Δ coordinates,
Δ ∈ {1,1,1,2,2,4,4,6,6}; a width-8 residual fits ⌊7Δ/8⌋ = Δ − 1 per window, so capacity = 4vL per replica (the wall
R ≤ 4v). R = 439 leaves 45 slots = 9 empty 42-windows per replica, and the remainder is forced:
REM = 32 × {7:666, 6:444, 4:444, 2:435, 42:9} — equal to the pinned `fastc.py` REM. Template idle mass =
288·8R + mass(REM) holds exactly, so **helpers add no mass and no deficit; they only re-price idle windows.**

Per-family partition of Λ = 207.92 (exact, from the replayed word):

| part | log-mass | share |
|---|---:|---:|
| data climbs (222 registers × 7 per invocation) | 50.94 | 24.5% |
| decoder copies (9 × rank 7) | 1.41 | 0.7% |
| template idle windows {7,14,28,42} (what R = 0 would pay) | 42.68 | 20.5% |
| **helper fragmentation** (430 helpers + 9 totals, 8 each, plus the 6/4/2 remainders, minus the windows they replace) | **112.90** | **54.3%** |

Helpers alone carry 60.3% of Λ (969 rank-1 and 1,041 rank-2 climbs per invocation; 4–6 climbs per helper). With R = 0,
a ≈ a₀/(1 − 0.543) = ×2.19 (exact ceiling ×2.185).

## 4. Lever map (exact ceilings; slot-exact hosting baseline; surcharge kept)

"×" is a/a₀. Targets: a = 0.0025 (×1.442), 0.0034 (×1.961); complex κ×2 is a = 0.0033565 (×1.936).

| lever | ideal × | needed for 0.0025 / 0.0034 | proven caps and obstructions | status |
|---|---:|---|---|---|
| idle packing | 1.128 (free concentration at rank 45) | impossible | PR1023: width-8 slot hosting is already optimal; recipe class closed (PR1487, open); rank ≤ 45 forced by row < 20161 (HD = 2) | **capped** |
| frame re-ranking, fixed climb counts | 1.037 | impossible | bit589-frame-optimum: exactly frame-optimal; r439 dirty-workspace frame bound | **capped** |
| one-shot data | 1.084 | impossible alone | B-block exits forced to rank-4 frames (PR1810) | capped |
| one-shot helpers / both | 1.348 / 1.505 | both one-shot clears 0.0025 only | cleanup forces all 439 helpers to full frame after B (bit-cleanup-frame-obstruction) | mostly unavailable |
| helper width 7 (all 439) | 1.234 | impossible alone | 433 forced to rank 8 after scatter, all after B (bit-width-seven-boundary); one constructed (helper 470) | open only with a changed prefix |
| line decoder γ | 1.119 at the 9/17 floor (2.771 at cst = 0) | cst ≤ 47.2 / ≤ 28.8 | γ ≥ 9/17 (label decoders); P111 optimum 63 (PR757/761); Parseval blocks cost 73 as scalar copies (PR1936) | **capped** |
| helper inventory | 2.185 (R = 0) | keep ≤ 43% / **≤ 9.4%** of helper climb mass | hosting wall; touch-forced frames in today's word | open, needs a new local algebra |
| plane ports (u = 2) | 3.8 at cst₂ = 63; 5.6 at 0 | **C₂ ≤ 146 / ≤ 127.5** per call; partial: net ports \|S\| − (9/8)ΔC ≥ 17.6 / 38.5 | additive common-r repair fails (step 0 below); Parseval plane decoder (PR1900) decodes the *line* target, not the masked plane target | **open** |
| uniform local rank (data+helpers) | — | ≥ 7 / ≥ 15 | — | — |
| universal cap (all mass at rank 45) | 9.3 (0.01608); 25.8 at cst = 0 | — | coverage law + halving guard | — |

Stacking all four line-port levers (one-shot data, width-7 one-shot helpers, 9/17 decoder) gives ×2.077 — enough for
0.0034, but only by combining four individually obstructed changes. Width-7 + 9/17 decoder gives ×1.381 (short of 0.0025);
adding one-shot data gives ×1.527.

## 5. Recommendation and first experiment

**Cheapest route: the shared plane-port lever.** The bit exchange law is the complex one: upgrading a set S of ports earns
(8/9)|S| deficit per invocation, each copy unit costs 1, so the bit side clears a = 0.0034 iff |S| − (9/8)ΔC ≥ 38.5
(17.6 for 0.0025), versus 40 for complex κ×2. Both sides need the *same kind* of masked E8 label decoder, so one mechanism
study serves both, and no Λ lever on either side has comparable headroom.

**Step 0 (done here, `probe_plane.py`, 25 s).** Rebuild the 120 lines from the 111 labels (B0 = (1+G) mod 2 has F2 rank 9).
For each root line r and each pairwise-odd set S of live ports in r^⊥, give S the planes ⟨q_a, r⟩. Every newly forced
zero lies in S×T ∪ T×S, and C = [[0,M],[Mᵀ,0]] has rank 2·rank(M). Over all 120 r and all 226,903 maximal sets
(|S| = 3…7), rank(M) ≥ |S| − 1, so the additive correction needs ≥ 1.6|S| rank-one terms, against a gain of 0.89|S|. Nullity is monotone under
subsets, so this covers every pairwise-odd common-r set, not only maximal ones.
Best gain − rank = −10/3. **Additive repair never pays**, matching the complex step 0 (16 > 7).

**Step 1 (proposed, small and exact).** Take the 11,950 most favourable sets (|S| = 5, rank(M) = 4, a one-dimensional
dependency among the T-rows). Ask for the minimum copy charge of a re-chosen P111 decoder for the masked target. Each
of these targets has 111 ports and 5 upgraded planes. Use the label-term catalogue (255 × 255 functional pairs), priced
by rational pairing rank on the plane-extended supports.
- **Pass iff ΔC ≤ 4**, because (8/9)·5 = 4.44.
- A pass gives the bit side its first net-positive port set, worth 5 − (9/8)ΔC net ports. a = 0.0025 then needs
  17.6 net ports in total (at least 4 disjoint repeats at ΔC = 0); 0.0034 needs 38.5.
- A fail with a min-rank certificate (for example, the masked completion needs ≥ 10 terms of cost ≥ 6 beyond the
  free entries) closes common-r upgrades on the bit side. The live question then becomes whole-assignment
  (per-port r_a or paired planes, as in PR1825).

## Reproduce

Inputs come from origin/main `0fdd7d82`, exported with `git archive` into `./main`, and are sha256-pinned. Standard library only.

```sh
python3 -I -B verify_bit.py      # ~90 s, PASS_BIT_BUDGET_LAW (output: full.out); --fast ~35 s
python3 -I -B probe_plane.py     # ~25 s, PASS_BIT_PLANE_UPGRADE_STEP0 (output: probe.out)
```

Controls reproduced exactly: the bit endpoint (own engine and pinned engine), the eleven-stage endpoint, the global
profile from 288 × replayed local + slot-derived REM, D = T(2Kv − L·cst) = 10,272, the PR1023 width-7 row (0.0021391622),
and the free-concentration row (0.0019559188). The plane-port ceilings here (cst₂ = 63: 0.00655–0.00676; cst₂ = 0:
0.00964–0.00995) sit below plane-port-copy-budgets' (0.00731–0.00767; 0.01076–0.01128) because REM is kept slot-exact
instead of ideally packed at rank 45.

Scope: ceilings are relaxations with the stated idle conventions, so a fail excludes and a pass is only a target.
Plane-port pricing of terms is not yet formally specified upstream; the step-1 design states its convention.

Prepared by Philippe Van Velzen with substantial Claude Code (Claude Opus 5.5) assistance; no independent human review is claimed. Apache-2.0.
