# P104 `I − G/2` correction: structure, floors and pricing screen

**Scoped screen; no κ claim and no supplier.** This follows the native channel
(research/p104-native-channel) and the existing floors #2082/#2127/#2136 (cut ≥ 96),
#2083/#2089 (terminal 104/208), and #2173/#2251 (shared space 96).

Run `python3 -I -B research/p104-correction-screen/verify.py`. It uses only the standard library,
takes about 55 s, and prints `PASS_EXPAZ_P104_CORRECTION_SCREEN`. Full output is in `RESULTS.json`.

## Exact structure of C = I − G/2 on the 104 ports
- **Spectrum and inverse.** G² = 26G, and C = −S/2, where S = G − 2I is the signed adjacency matrix
  of the 48-regular non-orthogonality graph. C² = 12I − 11C, the eigenvalues are 1 (×96) and −12 (×8),
  and C⁻¹ = I − G/24.
- **Smith form.** S = −2C has Smith form 1⁸·2⁸⁸·48⁸, and G has 1⁸. Over F2, F3, F5 and F13, S has
  rank 8, 96, 104 and 104.
- **What carries denominator 2.** Only an 8-dimensional part of C has denominator 2, the same
  part as G/2. The 96-dimensional remainder is integral.
- **Orbits.** The stabiliser of F0 ∪ C has orbits of size 24, 32 and 48 on the ports.
  Diagonal orbit blocks of C have full rank; off-diagonal blocks have rank 8.
- **Irreducible residual.** For any E of rank ≤ 8, rank(C − E) ≥ 96, with equality at E = −12P
  (P = G/26). After optimal reuse of the channel, the irreducible correction is the projector onto ker G.

## Floors with the native channel's frames kept fixed
- **Line-port pairing-rank model.** The support of C is exactly the off-diagonal set where the F2
  pairing is zero, so every entry is a free cross-addition. All the cost comes from timing and
  frame availability.
- **Timing.** Every source finishes climbing by event 1,533, and targets start at 2,526.
- **Carriers with no relays.** Each carrier covers at most 12 pairs, found by exact enumeration of
  frame chains. That forces at least 4,480 carriers, which is ruled out under any pricing.
- **Upward-only relays.** Each target needs two 7-dimensional source subspaces, so at least
  **208** new helpers are needed. This assumes the channel's own registers are not reused as relays.

## Pricing (`ceiling.py`, exact intervals)
| case | value |
|---|---|
| allowed added correction units per call, repo's "on-top" rule, at a = 0 / a_sel | 36.44 / 19.69, against floors of 96 or more: **ruled out** |
| free-correction ceiling (P104 channel) | a < 0.0036634; the bit side 0.00173343 binds first, so κ gain ≤ +3.46% |
| idealised minimal P104 word | 0.0064730 |
| helpers counted as coverage progress (#2152 accounting) | possible only with 208 ≤ new helpers per call ≤ 482 (idle mass binds) |

**Open route.** Correction carried by 208–482 new upward-only helpers hosted on idle coordinates.
That route still needs an explicit relay network, timing, hosting and precision.

Prepared by Philippe Van Velzen with substantial Claude Code (Claude Opus 5.5) assistance;
no independent human review is claimed. Apache-2.0.
