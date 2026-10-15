# Joint query: all four banks lose a query (follow-up to the first-pair catalogue)

**Scoped obstruction and pricing; no κ claim.** This continues discussion #2074 and the parent
`../README.md` (none of the 252 first pairs avoids the loss).

## Results
- **Theorem A (four-bank erasure).** In the retained nine-call replica, each bank k = 0..3 has its own
  erased family. Take h ∈ ker B_{2k}, enter it on X_kU, and add the derived T/A components. Then
  every observation in the generous contract is zero (f-values, selected registers at every cut
  after the bank opens, final Y_U), while the selected input is nonzero. PR2174's pair-0 records
  are identically zero on banks 1–3.
  - `controller.py` replays the literal maps with exact rationals and noncommuting permutation routes.
    It prints `PASS_QQ16_FOUR_PAIR_ERASURE`.
  - Four mutation controls are each rejected (`ctl_*.log`): h ∉ ker, 1/12 in place of 1/11,
    the wrong pair map, and a dropped bridge.
- **Moments.** All 84 stage triples give moments (24, 28, 109/3, 1823/36), so every bank's defect is
  at least ⌈77N/5⌉ (`moments.py`). This holds under both the homomorphism and anti-homomorphism conventions.
- **Theorem B (perturbation rank).** If Σ rank(M_i' − M_i) < dim J(⊕V_k), a selected input is lost.
  The proof telescopes the observation map. Any repair, including paid reads, writes, helpers,
  parking or low-rank word changes, therefore needs **≥ 62 coordinates per owner per replica**.
- **Theorem C (reordering).** Every assignment of charts to stages still loses at least ⌈77N/5⌉ per bank.
- **Finite-model sufficiency** (`modelcheck.py`, evidence only). With 16 records per bank, nothing is lost
  under the generous contract. Under the honest final-only contract the same 64 records still lose
  [32, 16, 16, 0] (m = 3).

## Pricing (`price3.py`, pinned `ceiling.py`)
| scenario | complex ceiling vs selected |
|---|---|
| QQ16 with centre charge C = 64 (PR2079's convention) | at most +0.28% even with a free controller |
| QQ16 with C = 0 (coverage law; the cores are centre-free) | gain survives only for in-place or single-call storage; the bit side caps it at +3.45% κ |

**Which centre charge applies (C = 0 or C = 64) decides the route and needs maintainer confirmation.**

## Run
```sh
python3 controller.py ../../community-frontier1042-audit/selected-replay.tar.gz 5 1   # needs python3.12 + the package's deps
python3 moments.py
QQ16_LIB=1 python3 modelcheck.py            # QQ16_CONTRACT=final|finalf for the honest variants
```

Prepared by Philippe Van Velzen with substantial Claude Code (Claude Opus 5.5) assistance;
no independent human review is claimed. Apache-2.0.
