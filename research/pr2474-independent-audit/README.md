# Independent audit of PR #2474 (head 5298f86c)

Adversarial audit, no κ claim. **Verdict: the claim holds.** κ = 0.0017036763879936764890339…
is reproduced and independently re-derived. No soundness defect was found in the new 637←31 alias,
the 429-helper owner fill or the pricing. On Linux the package's release verification fails on line
endings, and passes once they are normalised.

The scripts take `PKG=<checkout>/research/gmx-t74-borrow43-handoff`, and `WORK=<verify work dir>`
for `finish_verify.py`.

## Independent checks (own code)

- `trace.py`, `indep.py`, `ctl.py`: **alias 637←31.**
  - **Edit scope.** In both words the edit is exactly the 7 declared B-gate renames 637→31 (B 1044/1045/1046/1058/1629/1630/1632). There are no frame-only changes and the 50 old warps are unchanged.
  - **The new warp.** It sits at B1041 on lender 31, with G = PR2436's warp on 31. It maps [388,75,41,26] onto span [257,129,120,30], which equals 637's old entrance span (120 = 102⊕30).
  - **Retirement and cold ends.** 31's last original use is B574. 637 had no cold end, and 31 ends at FULL.
  - **Value replay.** A randomized non-commutative replay (random 3×3 matrices mod 2⁶¹−1 per warp, independent dual symbols, cleanup derived independently) passes forward and dual on both words.
  - **Counts.** The frame chain gives profile {1:2363,2:632,3:326,4:120,5:33,7:4,8:4}, inventory {4:33,5:8,7:62,8:86,9:240} and 429 helpers. Against PR2462 the gain is exactly one width-5 helper and one paid 5-climb.
  - **Mutants.** Sign flip, doubling, dropped warp, flipped write, dropped last read and early warp are all rejected, with the cleanup re-derived.
- `fill.py`: **owner fill.**
  - Widths come from the shipped words.
  - ORDER is a permutation of the 666 pinned calls.
  - There are 285,714 rows. In each call the helper set is exactly range(429), with distinct families and distinct bins.
  - Every role lies inside its window and within its bin's capacity, and bins are filled as gap-free prefixes.
  - Idle is {1:21857,3:7308,5:35,48:25752}, and the deficit identity holds (23680).
- `kappa.py B.json`: **κ arithmetic.** local×666 + idle reproduces the global profile (rank mass 4,816,512). A 90-digit Decimal root of Σ n_r (r/73)^(1−s) = 66304 matches the claimed κ to 79 digits. The merged `ceiling.py` bound is κ ≤ 0.00172850.

## Defects

1. **CRLF breaks Linux reproduction.** This is high impact for reproducibility and has no impact on κ.
   - `git ls-files --eol` shows 66 `i/crlf` and 2 `i/mixed` files.
   - `outer/run_outer.py:14` compares raw bytes, so it fails with `regenerated artifact differs from committed: LEDGER.json`. The files are identical modulo CR.
   - The `provider_registry.{gauge44,extra}.artifact_sha256` pins for LOCAL-CHECKS.json / RESIDUAL-ADAPTERS.json in `new/PLAN.json` are hashes of the CRLF files. On Linux the regenerated natives hash to b7d12b86… / bbfe9033… / dd735bb6…; converting them to CRLF gives the pinned 3dad2d64… / 4b2348f0… / 09f5320f….
   - So `verify.py:116`, and then `code/bind_providers_portable.py:54`, fail. As a result `verify_release.py` cannot pass on Linux.
   - With line endings normalised, everything passes: `run_outer_crlf.py` gives PASS_OUTER, and `finish_verify.py` gives PASS strict-t74 / geometry-controls / PASS_STRICT_MIXED_T74_EXACT_FILL 0.0017007788112916789…, plus the binding checks of `verify_release.py:65-153`.
2. **Verification mutates the package.**
   - `new/rebuild_borrow.py` (run_audit, kernel_controls) and `borrow_controls.py` rewrite 12 committed files in `new/borrow/receipts/`, which include nondeterministic `elapsed_seconds` fields.
   - A second run then fails the manifest and entry-hash checks.
   - `verify_release.py:154` re-hashes only MANIFEST.json, so it does not catch this.
3. **The alias-637 tamper controls are vacuous** (`new/rebuild_borrow.py:96-99`).
   - The sign-flip and doubling controls do not re-derive the cleanup, unlike the handoff controls at lines 107-108. They are rejected by "explicit cleanup differs …". The dropped-warp control is rejected by "warp operator id".
   - Re-run with a re-derived cleanup, they are rejected for real reasons ("target map failed 32", "frame not monotone"). So there is no soundness gap.
4. **Minor.**
   - `lowfour_controls.py:35-38` still indexes `warps[-4:]` (now hits lender 55, not 33).
   - The credits still say "430-helper".
   - `verify_release.py:128-135` repeats the same gate three times, and line 144 says PR2322 where PR2363 is meant.
   - Dependencies are undeclared: numpy, gmpy2, sympy, mpmath, psutil.
5. **"653 does not stack"** is argued in prose and not checked by code, but it is structurally confirmed. 533's cold end is restored by `restore[56]: 533 += -x73`, which needs a clean 73. In PR2436's FOUR word, which has 653←73, 533 indeed has no `helper_final`.

Prepared with substantial Claude Code assistance. A sub-agent did the audit; the main session
re-ran `fill.py` and separately confirmed the `run_outer.py` CRLF failure.
