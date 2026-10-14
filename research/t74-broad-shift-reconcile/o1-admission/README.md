# O1: independent replay of the inherited T74/#1319 supplier admission

**Scope.** Review evidence only; no κ claim. Covers the T74 supplier used by PR2044 (strict
κ = 0.0016769381110224…). The outer layers behind PR2044's 0.0016797550 are **not** covered.

- `t74_final_binding.py STRICT GEOMETRY NATIVE_GAUGE44 NATIVE_EXTRA PARENT_KAPPA OUT.json` is the T74
  analogue of the cross-receipt final binding that `verify_t16.py` performs and PR1319's `replay.py`
  omits. It checks role, call and chart equalities, the owner and controller audit hashes inside the
  invoice, every strict_check, 47 positive outer margins, provider content binding, the template
  and profile hashes, and the control counts (114/12/8, 6 geometry, controller 7 + 8 per provider).
- `t74_semantic_regen.py` compares regenerated OWNERS/ORDER/DISPATCH/PROFILE/TEMPLATE/CONSTRUCTION/PLAN
  by semantic (decompressed) hash rather than container bytes.

| receipt | target | result |
|---|---|---|
| FB-1319 | PR1319 ancestor strict | PASS |
| FB-1625 | PR1625 strict | PASS |
| FB-2044 | PR2044 strict (verify.py: all 11 stages PASS) | PASS |
| FB-r37 | 37-replica identity regression through T74 code | PASS; reproduces PR1297 κ 0.001676854978651125 and the pinned REGRESSION37 controls |
| SEM-* | 1319 and 2044, Python 3.12 and 3.14 | PASS (semantic) |

**Negative controls.**
- Swapping the gauge44 and extra natives fails provider content binding.
- A parent κ that the candidate does not improve on fails the binding.

**Gaps found.**
1. PLAN pins `assignment_sha256` and `owners_file_sha256` to gzip container bytes. On Python 3.12 the
   gzip OS byte is 03, versus ff on 3.13+, so the container check fails while the payload is
   identical. PR2044 patches byte 9, but a different deflate build is untested.
2. `replay.py` never asserts the controller mutation count (23 = 7 + 8×2), and its RESULT label says 15.
3. GEOMETRY-CONTROLS.json records no hash of the strict files it mutated.
4. No hosted CI and no Python 3.12 maintainer run.
