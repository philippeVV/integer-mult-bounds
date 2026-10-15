# PR2044 outer-layer audit (head 322949a5)

Review evidence; no κ claim. `kappa_exact.py` re-derives the outer κ from `new/SPEC.json`.
`knockout.py` solves the exact knockout LP: sup κ over (ε, β, ξ) with one layer's row reverted at a time; HiGHS picks the active set, then every constraint is re-checked in Fractions.
Results are in `knockout.json`, and the `log-*.txt` files are the PR's own outer checks run directly, since `run_outer.py` stops at the CRLF byte compare on LF hosts.
Prepared with substantial Claude Code assistance.
