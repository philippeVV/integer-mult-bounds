# Support-local global H call at PR2044 parameters

`sim_deferred.py OUT.json` runs an exact bookkeeping model: support-limited exposure, deferred phases, row padding and role split, and an adversarial dirty-control supplier. It compares the result with H^⊗g. `thresholds.py` computes the onset (log₂ p) of each "eventually" side condition at ε = 1−2e-30010, β = 1−6e-30010. Review evidence only.
