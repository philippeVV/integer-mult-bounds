# Pending (placeholder): the record column-operation ILP

Status: PLACEHOLDER. The job is still running at the time of writing; nothing here is final, nothing is replayed, and nothing
here enters the verdicts of the other lanes. No κ claim. Base: main `0fdd7d82`. (The full PG(2,16) schedule search, expXX,
was pending in an earlier draft of this package; it has finished and is now the lane `../full-pg16/`.)

## Record column-operation ILP — RUNNING

- Scope as recorded in expXO's next step: a frame-constrained column-operation solver, i.e. reads of helpers into each other
  (like #2857's `q3 += q2 − q0`), not only reads of k, formulated as an ILP over insertion slots with exact floor costs from
  the DP. Details, inputs and results to be filled in when the job reports.
- It works inside the current architecture (helper/word level) and is listed here only so that the round's open jobs are in
  one place. No sources are copied yet.

Prepared with Claude Code (Claude Opus 5.5) assistance. Not independently reviewed.
