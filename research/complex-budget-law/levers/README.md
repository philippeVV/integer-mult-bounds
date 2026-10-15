# Λ levers on the selected word: what they mean, and why hyperplane ends are blocked

Follow-up to `../README.md` (lever map: one-shot ×1.62, helper inventory ×2.46). Scope: the
selected #1094 k8 word (`research/k8-rank7-composition/sources/K8-SLOT-WORD.json.gz`). The
end-of-B frames are computed with the repo's `hpe_audit_and_transport.py` (complex-571-cold14).

## What the two levers are

- **Climbs are pure cost.** Frames only grow, and every scalar gate runs at one common frame. If
  you store the untransformed value, a climb F → F′ does nothing to the content. It is simply a paid
  child of rank dim F′ − dim F.
- **Helpers.** Helpers are borrowed dirty streams hosted on idle families. In the selected word,
  each helper's climb count is the number of distinct frames it is touched at, minus one. So a helper
  costs the log-mass of its frame chain between its first and last touch, and a helper used at a
  single frame costs nothing. Helpers carry about 3,730 of the 5,586 local mass per call.
  - **Helper inventory ×2.46** means every helper is frame-local.
  - **One-shot ×1.62** means every register is touched only at its entry and exit frames.
- **One-shot status.** Forbidden for fixed scalar words: `complex-x4-frame-bound`,
  `frame-domain-ceiling`, `complex-k8-singleton-frame-cap`, and on the bit side PR827 /
  `bit589-frame-optimum`. With a rewritten program it is capped at ×1.62.
- **Helper inventory status.** Open, with no general theorem. On this word it is blocked as
  described below.

The estimates below are first-order lever prices (ideal-idle ceiling) from a scout pass. Their
scripts are not shipped, so treat them as estimates:
- ×1.05 would need about 8% less helper climb mass.
- Ending all 402 full-ended helpers on a hyperplane would give ×1.066.
- Delaying the entry of the 240 width-9 helpers by one dimension would give ×1.047. That is blocked
  by the zero-frame precondition, since only 4 of the 240 commute.

## Hyperplane ends are blocked on this word (`census.py`, `cert.py`)

**Census.** For each of the 402 helpers that end at FULL, `census.py` finds the minimum, over all 511
hyperplanes y^⊥, of the number of cone offenders. These are the obstacles from the audit's necessary
cone pre-check: row support of P⁻¹, column support of P, and h's own end-of-B frame (P⁻¹[h][h] = ±1,
so h is always in its own cone).
- **No helper is at distance 1.** The distribution is 42 at distance 2, 11 at 4, 1 at 5, 3 at 6,
  22 at 7, 103 at 8, and the rest higher (`out/census-base.out`).
- **Positive control.** All 78 helpers that already end low pass (32 at 4D, 46 at 8D).
- **The 34 apparent single-offender helpers.** Leaving out h's own frame gives 34 helpers that look
  like single-offender near misses (`singles.json`). Every one of them is blocked by its own frame.

**Certificate.** For those 34, `cert.py` certifies impossibility in the following rewrite class:
- (i) any re-association of h's gates within their exact commutation windows;
- (ii) then **any number** of compensating gates h += c·g′ after h's last gate. Each must sit at a
  frame inside U, and the frame raise it forces on g′ must respect every later pin:
  - target port hyperplanes;
  - the terminal frames of already low-ended helpers;
  - their early-cleanup gates.

The test is whether row h of P⁻¹, restricted to registers ending outside U, lies in the span of what
the admissible compensations can cancel. It is computed mod p. All coefficients are dyadic, so
non-membership mod p implies non-membership over Q.

| pins kept | helpers with a repairable hyperplane |
|---|---|
| all (`out/cert.out`) | **0 / 34** |
| target port hyperplanes only (`out/cert--onlytargets.out`) | 0 / 34 |
| low-end pins only (`out/cert--onlyhpe.out`) | 32 / 34 |
| none, positive control (`out/cert--nopins.out`) | 32 / 34 |

So the target port pins alone block every case. This is the complex-side analogue of PR827's bit
cleanup-frame obstruction, restricted to this rewrite class. Two of the 34 (506 and 518) end in a read
into a target at its port-level frame, so no window remains for any repair.

**PR2044's T74 words.** The same census on gauge44-word-566 and extra-word-566 finds the same
structure: 402 full-ended helpers, none at distance 1, 42 at distance 2, and 0/34 repairable
(`out/*t74*`).

**Composition note.** The PR1673 outer layer applied to the 0.0016733419 T16 map (complex
a ≈ 0.0016761467) gives κ ≈ 0.0016761. That is below PR2044. T16-style packing does not stack on
T74, whose idle remainder is an exact fill `{48: 22748}`.

## Run

```sh
W=research/k8-rank7-composition/sources/K8-SLOT-WORD.json.gz
python3 census.py $W /tmp/census.json            # exact Fractions, a few minutes
python3 cert.py $W /tmp/census.json.pkl          # 0 of 34
python3 cert.py $W /tmp/census.json.pkl --nopins # 32 of 34 (control)
```

Needs numpy, run from this directory. Both outputs were re-run from this location in the main session
and are byte-identical to `out/`.

Prepared with substantial Claude Code assistance. Sub-agents did the scout and the certificate.
