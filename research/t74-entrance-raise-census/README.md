# Entrance-raise census on PR2474's T74 words: no split-free raise exists

**Scope.** This is a FINDING with no κ claim. It looks at the borrowed words of PR #2474 (head `5298f86c`,
package `research/gmx-t74-borrow43-handoff`, `new/borrow/words/BORROW-WORD-43{base,extra}.json.gz`).
It tests the hypothesis in PR #2453 that split-free entrance raises are worth up to +4.3e-6.

**What a raise is.** An entering helper h has `start_extra` frame S. A raise replaces S by S+x, which
lowers h's width by 1. `audit_local.run_borrow` checks only that frames nest gate by gate, and a gate
frame is shared by both registers it touches. So after the raise, every gate at the old start frame
must contain x. That forces x into each partner's later frames, up to the point where the partner's
own frame already contains x. The minimal forward closure (`frames.py`, `Word.lift`) is therefore
forced for any frame re-solve. The raise is illegal if the closure reaches a register whose end
frame does not contain x: a target, whose end is the hyperplane ⊥ its port, or a cold-ended helper.
A partner that already sat at the old frame must now climb by 1 mid-run, which is a split.

**Control.** On the rebuilt pre-entrance chain state (`chainstates.py`, `t4.py`), each of PR2435's
24 raises closes in 3 events, and together they change exactly the 72 events PR2435 changed. The
same code recounts PR2474's local profile `{1:2363, 2:632, 3:326, 4:120, 5:33, 7:4, 8:4}` exactly.

**Census** (`enum1.py`, about 6 min; `an3.py`). The candidates are all 141 entering helpers × every
one-dimension raise, giving 30,851 per word. The base and extra words agree, since only 4 warp G
matrices differ.
- 30,833 candidates reach an end frame that does not contain x. 14,204 of these hit a target's port
  hyperplane or a cold end directly; 16,629 do so through warps.
- 10 fail at a non-target end frame.
- **8 are legal**, all width 5→4 helpers in the PR2436 alias region: 606, 635, 638, 639, 667, 669,
  670 and 671.
- Every one of the 8 needs partner splits, adding rank-1 climbs:
  - 606: {1:+9, 2:−5}
  - 639 and 671: {1:+11, 2:−6}
  - 638 and 670: {1:+15, 2:−8}
  - 635 and 667: {1:+19, 2:−12, 3:−4, 4:+4}
  - 669: {1:+20, 2:−11, 3:−1, 4:+1}
- Priced with PR2453's vendored PR2387 width-aware LP ceiling (`lp.py`, `price_feas.py`, `joint.py`),
  each one lowers the ceiling, by 5.7e-7 to 1.07e-6. Joint lifts also lower it.
- **Stacking and class lifts.** The forced requirements of a joint raise contain the union of the
  single raises' requirements, so an infeasible single raise stays infeasible in any stack.
  `enum2.py` and `an2.py` lift whole frame classes (all 52 seed classes × all x), and every one fails.

**Conclusion.** On PR2474's words no split-free entrance raise exists, and every legal single raise
loses. PR2453's +4.3e-6 assumed 12 raises costing only {3:−1, 2:+1} each, with no partner effects.
Here every candidate site either floods into a target port hyperplane or needs splits. **Not
explored:** turning non-entering helpers into entering ones, which changes the entrance contract
rather than only the frames.

## Run

```sh
export PKG=/path/to/pr2474/research/gmx-t74-borrow43-handoff   # PR #2474 head 5298f86c
python3 -B enum1.py      # writes enum1.pkl (30,851 per word)
python3 -B an3.py        # obstruction breakdown + the 8 legal raises
# pricing (optional) needs PR2453's package: export P2453=/path/to/pr2453/tree
python3 -B price_feas.py; python3 -B joint.py
```

Standard library only, except the optional pricing step. A sub-agent did the work. The main session
re-ran `enum1.py` from this directory (output byte-identical) and `an3.py`. It did not re-run the LP
pricing. The sign of each raise is clear from the added rank-1 climbs.

Prepared with substantial Claude Code assistance.
