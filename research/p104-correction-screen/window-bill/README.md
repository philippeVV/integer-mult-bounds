# P104 side correction: global reachability windows and a favorable complete bill

FINDING, no κ claim. Card `p104-side-correction`, priority revision 4, base main `0fdd7d82`.
Model: the X-only correction model of PR2353/PR2393/PR2445. New helpers only climb to larger
frames, read only X, and write Y only by addition. Template as in #2082/#2083: m=73, W=832, N=9,
centre 56, idle mass 4344 per call.

## Results

1. **The line-only acquisition premise holds for Y0/Y1 but not globally.** In PR2393, each needed
   source of Y0 and Y1 is readable only on its initial line, which forces 48 separate lifetimes.
   For **48 of the 104 targets** that is false: **336 of the 4,992 correction entries** have a later
   source frame, of dimension 3–6, inside the target's final frame H_t. Up to 6 needed sources pass
   through one such frame, and one target has 12 sources readable on shared frames. So the
   per-target 48-lifetime injection does not carry over to those targets. The global 208 floor
   (#2353/#2383) is unaffected.
2. **Storage-window floor.** Every useful X read happens before cut 1104, and no target leaves the
   zero frame before cut 2526. A zero-frame register holds no source information. So at every cut in
   1104..2525 all correction information must sit in live new helpers. Because I−G/2 is
   nonsingular, **at least 104 new helpers must be live simultaneously** throughout that window, on
   top of the channel's 72 helper roles. This is a register-occupancy obligation for hosting.
   **Boundary:** at cut 1104 one original register (G6, inside H_11) still carries usable source
   information; the floor is stated for the X-only model.
3. **Favorable complete bill.** A lifetime arriving at a target needs at least three steps with
   ranks {7,1,1}: in from below, into H_t, and up to FULL. This follows from subadditivity and
   concavity of the per-step cost.
   - With 208 such lifetimes the channel+correction ceiling is a < **0.0025174**, so the optimistic
     total does not miss the target and the route is not closed.
   - The maximum number of lifetimes that still beats PR #2474 (κ 0.001703676387993676):

     | per-lifetime steps | max lifetimes |
     |---|---:|
     | {1,1,7} | 482 (idle-mass bound) |
     | {1,1,1,6} | 464 |
     | {1,1,1,1,5} | 418 |
     | nine rank-1 steps | 322 (312 at the bit cap) |

     That is about 3.1–4.6 lifetimes per target. The best published sharing is 39 per target
     (#2445).
   - Under the older rule, where added ranks earn no coverage credit, the allowance at #2474 is
     only 19.4 rank units.
4. **Value of a completed P104.**
   - The channel alone, with a free correction, has ceiling 0.0036634.
   - The bit side binds at **0.0017334270**, which is #2474's κ plus its `actual_bit_above_effective`
     margin, read from the #2474 ledger.
   - So even a perfect P104 lifts κ by at most **+1.746%** over #2474. Going further also needs a
     better bit supplier.

## Run

```sh
python3 -I -B check_p104_window_bill.py /path/to/repo-main \
  --ledger /path/to/pr2474/research/gmx-t74-borrow43-handoff/outer/LEDGER.json
# PASS_P104_WINDOW_BILL pairs=4992 window=1104..2525 simultaneous_min=104 floor_ceiling=0.002517 cap=0.001733
```

The checker uses the standard library only and takes about 7 s. It pins the witness, fixture and
EXPECTED.json by sha256 and uses the merged `ceiling.py`. Its controls:
- both window edges are tight;
- a corrupted source label breaks the support match;
- the idle cap of 52 is rejected;
- the 322 limit is exact (323 fails);
- the G6 boundary case is detected.

**Evidence status.** The finite facts are checked exactly, plus a short written argument for the
storage floor and the cheapest step shape. The work is author-checked: a sub-agent did it, and the
main session re-ran the checker. It has not been independently reviewed. This is not a
construction, a supplier or a new κ.

**Next.** Build a complete 104-target correction that uses the shared later source frames, with at
most about 482 lifetimes and at least 104 live over cuts 1104–2525. Alternatively, prove a lifetime
floor above 482 in this model. Then price hosting for the extra 104 simultaneous registers.

Prepared with substantial Claude Code assistance.
