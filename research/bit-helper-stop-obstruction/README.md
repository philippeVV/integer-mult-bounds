# Bit helper-contract: no rank-1 frame stop is removable by a single gate

**Scoped exact obstruction. No κ change, and the paid profile is unchanged.**
Card `bit-helper-contract` (priority revision 4, chooser seed
`d818a9bc64e148d10e8f548dd5bc0123`). The object is PR1122's 10,652-ADD bit word,
`research/r439-abc-semantic-bridge/inputs/candidate.json.gz` (PR589 lineage).

## Question

Can a replacement H2/K boundary implication remove or merge a rank-1 frame stop
that only the syntactic `H2_D_frame0_targets_from_helpers` needs? A successful
merge would change the paid profile `{1:1323, 2:1342, 3:79, 4:44, 5:102, 6:11, 7:19}`
(mass 5,129).

## Result

1. **D-only replacements cannot move the profile.**
   - D gates act at the zero frame and never touch sources; the D-operand check
     excludes indices below v.
   - Every target and helper starts at the zero frame.
   - So any H2/K replacement that concerns D alone leaves all 661 frame paths
     and all 9 COPY ranks unchanged.
   - This explains why about 70 R439 packages (for example PR1937 and PR2033)
     report the same profile.
2. **Census of stops.** There are exactly 205 stops in 0→1→2 climbs, all on
   helpers.
   - Each one is a single A gate `h ^= x_s` at frame `{s}` with one source donor.
   - They come from 111 sources: 94 sources give two stops each, 17 give one.
3. **Exhaustive single-gate relocation finds 0 improvements** (`relocate.py`).
   - Exact-value donors: 1,037 insertion times and 704 legal frame choices.
   - Dirty donors with D/K recompensation (`--dirty`): 1,328 legal choices.
   - Frames outside the catalogue cannot help, because the per-step cost
     g(d) = d(64/d)^c is concave.
4. **Dual witnesses, one per stop** (`witness.py`).
   - **185 lateral blocks.** After `{s}`, x_s climbs to a different rank-2
     frame G, and two distinct rank-2 frames that contain `{s}` meet only in
     `{s}`. Example: helper 459, source 109, F2 `[4,109]`.
   - **20 descent blocks.** x_s still has a B target write at `{s}` after h's
     first read. Example: helper 222, source 0, F2 `[0,3]`. Source 0 is still
     written by B at time 2,910, but h222 is first read at 1,471.
5. **Replaying the 205 natural merges** (`replay_merges.py`).
   - Each merge moves the stop gate onto F2 just before h's first F2 touch.
   - Every merge passes the 661-column value endpoint (dirty inputs included)
     and the semantic D/K boundary.
   - Every merge fails pinned `xcheck_p` with `climb not nested`, and fails the
     COPY/ERASE `local_bill` with `frame climb`.
   - So the interface never blocks a merge. Source frame nesting does.

Before and after, the profile is `{1:1323, 2:1342, 3:79, 4:44, 5:102, 6:11, 7:19}`.
A legal merge would have changed it by {1:−2, 2:+1}, worth −0.002419 per stop in
Σ n_d·d·(64/d)^c, but no merge is reachable.

## Not covered

Each of these is an A/B rewrite, not an H2/K question:
- rewrites over several gates that change a source's own path, for example
  moving the 20 late B writes off `{s}`, or giving sibling helpers a shared
  rank-2 frame;
- changing which helper sums the A tree forms;
- gate reorderings that go beyond inserting one gate.

## Reproduce

From the repository root: `sh research/bit-helper-stop-obstruction/reproduce.sh`
(Python 3.10+, standard library only, a few minutes). It imports the pinned
`research/e8-bit-p111-tuned/code/xcheck_p.py`. Expected output:
- `stops 205`
- `improving 0` (exact and dirty)
- `sources with stops 111`
- witness split `20 A / 185 B`
- 205 `reject:climb not nested` merge replays

The frame-only optimum (`frameopt.py` from `bit589-frame-optimum`, run on this
word) reports `improving_components 0`.

## Open next questions

- A multi-gate A/B rewrite that removes a source's lateral step, for example
  by assigning sibling helpers a common rank-2 frame.
- Whether PR1937's opposite split (2,920 → 2,921 children) can be reversed.

Credit: the PR1122, PR1194, PR1311, PR1806, PR1937, PR2033 and PR2037 authors,
and the PR589/PR555 bit word.
Prepared by Philippe Van Velzen with substantial Claude (Anthropic) assistance.
A sub-agent ran the search, and the main session re-ran every script from a
clean main checkout. Apache-2.0.
