# Followers 670/671 on PR2802's words: a lender pigeonhole (FINDING, no κ)

`movecen.py` runs a hostee-move / augmenting-path census on a BORROW-WORD. Every hosted life is
treated as movable: all 72 current hostees (B-warp anticipatory borrows, pre-life hostings and
retired-lender aliases), plus the two free followers 670 and 671. The census uses exact lattice frame
bounds over every refit of the fixed gate list (`bounds.py`). An edge from a life to a register Y
exists when:
- Y is not restored;
- Y's own life ends before the hosted life starts;
- dim lo(Y's last gate) ≤ dim hi(the life's first gate).

**Result on PR2802** (`research/t74-joint-twelve-retired7/inputs/WORD-{base,extra}.json.gz`; base and
extra give the same result):
- every life's first-gate ceiling is ≤ 5;
- exactly 72 registers have floor ≤ 5 (62 at floor 4, 10 at floor 5), and all of them already host a life;
- every free register has floor ≥ 6 (124 at 6, 14 at 7, 192 at 8);
- the 74 lives overlap pairwise in time (they start at B11292–B12281 and end at B12548–B12976), so no
  register can host two of them one after the other.

So 74 lives cannot fit on 72 registers: followers 670/671 have no lender under any refit of the fixed
gate list (0 augmenting paths). Opening a slot would need a gate-level (kernel) edit. Pre-life and
refit censuses on the same words also find 0 new legal hostings or aliases. PR2802's fill is at the
class-level ceiling (to about 1e-14).

Run: `python3 movecen.py <WORD.json.gz> <this dir>`. It uses `bounds.py` and `wordmodel.py`, both in
this directory.

Prepared with substantial Claude Code assistance (a sub-agent did the census; the main session re-ran it).
