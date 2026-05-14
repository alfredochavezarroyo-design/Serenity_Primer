"""
dimers.py — Antiparallel Watson-Crick dimer analysis for GibsonFinder.

Model
-----
seqA runs 5'→3': A[0] A[1] ... A[n-1]
seqB runs 5'→3': B[0] B[1] ... B[m-1]

For antiparallel alignment, B is read 3'→5'. Define:
    rcB = revComp(seqB)          rcB[k] = complement( B[m-1-k] )

At shift s, rcB[0] aligns with A[s]:
    A[i] pairs with rcB[i-s]    for max(0,s) <= i < min(n, s+m)

MATCH CONDITION: A[i] == rcB[i-s]
    Proof: need complement of B read 3'->5' at that column.
    B read 3'->5' at column i is B[m-1-(i-s)].
    rcB[i-s] = complement( B[m-1-(i-s)] ).
    So A[i] == rcB[i-s]  <=>  A[i] == complement(B[m-1-(i-s)])  ✓ WC.

DISPLAY (conventional antiparallel):
    Top:    A[i]                      (seqA 5'→3')
    Middle: '|' where A[i]==rcB[i-s]  (Watson-Crick pair)
    Bottom: B[m-1-(i-s)]              (seqB 3'→5', actual base)

Note: bottom = complement(rcB[i-s]), NOT rcB[i-s] itself.
Showing rcB on the bottom was the critical bug in the JS version —
it caused A-A / T-T display for matched pairs instead of A-T / G-C.

SELF-DIMER:
    seqA == seqB. No shift is excluded. s=0 means A[i] pairs with
    A[n-1-i], a valid intramolecular fold (hairpin geometry).
"""

from __future__ import annotations
from dataclasses import dataclass, field
from .thermo import rev_comp


# ── Data classes ─────────────────────────────────────────────────────────────

@dataclass
class DimerResult:
    max_run: int          # length of longest contiguous complementary run
    three_prime_risk: bool  # run of >=4 bp overlaps 3' end of either primer
    line_a: str           # display: seqA (5'→3')
    line_m: str           # display: match marks '|'
    line_b: str           # display: seqB (3'→5', actual bases)
    best_shift: int       # shift that produced the longest run

    @property
    def risk_level(self) -> str:
        if self.max_run >= 8 or (self.max_run >= 4 and self.three_prime_risk):
            return "high"
        if self.max_run >= 6 or self.three_prime_risk:
            return "low"
        return "none"

    @property
    def is_clean(self) -> bool:
        return self.risk_level == "none"

    def summary(self) -> str:
        if self.max_run == 0:
            return "Clean"
        tp = " 3′!" if self.three_prime_risk else ""
        return f"{self.max_run}bp{tp} [{self.risk_level}]"


# ── Core sliding-window alignment ────────────────────────────────────────────

def dimer_analysis(seq_a: str, seq_b: str, max_display_cols: int = 80) -> DimerResult:
    """
    Full antiparallel dimer analysis between seq_a and seq_b.

    Parameters
    ----------
    seq_a, seq_b     : primer sequences 5'→3' (case-insensitive)
    max_display_cols : truncate alignment display to this many columns

    Returns
    -------
    DimerResult with run statistics and alignment display lines.
    """
    a   = seq_a.upper()
    b   = seq_b.upper()
    rc_b = rev_comp(b)      # rcB[k] = complement(b[m-1-k])
    n, m = len(a), len(b)

    # ── Find longest contiguous WC run at any shift ───────────────────────
    max_run      = 0
    best_shift   = 0
    best_run_start = 0   # position in `a` where the best run begins

    for s in range(-(m - 1), n):
        run = 0
        run_start = 0
        i_start = max(0, s)
        i_end   = min(n, s + m)

        for i in range(i_start, i_end):
            rc_idx = i - s          # index into rc_b
            if a[i] == rc_b[rc_idx]:
                if run == 0:
                    run_start = i
                run += 1
                if run > max_run:
                    max_run       = run
                    best_shift    = s
                    best_run_start = run_start
            else:
                run = 0

    # ── 3' risk: any run >=4 overlapping last 5 bases of A or B ──────────
    tail_3a = n - 5      # A positions >= tail_3a are "3' end of A"
    tail_3b = m - 5      # rcB positions >= tail_3b correspond to "3' end of B"
    three_prime_risk = False

    outer_done = False
    for s in range(-(m - 1), n):
        if outer_done:
            break
        run = 0
        run_start = 0
        i_start = max(0, s)
        i_end   = min(n, s + m)

        for i in range(i_start, i_end):
            rc_idx = i - s
            if a[i] == rc_b[rc_idx]:
                if run == 0:
                    run_start = i
                run += 1
                if run >= 4:
                    rc_run_end = i - s
                    # 3' of A: current position i is in last 5
                    if i >= tail_3a:
                        three_prime_risk = True
                        outer_done = True
                        break
                    # 3' of B: rc_b positions >= tail_3b → 3' of B
                    if rc_run_end >= tail_3b:
                        three_prime_risk = True
                        outer_done = True
                        break
            else:
                run = 0

    # ── Build display alignment at best_shift ────────────────────────────
    # Window spans all positions with at least one strand present.
    disp_start = min(0, best_shift)
    disp_end   = max(n, best_shift + m)

    # Centre around best run if window is wider than max_display_cols
    if disp_end - disp_start > max_display_cols:
        centre   = best_run_start + max_run // 2
        disp_start = centre - max_display_cols // 2
        disp_end   = disp_start + max_display_cols

    line_a, line_m, line_b = [], [], []

    for pos in range(disp_start, disp_end):
        ai  = pos
        rci = pos - best_shift          # index into rc_b AND into b (via m-1-rci)

        # seqA character at this column
        ac = a[ai]   if 0 <= ai  < n else " "
        # rcB character (for match test)
        rc = rc_b[rci] if 0 <= rci < m else " "
        # seqB character read 3'→5' = b[m-1-rci]   ← THE CORRECTED LINE
        bc = b[m - 1 - rci] if 0 <= rci < m else " "

        is_match = (ac != " " and rc != " " and ac == rc)

        line_a.append(ac)
        line_m.append("|" if is_match else " ")
        line_b.append(bc)

    return DimerResult(
        max_run        = max_run,
        three_prime_risk = three_prime_risk,
        line_a         = "".join(line_a),
        line_m         = "".join(line_m),
        line_b         = "".join(line_b),
        best_shift     = best_shift,
    )


def self_dimer(seq: str) -> DimerResult:
    """Convenience wrapper: dimer of a primer with itself."""
    return dimer_analysis(seq, seq)
