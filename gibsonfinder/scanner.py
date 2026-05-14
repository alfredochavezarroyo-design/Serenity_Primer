"""
scanner.py — Primer scanning engine for GibsonFinder.

Scans one or two flanking windows around a theoretical insert and
enumerates all Gibson Assembly primer candidates that satisfy the
user-specified Tm, GC%, and binding-length constraints.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional

from .thermo import calc_tm, gc_percent, score_primer, rev_comp, is_valid_dna
from .dimers import DimerResult, self_dimer


# ── Primer data class ─────────────────────────────────────────────────────────

@dataclass
class Primer:
    direction:    str           # 'fwd' or 'rev'
    region:       str           # 'upstream', 'downstream', or 'window'
    pos:          int           # 1-based position on the template
    bind_seq:     str           # binding region (5'→3', no overlap)
    overlap_seq:  str           # overlap tail (5'→3', prepended)
    full_primer:  str           # overlap + binding (5'→3', ordered to order)
    tm:           float         # Tm of binding region only
    gc:           float         # GC% of binding region
    score:        int           # 0–100 quality score
    bind_len:     int           # length of binding region
    total_len:    int           # total primer length
    self_dimer:   DimerResult = field(repr=False)

    @property
    def label(self) -> str:
        dir_tag = "FWD" if self.direction == "fwd" else "REV"
        return f"{dir_tag} pos={self.pos} [{self.region}]"

    def fasta_header(self) -> str:
        sd = self.self_dimer
        return (
            f">Primer_{self.direction.upper()}_pos{self.pos}"
            f"_Tm{self.tm}_GC{self.gc}_score{self.score}"
            f"_sd{sd.max_run}bp"
        )


# ── Scan engine ───────────────────────────────────────────────────────────────

@dataclass
class ScanParams:
    bind_min:    int   = 18
    bind_max:    int   = 25
    overlap_len: int   = 25
    tm_target:   float = 60.0
    tm_tol:      float = 5.0
    gc_min:      float = 35.0
    gc_max:      float = 65.0
    step:        int   = 3
    circular:    bool  = True
    max_results: Optional[int] = 100


def _scan_window(
    seq: str,
    raw: str,
    win_start: int,
    win_end: int,
    params: ScanParams,
    region_label: str,
) -> list[Primer]:
    """
    Scan [win_start, win_end) on `seq` and return all valid primers.

    `seq` is the working sequence (may be raw + suffix for circular).
    `raw` is the original sequence (used for circular overlap wrapping).
    """
    p = params
    primers: list[Primer] = []

    for pos in range(win_start, win_end, p.step):

        # ── Forward primer ────────────────────────────────────────────────
        for b_len in range(p.bind_min, p.bind_max + 1):
            if pos + b_len > len(seq):
                break
            bind_seq = seq[pos: pos + b_len]
            if not is_valid_dna(bind_seq):
                continue
            tm = calc_tm(bind_seq)
            gc = gc_percent(bind_seq)
            if abs(tm - p.tm_target) > p.tm_tol or not (p.gc_min <= gc <= p.gc_max):
                continue

            # Overlap: region immediately upstream of the binding site
            ov_start = pos - p.overlap_len
            if ov_start < 0:
                if p.circular:
                    # Wrap around the plasmid
                    ov_seq = raw[ov_start:] + raw[:pos]
                else:
                    break
            else:
                ov_seq = seq[ov_start: pos]

            if len(ov_seq) < p.overlap_len:
                break

            full = ov_seq + bind_seq
            sd   = self_dimer(full)
            primers.append(Primer(
                direction   = "fwd",
                region      = region_label,
                pos         = pos + 1,
                bind_seq    = bind_seq,
                overlap_seq = ov_seq,
                full_primer = full,
                tm          = tm,
                gc          = gc,
                score       = score_primer(bind_seq, p.tm_target),
                bind_len    = b_len,
                total_len   = len(full),
                self_dimer  = sd,
            ))
            break   # use the shortest binding length that satisfies constraints

        # ── Reverse primer ────────────────────────────────────────────────
        for b_len in range(p.bind_min, p.bind_max + 1):
            if pos + b_len > len(seq):
                break
            bind_top = seq[pos: pos + b_len]   # template strand
            if not is_valid_dna(bind_top):
                continue
            bind_seq = rev_comp(bind_top)       # primer reads 5'→3'
            tm = calc_tm(bind_seq)
            gc = gc_percent(bind_seq)
            if abs(tm - p.tm_target) > p.tm_tol or not (p.gc_min <= gc <= p.gc_max):
                continue

            # Overlap: region immediately downstream of binding site on template
            ov_top = seq[pos + b_len: pos + b_len + p.overlap_len]
            if len(ov_top) < p.overlap_len:
                break
            ov_seq = rev_comp(ov_top)          # overlap tail on the primer

            full = ov_seq + bind_seq
            sd   = self_dimer(full)
            primers.append(Primer(
                direction   = "rev",
                region      = region_label,
                pos         = pos + 1,
                bind_seq    = bind_seq,
                overlap_seq = ov_seq,
                full_primer = full,
                tm          = tm,
                gc          = gc,
                score       = score_primer(bind_seq, p.tm_target),
                bind_len    = b_len,
                total_len   = len(full),
                self_dimer  = sd,
            ))
            break


    return primers


def scan_single(
    raw: str,
    win_start: int,
    win_end: int,
    params: ScanParams,
) -> list[Primer]:
    """
    Single-window scan: forward AND reverse primers across [win_start, win_end).
    Positions are 1-based.
    """
    ws = win_start - 1   # convert to 0-based
    we = win_end
    seq = (raw + raw[: params.overlap_len + params.bind_max]) if params.circular else raw
    primers = _scan_window(seq, raw, ws, we, params, "window")
    primers.sort(key=lambda p: p.score, reverse=True)
    if params.max_results:
        primers = primers[: params.max_results]
    return primers


def scan_paired(
    raw: str,
    up_start: int, up_end: int,
    dn_start: int, dn_end: int,
    params: ScanParams,
) -> tuple[list[Primer], list[Primer]]:
    """
    Paired flanking scan.

    Upstream window  → forward primers only (binding into backbone, overlap into insert).
    Downstream window → reverse primers only.

    Returns (upstream_primers, downstream_primers) sorted by score descending.
    Positions are 1-based.
    """
    seq = (raw + raw[: params.overlap_len + params.bind_max]) if params.circular else raw

    all_up = _scan_window(seq, raw, up_start - 1, up_end, params, "upstream")
    all_dn = _scan_window(seq, raw, dn_start - 1, dn_end, params, "downstream")

    up_fwd = [p for p in all_up if p.direction == "fwd"]
    dn_rev = [p for p in all_dn if p.direction == "rev"]

    up_fwd.sort(key=lambda p: p.score, reverse=True)
    dn_rev.sort(key=lambda p: p.score, reverse=True)

    n = params.max_results or 9999
    return up_fwd[:n], dn_rev[:n]
