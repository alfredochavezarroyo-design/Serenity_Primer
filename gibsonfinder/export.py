"""
export.py — CSV, FASTA, and text export for GibsonFinder.
"""

from __future__ import annotations
import csv
import io
from pathlib import Path
from typing import Optional

from .scanner import Primer
from .dimers import dimer_analysis


# ── FASTA ─────────────────────────────────────────────────────────────────────

def to_fasta(primers: list[Primer], path: Optional[Path] = None) -> str:
    lines = []
    for i, p in enumerate(primers, 1):
        sd = p.self_dimer
        header = (
            f">Primer_{i}_{p.direction.upper()}_pos{p.pos}"
            f"_Tm{p.tm}_GC{p.gc}_score{p.score}"
            f"_sd{sd.max_run}bp"
        )
        lines.append(header)
        lines.append(p.full_primer)

    fasta = "\n".join(lines) + "\n"
    if path:
        path.write_text(fasta)
    return fasta


# ── CSV: all primers ──────────────────────────────────────────────────────────

PRIMER_FIELDS = [
    "#", "Direction", "Region", "Position",
    "Full_Primer_5to3", "Overlap_Seq", "Binding_Seq",
    "Tm_C", "GC_pct", "Bind_len_bp", "Total_len_bp",
    "Score", "SelfDimer_max_run_bp", "SelfDimer_3prime_risk", "SelfDimer_risk_level",
]

def to_csv(primers: list[Primer], path: Optional[Path] = None) -> str:
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=PRIMER_FIELDS, lineterminator="\n")
    writer.writeheader()
    for i, p in enumerate(primers, 1):
        sd = p.self_dimer
        writer.writerow({
            "#":                      i,
            "Direction":              p.direction,
            "Region":                 p.region,
            "Position":               p.pos,
            "Full_Primer_5to3":       p.full_primer,
            "Overlap_Seq":            p.overlap_seq,
            "Binding_Seq":            p.bind_seq,
            "Tm_C":                   p.tm,
            "GC_pct":                 p.gc,
            "Bind_len_bp":            p.bind_len,
            "Total_len_bp":           p.total_len,
            "Score":                  p.score,
            "SelfDimer_max_run_bp":   sd.max_run,
            "SelfDimer_3prime_risk":  "yes" if sd.three_prime_risk else "no",
            "SelfDimer_risk_level":   sd.risk_level,
        })

    result = buf.getvalue()
    if path:
        path.write_text(result)
    return result


# ── CSV: primer pairs ─────────────────────────────────────────────────────────

PAIR_FIELDS = [
    "Pair", "FWD_pos", "FWD_full", "FWD_Tm", "FWD_GC", "FWD_score",
    "FWD_sd_bp", "FWD_sd_risk",
    "REV_pos", "REV_full", "REV_Tm", "REV_GC", "REV_score",
    "REV_sd_bp", "REV_sd_risk",
    "Hetero_max_run_bp", "Hetero_3prime", "Hetero_risk",
    "DeltaTm_C", "Pair_score",
]

def pairs_to_csv(
    fwds: list[Primer],
    revs: list[Primer],
    top_n: int = 5,
    path: Optional[Path] = None,
) -> str:
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=PAIR_FIELDS, lineterminator="\n")
    writer.writeheader()

    n = 0
    for fwd in fwds[:top_n]:
        for rev in revs[:top_n]:
            n += 1
            hd   = dimer_analysis(fwd.full_primer, rev.full_primer)
            delt = abs(fwd.tm - rev.tm)
            writer.writerow({
                "Pair":              n,
                "FWD_pos":           fwd.pos,
                "FWD_full":          fwd.full_primer,
                "FWD_Tm":            fwd.tm,
                "FWD_GC":            fwd.gc,
                "FWD_score":         fwd.score,
                "FWD_sd_bp":         fwd.self_dimer.max_run,
                "FWD_sd_risk":       fwd.self_dimer.risk_level,
                "REV_pos":           rev.pos,
                "REV_full":          rev.full_primer,
                "REV_Tm":            rev.tm,
                "REV_GC":            rev.gc,
                "REV_score":         rev.score,
                "REV_sd_bp":         rev.self_dimer.max_run,
                "REV_sd_risk":       rev.self_dimer.risk_level,
                "Hetero_max_run_bp": hd.max_run,
                "Hetero_3prime":     "yes" if hd.three_prime_risk else "no",
                "Hetero_risk":       hd.risk_level,
                "DeltaTm_C":         round(delt, 2),
                "Pair_score":        round((fwd.score + rev.score) / 2),
            })

    result = buf.getvalue()
    if path:
        path.write_text(result)
    return result
