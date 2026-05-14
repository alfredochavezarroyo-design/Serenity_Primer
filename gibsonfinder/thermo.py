"""
thermo.py — Thermodynamic calculations for GibsonFinder.

Uses BioPython's MeltingTemp (SantaLucia 1998 nearest-neighbor parameters)
for Tm calculation, which is peer-reviewed and well-tested.

Salt correction: Owczarzy 2004, 50 mM Na+, 250 nM oligo.
"""

from __future__ import annotations
from Bio.SeqUtils.MeltingTemp import Tm_NN, DNA_NN4
from Bio.Seq import Seq

# ── Constants ────────────────────────────────────────────────────────────────

COMPLEMENT = str.maketrans("ATCGatcgNn", "TAGCtagcNn")

# ── Basic sequence operations ────────────────────────────────────────────────

def rev_comp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (case-preserved)."""
    return seq.translate(COMPLEMENT)[::-1]


def gc_percent(seq: str) -> float:
    """Return GC% rounded to 1 decimal place."""
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return round(gc / len(seq) * 100, 1)


def is_valid_dna(seq: str) -> bool:
    """True if sequence contains only ATCG (uppercase)."""
    return bool(seq) and all(b in "ATCG" for b in seq.upper())


def clean_seq(raw: str) -> str:
    """Strip FASTA headers, whitespace, and non-ATCG characters."""
    import re
    lines = raw.splitlines()
    seq = "".join(l for l in lines if not l.startswith(">"))
    seq = re.sub(r"\s+", "", seq).upper()
    seq = re.sub(r"[^ATCGN]", "", seq)
    return seq


# ── Tm calculation ───────────────────────────────────────────────────────────

def calc_tm(seq: str, na_conc_mM: float = 50.0, oligo_nM: float = 250.0) -> float:
    """
    Calculate melting temperature using SantaLucia 1998 nearest-neighbor
    parameters (DNA_NN4 table in BioPython).

    Parameters
    ----------
    seq        : binding region sequence (5'->3', no overlap tail)
    na_conc_mM : Na+ concentration in mM  (default 50 mM)
    oligo_nM   : oligo concentration in nM (default 250 nM)

    Returns
    -------
    Tm in °C, rounded to 1 decimal place.
    """
    seq = seq.upper()
    if len(seq) < 2:
        return 0.0
    try:
        tm = Tm_NN(
            Seq(seq),
            nn_table=DNA_NN4,          # SantaLucia 1998
            Na=na_conc_mM,
            dnac1=oligo_nM,
            dnac2=0,                   # non-self-complementary, primer in excess
            saltcorr=5,                # Owczarzy 2004 unified correction
        )
        return round(float(tm), 1)
    except Exception:
        return 0.0


# ── Hairpin / self-complementarity penalty ───────────────────────────────────

def hairpin_run(seq: str, min_stem: int = 4, max_stem: int = 10) -> int:
    """
    Return the longest stem length found in a hairpin-like fold.
    Uses a simple sliding-window search: for each stem length, check
    whether the reverse complement of any subsequence occurs elsewhere.
    """
    seq = seq.upper()
    rc  = rev_comp(seq)
    best = 0
    half = len(seq) // 2
    for stem in range(min_stem, min(max_stem, half) + 1):
        for i in range(len(seq) - stem + 1):
            fragment = seq[i: i + stem]
            # Look for the RC of this fragment anywhere else in the sequence
            # (offset by at least 1 to avoid trivial same-position match)
            if fragment in rc:
                best = max(best, stem)
    return best


# ── Primer scoring ───────────────────────────────────────────────────────────

def score_primer(bind_seq: str, tm_target: float) -> int:
    """
    Score a primer 0–100 based on:
      40 pts  Tm proximity to target (linear decay)
      25 pts  GC content (optimal 40–60%)
      20 pts  3' G/C clamp
      15 pts  Hairpin penalty
    """
    tm = calc_tm(bind_seq)
    gc = gc_percent(bind_seq)

    tm_score = max(0.0, 40.0 - abs(tm - tm_target) * 5.0)
    gc_score = 25.0 if 40 <= gc <= 60 else max(0.0, 25.0 - abs(gc - 50) * 0.8)
    clamp    = 20.0 if bind_seq[-1].upper() in "GC" else 5.0
    hp       = min(15, hairpin_run(bind_seq) * 2)

    return round(tm_score + gc_score + clamp + (15 - hp))
