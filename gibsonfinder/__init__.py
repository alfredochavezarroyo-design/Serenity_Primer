"""
GibsonFinder — Gibson Assembly primer design toolkit.
"""

from .thermo import calc_tm, gc_percent, rev_comp, clean_seq, score_primer
from .dimers import dimer_analysis, self_dimer, DimerResult
from .scanner import Primer, ScanParams, scan_single, scan_paired
from .export import to_fasta, to_csv, pairs_to_csv
from .reporter import (
    console,
    print_primers_table,
    print_paired_analysis,
    print_dimer_table,
    print_summary,
)

__version__ = "2.1.0"
__all__ = [
    "calc_tm", "gc_percent", "rev_comp", "clean_seq", "score_primer",
    "dimer_analysis", "self_dimer", "DimerResult",
    "Primer", "ScanParams", "scan_single", "scan_paired",
    "to_fasta", "to_csv", "pairs_to_csv",
    "console", "print_primers_table", "print_paired_analysis",
    "print_dimer_table", "print_summary",
]
