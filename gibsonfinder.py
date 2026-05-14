#!/usr/bin/env python3
"""
gibsonfinder — CLI for Gibson Assembly primer design.

Usage examples
--------------
# Single-window scan (prints table + dimer check):
  python gibsonfinder.py single sequence.fasta --win-start 1 --win-end 200

# Paired flanking scan (FWD upstream × REV downstream, paired analysis):
  python gibsonfinder.py paired sequence.fasta \\
      --up-start 1 --up-end 150 \\
      --dn-start 501 --dn-end 650 \\
      --insert-start 151 --insert-end 500

# Tune parameters:
  python gibsonfinder.py single seq.fa \\
      --tm-target 62 --tm-tol 4 \\
      --overlap 25 --bind-min 18 --bind-max 25 \\
      --gc-min 40 --gc-max 65 --step 3

# Export results:
  python gibsonfinder.py single seq.fa --out-csv primers.csv --out-fasta primers.fasta

# Show dimer alignments inline:
  python gibsonfinder.py single seq.fa --show-alignments
"""

import sys
from pathlib import Path

import click
from rich.console import Console
from rich.rule import Rule

# Allow running from the parent directory of the gibsonfinder package
sys.path.insert(0, str(Path(__file__).parent))

from gibsonfinder import (
    clean_seq,
    ScanParams,
    scan_single,
    scan_paired,
    print_primers_table,
    print_paired_analysis,
    print_dimer_table,
    print_summary,
    to_csv,
    to_fasta,
    pairs_to_csv,
    console,
)

# ── Shared parameter decorators ───────────────────────────────────────────────

_COMMON = [
    click.option("--bind-min",   default=18,   show_default=True, help="Min binding length (bp)"),
    click.option("--bind-max",   default=25,   show_default=True, help="Max binding length (bp)"),
    click.option("--overlap",    default=25,   show_default=True, help="Overlap tail length (bp)"),
    click.option("--tm-target",  default=60.0, show_default=True, help="Target Tm (°C)"),
    click.option("--tm-tol",     default=5.0,  show_default=True, help="Tm tolerance ± (°C)"),
    click.option("--gc-min",     default=35.0, show_default=True, help="Min GC% for binding region"),
    click.option("--gc-max",     default=65.0, show_default=True, help="Max GC% for binding region"),
    click.option("--step",       default=3,    show_default=True, help="Scan step size (bp)"),
    click.option("--circular/--linear", default=True, show_default=True, help="Topology"),
    click.option("--max-results",default=100,  show_default=True, help="Max primers to report (0=all)"),
    click.option("--out-csv",    default=None, help="Export primers to CSV file"),
    click.option("--out-fasta",  default=None, help="Export primers to FASTA file"),
    click.option("--show-alignments", is_flag=True, default=False, help="Print dimer alignments in dimer table"),
    click.option("--top-pairs",  default=5,    show_default=True, help="Top N from each window for paired analysis"),
]

def add_common(func):
    for decorator in reversed(_COMMON):
        func = decorator(func)
    return func


def _make_params(**kw) -> ScanParams:
    return ScanParams(
        bind_min    = kw["bind_min"],
        bind_max    = kw["bind_max"],
        overlap_len = kw["overlap"],
        tm_target   = kw["tm_target"],
        tm_tol      = kw["tm_tol"],
        gc_min      = kw["gc_min"],
        gc_max      = kw["gc_max"],
        step        = kw["step"],
        circular    = kw["circular"],
        max_results = kw["max_results"] or None,
    )


def _load_seq(seq_file: str) -> str:
    path = Path(seq_file)
    if not path.exists():
        console.print(f"[red]File not found: {seq_file}[/]")
        sys.exit(1)
    raw = path.read_text(errors="ignore")
    seq = clean_seq(raw)
    if len(seq) < 30:
        console.print("[red]Sequence too short (< 30 bp) or no valid ATCG bases found.[/]")
        sys.exit(1)
    return seq


# ── CLI root ──────────────────────────────────────────────────────────────────

@click.group()
def cli():
    """
    GibsonFinder v2.1 — Gibson Assembly Primer Designer

    Design primers with correctly computed Tm (SantaLucia 1998 via BioPython),
    validated antiparallel dimer detection, and rich terminal output.
    """
    console.print()
    console.rule("[bold bright_green]GibsonFinder v2.1[/]")


# ── SINGLE command ────────────────────────────────────────────────────────────

@cli.command()
@click.argument("seq_file")
@click.option("--win-start", default=1,   show_default=True, help="Window start (1-based bp)")
@click.option("--win-end",   default=500, show_default=True, help="Window end (1-based bp)")
@add_common
def single(seq_file, win_start, win_end, **kw):
    """
    Single-window scan: find FWD and REV Gibson primers within a window.

    SEQ_FILE can be a FASTA file or plain text with just the sequence.
    """
    seq    = _load_seq(seq_file)
    params = _make_params(**kw)

    # Clamp window to sequence length
    win_end = min(win_end, len(seq))

    console.print(f"[dim]Sequence: {len(seq):,} bp  ·  Window: {win_start}–{win_end}  ·  Topology: {'circular' if params.circular else 'linear'}[/]")
    console.print()

    with console.status("[cyan]Scanning window...[/]"):
        primers = scan_single(seq, win_start, win_end, params)

    if not primers:
        console.print("[yellow]No primers found. Try relaxing --tm-tol, --gc-min/max, or widening the window.[/]")
        return

    print_summary(primers, params, len(seq), mode="single")
    print_primers_table(primers, params, title=f"Single Window Scan — {win_start}:{win_end}")
    print_dimer_table(primers, show_alignments=kw["show_alignments"])

    # Exports
    if kw["out_csv"]:
        p = Path(kw["out_csv"])
        to_csv(primers, p)
        console.print(f"[green]CSV saved → {p}[/]")
    if kw["out_fasta"]:
        p = Path(kw["out_fasta"])
        to_fasta(primers, p)
        console.print(f"[green]FASTA saved → {p}[/]")


# ── PAIRED command ────────────────────────────────────────────────────────────

@cli.command()
@click.argument("seq_file")
@click.option("--up-start",      required=True, type=int, help="Upstream window start (bp)")
@click.option("--up-end",        required=True, type=int, help="Upstream window end (bp)")
@click.option("--dn-start",      required=True, type=int, help="Downstream window start (bp)")
@click.option("--dn-end",        required=True, type=int, help="Downstream window end (bp)")
@click.option("--insert-start",  default=None,  type=int, help="Insert region start (display only)")
@click.option("--insert-end",    default=None,  type=int, help="Insert region end (display only)")
@click.option("--out-pairs-csv", default=None,            help="Export primer pairs to CSV")
@add_common
def paired(seq_file, up_start, up_end, dn_start, dn_end,
           insert_start, insert_end, out_pairs_csv, **kw):
    """
    Paired flanking scan: FWD primers from upstream window, REV from downstream.
    All FWD × REV combinations reported with heterodimer analysis.

    SEQ_FILE can be a FASTA file or plain text with just the sequence.
    """
    seq    = _load_seq(seq_file)
    params = _make_params(**kw)

    ins_str = ""
    if insert_start and insert_end:
        ins_str = f"  ·  Insert: {insert_start}–{insert_end}"

    console.print(
        f"[dim]Sequence: {len(seq):,} bp  ·  "
        f"Upstream: {up_start}–{up_end}  ·  "
        f"Downstream: {dn_start}–{dn_end}{ins_str}  ·  "
        f"Topology: {'circular' if params.circular else 'linear'}[/]"
    )

    with console.status("[cyan]Scanning flanking windows...[/]"):
        fwds, revs = scan_paired(
            seq, up_start, up_end, dn_start, dn_end, params
        )

    all_primers = fwds + revs
    if not all_primers:
        console.print("[yellow]No primers found. Try relaxing parameters or widening windows.[/]")
        return

    print_summary(all_primers, params, len(seq), mode="paired")

    console.print()
    console.print(Rule("[bold cyan]Upstream FWD Primers[/]"))
    print_primers_table(fwds, params, title=f"Upstream FWD — {up_start}:{up_end}")

    console.print()
    console.print(Rule("[bold cyan]Downstream REV Primers[/]"))
    print_primers_table(revs, params, title=f"Downstream REV — {dn_start}:{dn_end}")

    print_paired_analysis(fwds, revs, params, top_n=kw["top_pairs"])
    print_dimer_table(all_primers, show_alignments=kw["show_alignments"])

    # Exports
    if kw["out_csv"]:
        p = Path(kw["out_csv"])
        to_csv(all_primers, p)
        console.print(f"[green]All primers CSV saved → {p}[/]")
    if kw["out_fasta"]:
        p = Path(kw["out_fasta"])
        to_fasta(all_primers, p)
        console.print(f"[green]FASTA saved → {p}[/]")
    if out_pairs_csv:
        p = Path(out_pairs_csv)
        pairs_to_csv(fwds, revs, top_n=kw["top_pairs"], path=p)
        console.print(f"[green]Pairs CSV saved → {p}[/]")


# ── DIMER command — standalone dimer check for user-supplied sequences ─────────

@cli.command()
@click.argument("seq_a")
@click.argument("seq_b", required=False, default=None)
def dimer(seq_a, seq_b):
    """
    Quick dimer check for one or two sequences.

    SEQ_A: first primer sequence (5'→3')
    SEQ_B: second primer sequence (optional; omit for self-dimer only)

    Examples:
      python gibsonfinder.py dimer ATCGATCGATCG
      python gibsonfinder.py dimer ATCGATCGATCG CGATCGATCGAT
    """
    from gibsonfinder import dimer_analysis, self_dimer
    from gibsonfinder.reporter import render_alignment

    seq_a = seq_a.upper().strip()
    console.print(f"\n[bold]Seq A[/]: [cyan]{seq_a}[/]  ({len(seq_a)} bp)")

    # Self-dimer
    sd = self_dimer(seq_a)
    console.print()
    console.print(Rule("[yellow]Self-dimer[/]"))
    console.print(render_alignment(sd, "5′→3′", "3′→5′"))

    # Heterodimer
    if seq_b:
        seq_b = seq_b.upper().strip()
        console.print(f"[bold]Seq B[/]: [cyan]{seq_b}[/]  ({len(seq_b)} bp)")
        hd = dimer_analysis(seq_a, seq_b)
        console.print()
        console.print(Rule("[cyan]Heterodimer[/]"))
        console.print(render_alignment(hd, "A  5′→3′", "B  3′→5′"))


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    cli()
