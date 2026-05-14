"""
reporter.py — Rich-formatted terminal output for GibsonFinder.
"""

from __future__ import annotations
from typing import Optional

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text
from rich.columns import Columns
from rich import box
from rich.rule import Rule
from rich.syntax import Syntax
from rich.padding import Padding

from .scanner import Primer, ScanParams
from .dimers import DimerResult, dimer_analysis

console = Console()


# ── Colour helpers ────────────────────────────────────────────────────────────

def risk_color(level: str) -> str:
    return {"none": "green", "low": "yellow", "high": "red"}.get(level, "white")

def score_color(score: int) -> str:
    if score >= 70: return "bright_green"
    if score >= 45: return "yellow"
    return "red"

def tm_color(tm: float, target: float, tol: float) -> str:
    if abs(tm - target) <= 2:  return "bright_green"
    if abs(tm - target) <= tol: return "yellow"
    return "red"


# ── Alignment block ───────────────────────────────────────────────────────────

def render_alignment(d: DimerResult, label_a: str = "A  5′→3′", label_b: str = "B  3′→5′") -> Text:
    """
    Build a Rich Text object for a dimer alignment.
    Displays:
        Seq A (5'→3'):  ATCG...
                        ||||
        Seq B (3'→5'):  TAGC...
    with | marks coloured red and matched bases highlighted.
    """
    t = Text()
    pad = max(len(label_a), len(label_b)) + 2

    t.append(label_a.ljust(pad), style="dim")
    t.append(d.line_a + "\n", style="cyan")

    t.append(" " * pad)
    for ch in d.line_m:
        if ch == "|":
            t.append("|", style="bold red")
        else:
            t.append(" ")
    t.append("\n")

    t.append(label_b.ljust(pad), style="dim")
    t.append(d.line_b + "\n", style="yellow")

    risk_col = risk_color(d.risk_level)
    t.append(f"\n  Max run: ", style="dim")
    t.append(f"{d.max_run} bp", style=f"bold {risk_col}")
    if d.three_prime_risk:
        t.append("  ⚠ 3′ end involved", style="bold red")
    t.append(f"  [{d.risk_level.upper()}]\n", style=risk_col)

    return t


# ── Primer sequence coloring ──────────────────────────────────────────────────

def colored_primer(p: Primer) -> Text:
    """Return primer sequence with overlap in cyan, binding in green."""
    t = Text()
    t.append(p.overlap_seq, style="cyan")
    t.append(p.bind_seq, style="bright_green")
    return t


# ── Single-window table ───────────────────────────────────────────────────────

def print_primers_table(
    primers: list[Primer],
    params: ScanParams,
    title: str = "Gibson Primer Candidates",
) -> None:
    table = Table(
        title       = title,
        box         = box.MINIMAL_HEAVY_HEAD,
        border_style= "bright_black",
        header_style= "bold cyan",
        show_lines  = True,
    )

    table.add_column("#",         style="dim",          width=4)
    table.add_column("Dir",       width=5)
    table.add_column("Pos",       style="bright_black", width=6)
    table.add_column("Primer 5′→3′ (overlap | binding)", no_wrap=False, min_width=42)
    table.add_column("Tm",        width=7)
    table.add_column("GC%",       width=6)
    table.add_column("Len",       style="dim",          width=7)
    table.add_column("Self-dimer",width=14)
    table.add_column("Score",     width=7)

    for i, p in enumerate(primers, 1):
        sd      = p.self_dimer
        sd_col  = risk_color(sd.risk_level)
        sd_text = Text(sd.summary(), style=sd_col)
        tm_col  = tm_color(p.tm, params.tm_target, params.tm_tol)

        dir_style = "green" if p.direction == "fwd" else "blue"
        dir_tag   = Text(p.direction.upper(), style=f"bold {dir_style}")

        table.add_row(
            str(i),
            dir_tag,
            str(p.pos),
            colored_primer(p),
            Text(f"{p.tm}°C", style=tm_col),
            f"{p.gc}%",
            f"{p.total_len}/{p.bind_len}",
            sd_text,
            Text(str(p.score), style=score_color(p.score)),
        )

    console.print()
    console.print(table)


# ── Paired analysis ───────────────────────────────────────────────────────────

def print_paired_analysis(
    fwds: list[Primer],
    revs: list[Primer],
    params: ScanParams,
    top_n: int = 5,
) -> None:
    console.print()
    console.print(Rule("[bold cyan]Paired Analysis — FWD × REV combinations[/]"))

    fwd_top = fwds[:top_n]
    rev_top = revs[:top_n]

    if not fwd_top or not rev_top:
        console.print("[yellow]Not enough primers in one or both windows.[/]")
        return

    pair_num = 0
    for fwd in fwd_top:
        for rev in rev_top:
            pair_num += 1
            hd        = dimer_analysis(fwd.full_primer, rev.full_primer)
            hd_col    = risk_color(hd.risk_level)
            tm_diff   = abs(fwd.tm - rev.tm)
            tm_ok     = tm_diff <= 5
            pair_score= round((fwd.score + rev.score) / 2)

            title_str = (
                f"Pair #{pair_num}  ·  "
                f"FWD pos {fwd.pos}  ×  REV pos {rev.pos}  ·  "
                f"ΔTm {tm_diff:.1f}°C  ·  PairScore {pair_score}"
            )

            # ── Header panel ──────────────────────────────────────────
            console.print()
            console.print(Panel(
                title_str,
                style   = "bold white on bright_black",
                box     = box.HEAVY,
                expand  = True,
            ))

            # ── Two-column primer detail ──────────────────────────────
            def primer_block(p: Primer, label: str) -> str:
                sd     = p.self_dimer
                sd_col = risk_color(sd.risk_level)
                color  = "green" if p.direction == "fwd" else "blue"
                lines  = [
                    f"[bold {color}]{label} (5′→3′)[/]",
                    f"[cyan]{p.overlap_seq}[/][bright_green]{p.bind_seq}[/]",
                    f"Tm [bold]{p.tm}°C[/]  GC [bold]{p.gc}%[/]  Score [bold]{p.score}[/]",
                    f"Self-dimer: [{sd_col}]{sd.summary()}[/]",
                    f"Total {p.total_len} bp  ({p.bind_len} bp binding)",
                ]
                return "\n".join(lines)

            console.print(Columns([
                Panel(primer_block(fwd, "Forward Primer"), border_style="green",  expand=True),
                Panel(primer_block(rev, "Reverse Primer"), border_style="blue", expand=True),
            ]))

            # ── Self-dimer alignments (if flagged) ────────────────────
            for p, lbl in [(fwd, "FWD self-dimer"), (rev, "REV self-dimer")]:
                if not p.self_dimer.is_clean:
                    console.print(
                        Panel(
                            render_alignment(p.self_dimer,
                                             f"{lbl} 5′→3′",
                                             f"{lbl} 3′→5′"),
                            title       = f"[yellow]{lbl}[/]",
                            border_style= "yellow",
                        )
                    )

            # ── Heterodimer alignment ─────────────────────────────────
            if hd.risk_level != "none":
                hd_title = f"[{hd_col}]Heterodimer — {hd.max_run} bp run[/]"
                console.print(
                    Panel(
                        render_alignment(hd,
                                         "FWD  5′→3′",
                                         "REV  3′→5′"),
                        title       = hd_title,
                        border_style= hd_col,
                    )
                )
            else:
                console.print(
                    Panel("[green]✓ No significant heterodimer detected[/]",
                          border_style="green")
                )


# ── Dimer check table ─────────────────────────────────────────────────────────

def print_dimer_table(
    primers: list[Primer],
    show_alignments: bool = False,
) -> None:
    top = primers[:20]
    console.print()
    console.print(Rule("[bold cyan]Dimer Analysis[/]"))

    table = Table(
        box         = box.MINIMAL_HEAVY_HEAD,
        border_style= "bright_black",
        header_style= "bold cyan",
        show_lines  = False,
    )
    table.add_column("Type",    width=8)
    table.add_column("Primer A",width=22)
    table.add_column("Primer B",width=22)
    table.add_column("Max run", width=9)
    table.add_column("3′ end",  width=7)
    table.add_column("Risk",    width=10)

    rows = []

    # Self-dimer rows
    for p in primers[:30]:
        sd = p.self_dimer
        rows.append(("Self", p, p, sd))

    # Heterodimer rows (top 20 × top 20)
    for i, pa in enumerate(top):
        for pb in top[i+1:]:
            hd = dimer_analysis(pa.full_primer, pb.full_primer)
            rows.append(("Hetero", pa, pb, hd))

    # Sort by max_run descending
    rows.sort(key=lambda r: r[3].max_run, reverse=True)

    for kind, pa, pb, d in rows:
        rc = risk_color(d.risk_level)
        b_name = "(self)" if kind == "Self" else pb.label
        table.add_row(
            Text(kind, style="dim"),
            pa.label,
            b_name,
            Text(f"{d.max_run} bp", style=rc),
            Text("⚠ Yes" if d.three_prime_risk else "No",
                 style="red" if d.three_prime_risk else "dim"),
            Text(d.risk_level.upper(), style=f"bold {rc}"),
        )

        if show_alignments and d.max_run > 0:
            la = pa.label + " 5′→3′"
            lb = ("(self) 3′→5′" if kind == "Self"
                  else pb.label + " 3′→5′")
            console.print(table)
            console.print(Padding(render_alignment(d, la, lb), (0, 4)))
            table = Table(box=box.MINIMAL_HEAVY_HEAD, border_style="bright_black",
                          header_style="bold cyan", show_lines=False)
            table.add_column("Type",    width=8)
            table.add_column("Primer A",width=22)
            table.add_column("Primer B",width=22)
            table.add_column("Max run", width=9)
            table.add_column("3′ end",  width=7)
            table.add_column("Risk",    width=10)

    console.print(table)


# ── Summary banner ────────────────────────────────────────────────────────────

def print_summary(
    primers: list[Primer],
    params: ScanParams,
    seq_len: int,
    mode: str,
) -> None:
    clean  = sum(1 for p in primers if p.self_dimer.is_clean)
    fwd_n  = sum(1 for p in primers if p.direction == "fwd")
    rev_n  = sum(1 for p in primers if p.direction == "rev")
    top_tm = [p.tm for p in primers[:10]]
    tm_range = f"{min(top_tm):.1f}–{max(top_tm):.1f}°C" if top_tm else "—"

    lines = [
        f"[bold]Mode:[/]          {mode}",
        f"[bold]Sequence:[/]      {seq_len:,} bp",
        f"[bold]Candidates:[/]    {len(primers)} found",
        f"[bold]  Forward:[/]     {fwd_n}",
        f"[bold]  Reverse:[/]     {rev_n}",
        f"[bold]  Clean SD:[/]    {clean} / {len(primers)}",
        f"[bold]Top-10 Tm range:[/] {tm_range}",
        f"[bold]Target Tm:[/]     {params.tm_target}°C ± {params.tm_tol}°C",
        f"[bold]Overlap:[/]       {params.overlap_len} bp",
        f"[bold]Binding:[/]       {params.bind_min}–{params.bind_max} bp",
    ]
    console.print()
    console.print(Panel(
        "\n".join(lines),
        title       = "[bold bright_green]GibsonFinder — Results Summary[/]",
        border_style= "bright_green",
    ))
