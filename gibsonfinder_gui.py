"""
gibsonfinder_gui.py — Tkinter GUI for GibsonFinder.

Drop this file next to your gibsonfinder/ package folder and gibsonfinder.py.
Run with:  python gibsonfinder_gui.py

Requires: biopython, rich, click  (same as the CLI)
"""

import sys
import os
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from pathlib import Path

# ── Make sure the package is importable ──────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))

try:
    from gibsonfinder import (
        clean_seq, ScanParams, scan_single, scan_paired,
        to_csv, to_fasta, pairs_to_csv,
        dimer_analysis,
    )
    from gibsonfinder.reporter import render_alignment
    BACKEND_OK = True
except ImportError as e:
    BACKEND_OK = False
    IMPORT_ERROR = str(e)


# ─────────────────────────────────────────────────────────────────────────────
#  Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _int(var, default):
    try:
        return int(var.get())
    except Exception:
        return default

def _float(var, default):
    try:
        return float(var.get())
    except Exception:
        return default

def _bool(var):
    return var.get() == 1


# ─────────────────────────────────────────────────────────────────────────────
#  Main Application
# ─────────────────────────────────────────────────────────────────────────────

class GibsonFinderApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("GibsonFinder v2.1 — Gibson Assembly Primer Designer")
        self.resizable(True, True)
        self.configure(bg="#f5f5f0")

        # ── Fonts ────────────────────────────────────────────────────────────
        self.FONT_HEAD  = ("Courier", 10, "bold")
        self.FONT_LABEL = ("Courier", 9)
        self.FONT_ENTRY = ("Courier", 9)
        self.FONT_MONO  = ("Courier", 9)
        self.FONT_BTN   = ("Courier", 9, "bold")

        # ── Style ────────────────────────────────────────────────────────────
        style = ttk.Style(self)
        style.theme_use("clam")
        style.configure("TNotebook",            background="#f5f5f0")
        style.configure("TNotebook.Tab",        font=self.FONT_LABEL, padding=[10, 4])
        style.configure("TFrame",               background="#f5f5f0")
        style.configure("TLabel",               background="#f5f5f0", font=self.FONT_LABEL)
        style.configure("TLabelframe",          background="#f5f5f0")
        style.configure("TLabelframe.Label",    font=self.FONT_HEAD)
        style.configure("TCheckbutton",         background="#f5f5f0", font=self.FONT_LABEL)
        style.configure("TEntry",               font=self.FONT_ENTRY)
        style.configure("TSpinbox",             font=self.FONT_ENTRY)
        style.configure("TCombobox",            font=self.FONT_ENTRY)

        if not BACKEND_OK:
            self._show_import_error()
            return

        self._build_ui()
        self.minsize(780, 600)

    # ── Import error screen ───────────────────────────────────────────────────
    def _show_import_error(self):
        f = ttk.Frame(self, padding=30)
        f.pack(fill="both", expand=True)
        ttk.Label(f, text="⚠  Could not import GibsonFinder package",
                  font=("Courier", 11, "bold"), foreground="red").pack(pady=(0,10))
        ttk.Label(f, text=IMPORT_ERROR, font=("Courier", 9),
                  wraplength=600).pack()
        ttk.Label(f, text="\nInstall dependencies:\n  pip install biopython rich click",
                  font=("Courier", 9), foreground="#555").pack(pady=10)

    # ── Build main UI ─────────────────────────────────────────────────────────
    def _build_ui(self):
        # Top: sequence input
        self._build_seq_section()

        # Mode tabs
        nb = ttk.Notebook(self)
        nb.pack(fill="x", padx=10, pady=(0, 6))
        self._tab_single = ttk.Frame(nb, padding=10)
        self._tab_paired = ttk.Frame(nb, padding=10)
        self._tab_dimer  = ttk.Frame(nb, padding=10)
        nb.add(self._tab_single, text="  Single Window  ")
        nb.add(self._tab_paired, text="  Paired Flanking  ")
        nb.add(self._tab_dimer,  text="  Quick Dimer Check  ")
        self._nb = nb

        self._build_single_tab()
        self._build_paired_tab()
        self._build_dimer_tab()

        # Shared parameters
        self._build_params_section()

        # Export section
        self._build_export_section()

        # Run button
        self._build_run_section()

        # Output log
        self._build_log_section()

    # ── Sequence input ────────────────────────────────────────────────────────
    def _build_seq_section(self):
        lf = ttk.LabelFrame(self, text="Sequence Input", padding=8)
        lf.pack(fill="x", padx=10, pady=(10, 6))

        # File picker row
        row = ttk.Frame(lf)
        row.pack(fill="x")
        ttk.Label(row, text="FASTA / sequence file:").pack(side="left")
        self.seq_file_var = tk.StringVar()
        ttk.Entry(row, textvariable=self.seq_file_var, width=55,
                  font=self.FONT_ENTRY).pack(side="left", padx=6)
        ttk.Button(row, text="Browse",
                   command=self._browse_seq).pack(side="left")

        # Or paste directly
        ttk.Label(lf, text="— or paste sequence below (FASTA or raw):").pack(
            anchor="w", pady=(6, 2))
        self.seq_text = tk.Text(lf, height=4, font=self.FONT_MONO,
                                bg="#fffffe", relief="solid", bd=1)
        self.seq_text.pack(fill="x")
        self.seq_text.bind("<KeyRelease>", self._on_seq_change)

        self.seq_info_var = tk.StringVar(value="No sequence loaded.")
        ttk.Label(lf, textvariable=self.seq_info_var,
                  foreground="#666").pack(anchor="w", pady=(2, 0))

    def _browse_seq(self):
        path = filedialog.askopenfilename(
            title="Select sequence file",
            filetypes=[("FASTA / text", "*.fa *.fasta *.fna *.txt"), ("All", "*.*")],
        )
        if path:
            self.seq_file_var.set(path)
            raw = Path(path).read_text(errors="ignore")
            self.seq_text.delete("1.0", "end")
            self.seq_text.insert("1.0", raw[:4000])   # preview
            self._update_seq_info(raw)

    def _on_seq_change(self, *_):
        self._update_seq_info(self.seq_text.get("1.0", "end"))

    def _update_seq_info(self, raw):
        seq = clean_seq(raw)
        if seq:
            from gibsonfinder.thermo import gc_percent
            gc = gc_percent(seq)
            self.seq_info_var.set(
                f"Length: {len(seq):,} bp  ·  GC: {gc}%  ·  "
                f"{'Plasmid-sized' if len(seq)>2000 else 'Fragment'}"
            )
        else:
            self.seq_info_var.set("No valid sequence detected.")

    def _get_seq(self):
        """Return cleaned sequence from file or text box."""
        file_path = self.seq_file_var.get().strip()
        if file_path and Path(file_path).exists():
            raw = Path(file_path).read_text(errors="ignore")
        else:
            raw = self.seq_text.get("1.0", "end")
        return clean_seq(raw)

    # ── Single window tab ─────────────────────────────────────────────────────
    def _build_single_tab(self):
        f = self._tab_single
        ttk.Label(f, text="Scan a single region for both FWD and REV primer candidates.",
                  foreground="#555").grid(row=0, column=0, columnspan=6, sticky="w", pady=(0,8))

        self.s_win_start = tk.IntVar(value=1)
        self.s_win_end   = tk.IntVar(value=500)

        self._lbl_entry(f, "Window Start (bp):", self.s_win_start, 1, 0, width=8)
        self._lbl_entry(f, "Window End (bp):",   self.s_win_end,   1, 2, width=8)

        ttk.Label(f, text="Topology:").grid(row=1, column=4, sticky="e", padx=(16,4))
        self.s_topo = tk.StringVar(value="circular")
        ttk.Combobox(f, textvariable=self.s_topo, values=["circular", "linear"],
                     width=10, state="readonly").grid(row=1, column=5, sticky="w")

    # ── Paired flanking tab ───────────────────────────────────────────────────
    def _build_paired_tab(self):
        f = self._tab_paired
        ttk.Label(f, text="FWD primers from upstream flank  ×  REV primers from downstream flank.",
                  foreground="#555").grid(row=0, column=0, columnspan=6, sticky="w", pady=(0,8))

        self.p_up_start  = tk.IntVar(value=1)
        self.p_up_end    = tk.IntVar(value=150)
        self.p_ins_start = tk.IntVar(value=151)
        self.p_ins_end   = tk.IntVar(value=500)
        self.p_dn_start  = tk.IntVar(value=501)
        self.p_dn_end    = tk.IntVar(value=650)
        self.p_top_n     = tk.IntVar(value=5)

        r = 1
        ttk.Label(f, text="Upstream (FWD):", font=self.FONT_HEAD).grid(
            row=r, column=0, sticky="w", pady=(0,4))
        self._lbl_entry(f, "Start:", self.p_up_start,  r, 1, width=7)
        self._lbl_entry(f, "End:",   self.p_up_end,    r, 3, width=7)

        r = 2
        ttk.Label(f, text="Insert region:", font=self.FONT_HEAD).grid(
            row=r, column=0, sticky="w", pady=(4,4))
        self._lbl_entry(f, "Start:", self.p_ins_start, r, 1, width=7)
        self._lbl_entry(f, "End:",   self.p_ins_end,   r, 3, width=7)
        ttk.Label(f, text="(display only)", foreground="#888").grid(
            row=r, column=5, sticky="w")

        r = 3
        ttk.Label(f, text="Downstream (REV):", font=self.FONT_HEAD).grid(
            row=r, column=0, sticky="w", pady=(4,4))
        self._lbl_entry(f, "Start:", self.p_dn_start,  r, 1, width=7)
        self._lbl_entry(f, "End:",   self.p_dn_end,    r, 3, width=7)

        r = 4
        ttk.Label(f, text="Topology:").grid(row=r, column=0, sticky="e", padx=(0,4), pady=(6,0))
        self.p_topo = tk.StringVar(value="circular")
        ttk.Combobox(f, textvariable=self.p_topo, values=["circular", "linear"],
                     width=10, state="readonly").grid(row=r, column=1, sticky="w", pady=(6,0))
        self._lbl_entry(f, "Top pairs to show:", self.p_top_n, r, 2, width=5, pady=(6,0))

    # ── Quick dimer tab ───────────────────────────────────────────────────────
    def _build_dimer_tab(self):
        f = self._tab_dimer
        ttk.Label(f, text="Paste one or two sequences (5'→3') to check for self- or heterodimer.",
                  foreground="#555").grid(row=0, column=0, columnspan=3, sticky="w", pady=(0,8))

        ttk.Label(f, text="Sequence A (5'→3'):").grid(row=1, column=0, sticky="w")
        self.d_seq_a = ttk.Entry(f, width=55, font=self.FONT_ENTRY)
        self.d_seq_a.grid(row=1, column=1, padx=6, sticky="ew")

        ttk.Label(f, text="Sequence B (5'→3', optional):").grid(row=2, column=0, sticky="w", pady=(4,0))
        self.d_seq_b = ttk.Entry(f, width=55, font=self.FONT_ENTRY)
        self.d_seq_b.grid(row=2, column=1, padx=6, sticky="ew", pady=(4,0))

        ttk.Button(f, text="Check Dimer", command=self._run_dimer_check).grid(
            row=3, column=1, sticky="w", padx=6, pady=8)
        f.columnconfigure(1, weight=1)

    # ── Shared parameters ─────────────────────────────────────────────────────
    def _build_params_section(self):
        lf = ttk.LabelFrame(self, text="Primer Parameters", padding=8)
        lf.pack(fill="x", padx=10, pady=(0, 6))

        # Row 0
        self.bind_min   = tk.IntVar(value=18)
        self.bind_max   = tk.IntVar(value=25)
        self.overlap    = tk.IntVar(value=25)
        self.tm_target  = tk.DoubleVar(value=60.0)
        self.tm_tol     = tk.DoubleVar(value=5.0)
        self.gc_min     = tk.DoubleVar(value=35.0)
        self.gc_max     = tk.DoubleVar(value=65.0)
        self.step       = tk.IntVar(value=3)
        self.max_res    = tk.IntVar(value=100)
        self.show_align = tk.IntVar(value=0)

        row0 = ttk.Frame(lf); row0.pack(fill="x", pady=(0, 4))
        row1 = ttk.Frame(lf); row1.pack(fill="x")

        def sp(parent, label, var, row_frame, width=6, fmt="int"):
            ttk.Label(row_frame, text=label).pack(side="left", padx=(10,2))
            if fmt == "float":
                w = ttk.Spinbox(row_frame, textvariable=var, from_=0, to=999,
                                increment=0.5, width=width, format="%.1f")
            else:
                w = ttk.Spinbox(row_frame, textvariable=var, from_=0, to=9999,
                                increment=1, width=width)
            w.pack(side="left")

        sp(row0, "Bind min (bp):",  self.bind_min,  row0)
        sp(row0, "Bind max (bp):",  self.bind_max,  row0)
        sp(row0, "Overlap (bp):",   self.overlap,   row0)
        sp(row0, "Tm target (°C):", self.tm_target, row0, fmt="float")
        sp(row0, "Tm tol ± (°C):",  self.tm_tol,   row0, fmt="float")

        sp(row1, "GC min (%):",     self.gc_min,    row1, fmt="float")
        sp(row1, "GC max (%):",     self.gc_max,    row1, fmt="float")
        sp(row1, "Step (bp):",      self.step,      row1)
        sp(row1, "Max results:",    self.max_res,   row1)
        ttk.Checkbutton(row1, text="Show dimer alignments",
                        variable=self.show_align).pack(side="left", padx=(16,0))

    # ── Export section ────────────────────────────────────────────────────────
    def _build_export_section(self):
        lf = ttk.LabelFrame(self, text="Export", padding=8)
        lf.pack(fill="x", padx=10, pady=(0, 6))

        row = ttk.Frame(lf); row.pack(fill="x")

        ttk.Label(row, text="Output folder:").pack(side="left")
        self.out_folder_var = tk.StringVar()
        ttk.Entry(row, textvariable=self.out_folder_var, width=40,
                  font=self.FONT_ENTRY).pack(side="left", padx=6)
        ttk.Button(row, text="Browse",
                   command=self._browse_out_folder).pack(side="left")

        row2 = ttk.Frame(lf); row2.pack(fill="x", pady=(6, 0))
        self.exp_csv    = tk.IntVar(value=1)
        self.exp_fasta  = tk.IntVar(value=1)
        self.exp_pairs  = tk.IntVar(value=1)
        self.out_prefix = tk.StringVar(value="gibson_primers")

        ttk.Checkbutton(row2, text="CSV", variable=self.exp_csv).pack(side="left")
        ttk.Checkbutton(row2, text="FASTA", variable=self.exp_fasta).pack(side="left", padx=8)
        ttk.Checkbutton(row2, text="Pairs CSV (paired mode)",
                        variable=self.exp_pairs).pack(side="left")
        ttk.Label(row2, text="  Filename prefix:").pack(side="left", padx=(16,4))
        ttk.Entry(row2, textvariable=self.out_prefix, width=20,
                  font=self.FONT_ENTRY).pack(side="left")

    def _browse_out_folder(self):
        folder = filedialog.askdirectory(title="Select output folder")
        if folder:
            self.out_folder_var.set(folder)

    # ── Run button ────────────────────────────────────────────────────────────
    def _build_run_section(self):
        row = ttk.Frame(self); row.pack(fill="x", padx=10, pady=(0, 6))
        self.run_btn = tk.Button(
            row, text="▶  Find Gibson Primers",
            font=("Courier", 10, "bold"),
            bg="#2a7a4f", fg="white",
            activebackground="#1e5c3a", activeforeground="white",
            relief="flat", padx=18, pady=8,
            cursor="hand2",
            command=self._run,
        )
        self.run_btn.pack(side="left")
        self.status_var = tk.StringVar(value="Ready.")
        ttk.Label(row, textvariable=self.status_var,
                  foreground="#555").pack(side="left", padx=16)

    # ── Output log ────────────────────────────────────────────────────────────
    def _build_log_section(self):
        lf = ttk.LabelFrame(self, text="Output", padding=6)
        lf.pack(fill="both", expand=True, padx=10, pady=(0, 10))

        btn_row = ttk.Frame(lf); btn_row.pack(fill="x", pady=(0,4))
        ttk.Button(btn_row, text="Clear", command=self._clear_log).pack(side="left")

        self.log = scrolledtext.ScrolledText(
            lf, font=self.FONT_MONO, height=18, wrap="word",
            bg="#0d1117", fg="#c9d1d9",
            insertbackground="white", relief="flat",
        )
        self.log.pack(fill="both", expand=True)

        # Text tags for coloring
        self.log.tag_config("green",  foreground="#3fb950")
        self.log.tag_config("cyan",   foreground="#79c0ff")
        self.log.tag_config("yellow", foreground="#e3b341")
        self.log.tag_config("red",    foreground="#f85149")
        self.log.tag_config("dim",    foreground="#8b949e")
        self.log.tag_config("bold",   font=("Courier", 9, "bold"), foreground="#c9d1d9")
        self.log.tag_config("head",   font=("Courier", 9, "bold"), foreground="#79c0ff")
        self.log.tag_config("match",  foreground="#f85149", font=("Courier", 9, "bold"))

    def _clear_log(self):
        self.log.configure(state="normal")
        self.log.delete("1.0", "end")
        self.log.configure(state="disabled")

    # ── Logging helpers ───────────────────────────────────────────────────────
    def _log(self, text, tag=None):
        self.log.configure(state="normal")
        if tag:
            self.log.insert("end", text, tag)
        else:
            self.log.insert("end", text)
        self.log.see("end")
        self.log.configure(state="disabled")

    def _log_rule(self, title=""):
        width = 72
        if title:
            pad = (width - len(title) - 2) // 2
            line = "─" * pad + f" {title} " + "─" * pad
        else:
            line = "─" * width
        self._log("\n" + line + "\n", "head")

    def _log_primer_table(self, primers, params):
        """Render a text table of primers into the log."""
        self._log(
            f"\n{'#':>3}  {'Dir':<4}  {'Pos':>5}  "
            f"{'Tm':>7}  {'GC%':>5}  {'Len':>7}  "
            f"{'Self-dimer':<16}  {'Score':>5}\n", "dim"
        )
        self._log("─" * 72 + "\n", "dim")

        for i, p in enumerate(primers, 1):
            sd = p.self_dimer
            sd_tag = "red" if sd.risk_level == "high" else \
                     "yellow" if sd.risk_level == "low" else "green"
            dir_tag = "green" if p.direction == "fwd" else "cyan"

            self._log(f"{i:>3}  ", "dim")
            self._log(f"{p.direction.upper():<4}  ", dir_tag)
            self._log(f"{p.pos:>5}  ", "dim")
            self._log(f"{p.tm:>6.1f}°C  {p.gc:>4.1f}%  "
                      f"{p.total_len:>3}/{p.bind_len:<3}  ")
            self._log(f"{sd.summary():<16}  ", sd_tag)
            sc_tag = "green" if p.score >= 70 else "yellow" if p.score >= 45 else "red"
            self._log(f"{p.score:>5}\n", sc_tag)

            # Sequence line: overlap in cyan, binding in green
            self._log(f"     5'─ ", "dim")
            self._log(p.overlap_seq, "cyan")
            self._log(p.bind_seq, "green")
            self._log(" ─3'\n")

    def _log_alignment(self, d, label_a="5'→3'", label_b="3'→5'"):
        """Print a dimer alignment block."""
        pad = max(len(label_a), len(label_b)) + 2
        self._log(f"\n{label_a.ljust(pad)}", "dim")
        self._log(d.line_a + "\n", "cyan")
        self._log(" " * pad)
        for ch in d.line_m:
            if ch == "|":
                self._log("|", "match")
            else:
                self._log(" ")
        self._log("\n")
        self._log(f"{label_b.ljust(pad)}", "dim")
        self._log(d.line_b + "\n", "yellow")
        risk_tag = "red" if d.risk_level in ("high",) else \
                   "yellow" if d.risk_level == "low" else "green"
        tp = "  ⚠ 3′ end involved" if d.three_prime_risk else ""
        self._log(f"  Max run: {d.max_run} bp{tp}  [{d.risk_level.upper()}]\n\n", risk_tag)

    def _log_summary(self, primers, params, seq_len, mode):
        from gibsonfinder.thermo import gc_percent
        clean  = sum(1 for p in primers if p.self_dimer.is_clean)
        fwd_n  = sum(1 for p in primers if p.direction == "fwd")
        rev_n  = sum(1 for p in primers if p.direction == "rev")
        self._log_rule("Results Summary")
        self._log(f"  Mode:          {mode}\n")
        self._log(f"  Sequence:      {seq_len:,} bp\n")
        self._log(f"  Candidates:    {len(primers)}\n", "bold")
        self._log(f"  Forward:       {fwd_n}\n", "green")
        self._log(f"  Reverse:       {rev_n}\n", "cyan")
        self._log(f"  Clean SD:      {clean} / {len(primers)}\n",
                  "green" if clean == len(primers) else "yellow")
        self._log(f"  Target Tm:     {params.tm_target}°C ± {params.tm_tol}°C\n")
        self._log(f"  Overlap:       {params.overlap_len} bp\n")
        self._log(f"  Binding:       {params.bind_min}–{params.bind_max} bp\n")

    # ── Parameter builder ─────────────────────────────────────────────────────
    def _make_params(self):
        return ScanParams(
            bind_min    = _int(self.bind_min,  18),
            bind_max    = _int(self.bind_max,  25),
            overlap_len = _int(self.overlap,   25),
            tm_target   = _float(self.tm_target, 60.0),
            tm_tol      = _float(self.tm_tol,    5.0),
            gc_min      = _float(self.gc_min,   35.0),
            gc_max      = _float(self.gc_max,   65.0),
            step        = _int(self.step,        3),
            circular    = True,   # handled per-mode below
            max_results = _int(self.max_res, 100) or None,
        )

    # ── Run dispatcher ────────────────────────────────────────────────────────
    def _run(self):
        tab = self._nb.index(self._nb.select())
        if tab == 2:   # dimer tab has its own button
            self._run_dimer_check()
            return
        seq = self._get_seq()
        if len(seq) < 30:
            messagebox.showerror("Error", "No valid sequence found (need ≥ 30 bp ATCG).")
            return

        self.run_btn.configure(state="disabled", text="Running…")
        self.status_var.set("Scanning…")
        self.update_idletasks()

        thread = threading.Thread(target=self._run_in_thread,
                                  args=(tab, seq), daemon=True)
        thread.start()

    def _run_in_thread(self, tab, seq):
        try:
            if tab == 0:
                self._do_single(seq)
            elif tab == 1:
                self._do_paired(seq)
        except Exception as e:
            self.after(0, lambda: messagebox.showerror("Error", str(e)))
        finally:
            self.after(0, self._run_done)

    def _run_done(self):
        self.run_btn.configure(state="normal", text="▶  Find Gibson Primers")

    # ── Single window run ─────────────────────────────────────────────────────
    def _do_single(self, seq):
        params = self._make_params()
        params.circular = self.s_topo.get() == "circular"
        win_start = _int(self.s_win_start, 1)
        win_end   = min(_int(self.s_win_end, 500), len(seq))

        self.after(0, lambda: self.status_var.set("Scanning window…"))
        primers = scan_single(seq, win_start, win_end, params)

        def update():
            if not primers:
                self._log("\n⚠  No primers found. Try relaxing Tm tolerance or GC limits.\n", "yellow")
                self.status_var.set("No results.")
                return
            self.status_var.set(f"Done. {len(primers)} primer(s) found.")
            self._log_summary(primers, params, len(seq), "single")
            self._log_rule(f"Single Window — {win_start}:{win_end}")
            self._log_primer_table(primers, params)
            self._log_dimer_section(primers, params)
            self._do_exports_single(primers, params)

        self.after(0, update)

    # ── Paired run ────────────────────────────────────────────────────────────
    def _do_paired(self, seq):
        params = self._make_params()
        params.circular = self.p_topo.get() == "circular"
        top_n = _int(self.p_top_n, 5)

        up_s = _int(self.p_up_start,  1)
        up_e = _int(self.p_up_end,  150)
        dn_s = _int(self.p_dn_start, 501)
        dn_e = _int(self.p_dn_end,  650)

        self.after(0, lambda: self.status_var.set("Scanning flanking windows…"))
        fwds, revs = scan_paired(seq, up_s, up_e, dn_s, dn_e, params)

        def update():
            all_p = fwds + revs
            if not all_p:
                self._log("\n⚠  No primers found. Try relaxing parameters or widening windows.\n", "yellow")
                self.status_var.set("No results.")
                return
            self.status_var.set(f"Done. {len(fwds)} FWD  +  {len(revs)} REV primers found.")
            self._log_summary(all_p, params, len(seq), "paired")

            self._log_rule(f"Upstream FWD — {up_s}:{up_e}")
            self._log_primer_table(fwds, params)

            self._log_rule(f"Downstream REV — {dn_s}:{dn_e}")
            self._log_primer_table(revs, params)

            self._log_paired_analysis(fwds, revs, params, top_n)
            self._log_dimer_section(all_p, params)
            self._do_exports_paired(fwds, revs, params, top_n)

        self.after(0, update)

    # ── Paired analysis log ───────────────────────────────────────────────────
    def _log_paired_analysis(self, fwds, revs, params, top_n):
        self._log_rule("Paired Analysis — FWD × REV combinations")
        fw_top = fwds[:top_n]
        rv_top = revs[:top_n]
        if not fw_top or not rv_top:
            self._log("  Not enough primers in one or both windows.\n", "yellow")
            return

        pair_num = 0
        for fwd in fw_top:
            for rev in rv_top:
                pair_num += 1
                hd       = dimer_analysis(fwd.full_primer, rev.full_primer)
                tm_diff  = abs(fwd.tm - rev.tm)
                ps       = round((fwd.score + rev.score) / 2)
                hd_tag   = "red" if hd.risk_level == "high" else \
                           "yellow" if hd.risk_level == "low" else "green"
                tm_tag   = "red" if tm_diff > 5 else \
                           "yellow" if tm_diff > 3 else "green"

                self._log(f"\n  ┌─ Pair #{pair_num}  FWD pos {fwd.pos}  ×  REV pos {rev.pos}", "bold")
                self._log(f"  ΔTm: ", "dim")
                self._log(f"{tm_diff:.1f}°C", tm_tag)
                self._log(f"  PairScore: {ps}  Heterodimer: ", "dim")
                self._log(f"{hd.summary()}\n", hd_tag)

                # FWD
                sd_f_tag = "red" if fwd.self_dimer.risk_level == "high" else \
                           "yellow" if fwd.self_dimer.risk_level == "low" else "green"
                self._log(f"  FWD  5'─ ", "dim")
                self._log(fwd.overlap_seq, "cyan")
                self._log(fwd.bind_seq, "green")
                self._log(" ─3'")
                self._log(f"  Tm {fwd.tm}°C  GC {fwd.gc}%  Score {fwd.score}"
                          f"  SD: ", "dim")
                self._log(f"{fwd.self_dimer.summary()}\n", sd_f_tag)

                # REV
                sd_r_tag = "red" if rev.self_dimer.risk_level == "high" else \
                           "yellow" if rev.self_dimer.risk_level == "low" else "green"
                self._log(f"  REV  5'─ ", "dim")
                self._log(rev.overlap_seq, "cyan")
                self._log(rev.bind_seq, "green")
                self._log(" ─3'")
                self._log(f"  Tm {rev.tm}°C  GC {rev.gc}%  Score {rev.score}"
                          f"  SD: ", "dim")
                self._log(f"{rev.self_dimer.summary()}\n", sd_r_tag)

                # Show self-dimer alignments if flagged
                if not fwd.self_dimer.is_clean:
                    self._log("  FWD self-dimer:\n", "yellow")
                    self._log_alignment(fwd.self_dimer,
                                        "  FWD 5'→3' ", "  FWD 3'→5' ")
                if not rev.self_dimer.is_clean:
                    self._log("  REV self-dimer:\n", "yellow")
                    self._log_alignment(rev.self_dimer,
                                        "  REV 5'→3' ", "  REV 3'→5' ")

                # Heterodimer
                if hd.risk_level != "none":
                    self._log("  Heterodimer alignment:\n", hd_tag)
                    self._log_alignment(hd, "  FWD 5'→3' ", "  REV 3'→5' ")
                else:
                    self._log("  ✓ No significant heterodimer detected\n", "green")

    # ── Dimer section ─────────────────────────────────────────────────────────
    def _log_dimer_section(self, primers, params):
        show_align = _bool(self.show_align)
        self._log_rule("Dimer Analysis")
        self._log(
            f"  {'Type':<7}  {'Primer A':<22}  {'Primer B':<22}  "
            f"{'Max run':>8}  {'3′ end':>6}  {'Risk'}\n", "dim"
        )
        self._log("  " + "─" * 68 + "\n", "dim")

        top = primers[:20]
        rows = []
        for p in primers[:30]:
            rows.append(("Self", p, p, p.self_dimer))
        for i, pa in enumerate(top):
            for pb in top[i+1:]:
                hd = dimer_analysis(pa.full_primer, pb.full_primer)
                rows.append(("Hetero", pa, pb, hd))
        rows.sort(key=lambda r: r[3].max_run, reverse=True)

        for kind, pa, pb, d in rows:
            rc = "red" if d.risk_level == "high" else \
                 "yellow" if d.risk_level == "low" else "green"
            b_name = "(self)" if kind == "Self" else pb.label
            self._log(f"  {kind:<7}  {pa.label:<22}  {b_name:<22}  "
                      f"{d.max_run:>5} bp  "
                      f"{'⚠ Yes' if d.three_prime_risk else 'No':>6}  ", "dim")
            self._log(f"{d.risk_level.upper()}\n", rc)

            if show_align and d.max_run > 0:
                la = pa.label + " 5'→3'"
                lb = ("(self) 3'→5'" if kind == "Self"
                      else pb.label + " 3'→5'")
                self._log_alignment(d, "  " + la, "  " + lb)

    # ── Dimer quick check ─────────────────────────────────────────────────────
    def _run_dimer_check(self):
        seq_a = self.d_seq_a.get().strip().upper()
        seq_b = self.d_seq_b.get().strip().upper()
        if not seq_a:
            messagebox.showerror("Error", "Enter at least Sequence A.")
            return

        from gibsonfinder.dimers import self_dimer as sd_fn
        self._log_rule(f"Quick Dimer Check")
        self._log(f"  Seq A: {seq_a}  ({len(seq_a)} bp)\n", "cyan")
        d = sd_fn(seq_a)
        self._log("  Self-dimer:\n", "bold")
        self._log_alignment(d, "  5'→3' ", "  3'→5' ")

        if seq_b:
            self._log(f"  Seq B: {seq_b}  ({len(seq_b)} bp)\n", "cyan")
            hd = dimer_analysis(seq_a, seq_b)
            self._log("  Heterodimer:\n", "bold")
            self._log_alignment(hd, "  A 5'→3' ", "  B 3'→5' ")

    # ── Exports ───────────────────────────────────────────────────────────────
    def _get_out_path(self, suffix):
        folder = self.out_folder_var.get().strip()
        if not folder:
            return None
        prefix = self.out_prefix.get().strip() or "gibson_primers"
        return Path(folder) / f"{prefix}{suffix}"

    def _do_exports_single(self, primers, params):
        saved = []
        if _bool(self.exp_csv):
            p = self._get_out_path(".csv")
            if p:
                to_csv(primers, p)
                saved.append(str(p))
        if _bool(self.exp_fasta):
            p = self._get_out_path(".fasta")
            if p:
                to_fasta(primers, p)
                saved.append(str(p))
        if saved:
            self._log_rule("Exports")
            for s in saved:
                self._log(f"  Saved → {s}\n", "green")

    def _do_exports_paired(self, fwds, revs, params, top_n):
        saved = []
        all_p = fwds + revs
        if _bool(self.exp_csv):
            p = self._get_out_path(".csv")
            if p:
                to_csv(all_p, p)
                saved.append(str(p))
        if _bool(self.exp_fasta):
            p = self._get_out_path(".fasta")
            if p:
                to_fasta(all_p, p)
                saved.append(str(p))
        if _bool(self.exp_pairs):
            p = self._get_out_path("_pairs.csv")
            if p:
                pairs_to_csv(fwds, revs, top_n=top_n, path=p)
                saved.append(str(p))
        if saved:
            self._log_rule("Exports")
            for s in saved:
                self._log(f"  Saved → {s}\n", "green")

    # ── Misc layout helper ────────────────────────────────────────────────────
    def _lbl_entry(self, parent, text, var, row, col,
                   width=8, pady=0):
        ttk.Label(parent, text=text).grid(
            row=row, column=col, sticky="e", padx=(12, 4), pady=pady)
        ttk.Spinbox(parent, textvariable=var, from_=0, to=999999,
                    increment=1, width=width).grid(
            row=row, column=col+1, sticky="w", pady=pady)


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    app = GibsonFinderApp()
    app.mainloop()
