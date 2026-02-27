"""
POET V2.0b — Graphical User Interface

A tkinter-based GUI for configuring, running, and visualising POET experiments.

Features
--------
- Full configuration editor with grouped settings
- Save / load config files
- One-click Run with live console output
- Live fitness chart updated every generation
- Experiment mode: pick multiple configs, set replicates & generations, run comparison
- Stop button to abort a running experiment

Launch
------
    python poet_gui.py
"""

import copy
import math
import os
import queue
import re
import sys
import threading
import time
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# ---------------------------------------------------------------------------
# POET imports
# ---------------------------------------------------------------------------
import eletility

# ---------------------------------------------------------------------------
# Constants / layout helpers
# ---------------------------------------------------------------------------
PAD = 6
ENTRY_W = 22
LABEL_W = 30

# Organised config sections — each tuple is (key, label, widget_type, tooltip)
# widget_type: "entry" | "bool" | "choice:a,b,c" | "file"
CONFIG_SECTIONS: dict[str, list[tuple[str, str, str, str]]] = {
    "General": [
        ("seed", "Seed", "entry", "Random seed for reproducibility"),
        ("learn_data", "Training data", "file", "CSV with sequences and fitness"),
        ("unseen_data", "Test data (optional)", "file", "Held-out test CSV for CGP"),
        ("fitness_alg", "Fitness algorithm", "choice:correlation,RMSE", ""),
        (
            "diversity_selection",
            "Diversity selection",
            "bool",
            "Cluster-based diversity in selection",
        ),
        (
            "diversity_weight",
            "Diversity weight",
            "entry",
            "Weight for diversity term (0–1)",
        ),
    ],
    "Population": [
        ("population_size", "Population size", "entry", ""),
        ("runs", "Generations (runs)", "entry", "Number of evolutionary generations"),
        ("tournament_size", "Tournament size", "entry", ""),
        ("maximum_rule_size", "Max rule size", "entry", "Max pattern length per rule"),
        ("maximum_rule_count", "Max rule count", "entry", "Max rules per individual"),
        ("rule_weight_min", "Weight min", "entry", ""),
        ("rule_weight_max", "Weight max", "entry", ""),
        (
            "crossover_unused_selection_chance",
            "Crossover unused chance",
            "entry",
            "Chance to carry unused rules during crossover",
        ),
        ("workers", "Workers (MP)", "entry", "Parallel worker processes (0 = off)"),
    ],
    "Output": [
        (
            "output_evo",
            "Evolution log CSV",
            "entry",
            "Path for per-generation stats CSV",
        ),
        ("output_model", "Model CSV", "entry", "Path for best model output"),
    ],
    "Mutation — Model": [
        ("mut_add_rule", "Add rule", "entry", ""),
        ("mut_remove_rule", "Remove rule", "entry", ""),
    ],
    "Mutation — Rule": [
        ("mut_change_weight", "Change weight", "entry", ""),
        ("mut_add_to_pattern", "Add to pattern", "entry", ""),
        ("mut_remove_from_pattern", "Remove from pattern", "entry", ""),
        ("mut_change_character", "Change character", "entry", ""),
        ("mut_insert_gap", "Insert gap", "entry", ""),
        ("mut_remove_gap", "Remove gap", "entry", ""),
    ],
    "Labels": [
        (
            "experiment_label",
            "Experiment label",
            "entry",
            "Label for this config in experiment plots",
        ),
    ],
    "Advanced": [
        ("parsimony_pressure", "Parsimony pressure", "entry", "Penalty per rule"),
        (
            "optimize_weights",
            "OLS weight optimisation",
            "bool",
            "Analytical weight setting via least squares",
        ),
        ("enable_gaps", "Enable gaps", "bool", "Allow gapped motifs (A___B)"),
        (
            "matching_mode",
            "Matching mode",
            "choice:substring,regex",
            "Rule representation",
        ),
    ],
    "Regex Mode": [
        ("regex_alphabet", "Regex alphabet", "file", ""),
        ("max_depth_tree", "Max tree depth", "entry", ""),
        ("min_braces", "Min braces", "entry", ""),
        ("max_braces", "Max braces", "entry", ""),
        ("init_method", "Init method", "choice:half,full,grow", ""),
        ("mut_replace_rule", "Replace rule", "entry", ""),
        ("mut_replace_subtree", "Replace subtree", "entry", ""),
        ("mut_add_aa", "Add AA", "entry", ""),
        ("mut_replace_node", "Replace node", "entry", ""),
    ],
    "Experimental Features": [
        ("exp_char_classes", "Character classes", "bool", "e.g. [AG]K"),
        (
            "exp_variable_gaps",
            "Variable-length gaps",
            "bool",
            "e.g. A_{2,5}B",
        ),
        (
            "exp_weighted_positions",
            "Weighted positions",
            "bool",
            "Per-position importance weights",
        ),
        ("exp_pw_threshold", "PW threshold", "entry", "Min match quality (0–1)"),
        ("exp_composition", "Rule composition", "bool", "AND/OR groups of rules"),
        (
            "exp_match_count",
            "Match count scaling",
            "bool",
            "Weight proportional to match count",
        ),
        ("exp_circular", "Circular matching", "bool", "Wrap-around sequences"),
        ("max_variable_gap", "Max variable gap", "entry", ""),
    ],
    "Exp. Mutation Rates": [
        ("mut_char_class", "Char class", "entry", ""),
        ("mut_variable_gap", "Variable gap", "entry", ""),
        ("mut_position_weight", "Position weight", "entry", ""),
        ("mut_composition", "Composition", "entry", ""),
    ],
    "CGP Post-Processing": [
        ("use_cgp", "Enable CGP", "bool", "Run CGP after POET evolution"),
        ("cgp_generations", "CGP generations", "entry", ""),
        ("cgp_model_size", "Model size", "entry", ""),
        ("cgp_parents", "Parents", "entry", ""),
        ("cgp_children", "Children", "entry", ""),
        (
            "cgp_mutation_type",
            "Mutation type",
            "choice:full,point",
            "",
        ),
        ("cgp_mutation_rate", "Mutation rate", "entry", ""),
        (
            "cgp_selection_type",
            "Selection type",
            "choice:paretoelite,elite,tournament",
            "",
        ),
        (
            "cgp_fitness_function",
            "Fitness function",
            "choice:correlation,correlation_complexity",
            "",
        ),
        ("cgp_n_elites", "N elites", "entry", ""),
        ("cgp_tournament_size", "Tournament size", "entry", ""),
        ("cgp_step_size", "Step size", "entry", "Reporting interval"),
    ],
}

# Default config values (used when creating a fresh config)
DEFAULT_CONFIG: dict[str, str] = {
    "seed": "333",
    "learn_data": "data/uPAR-binding-peptides-epoch0.csv",
    "unseen_data": "",
    "fitness_alg": "correlation",
    "diversity_selection": "True",
    "diversity_weight": "0.05",
    "population_size": "100",
    "runs": "50",
    "tournament_size": "5",
    "maximum_rule_size": "6",
    "maximum_rule_count": "100",
    "rule_weight_min": "0.0",
    "rule_weight_max": "10.0",
    "crossover_unused_selection_chance": "0.2",
    "workers": "0",
    "experiment_label": "",
    "output_evo": "output/evo.csv",
    "output_model": "output/model.csv",
    "mut_add_rule": "0.2",
    "mut_remove_rule": "0.2",
    "mut_change_weight": "0.2",
    "mut_add_to_pattern": "0.1",
    "mut_remove_from_pattern": "0.1",
    "mut_change_character": "0.1",
    "mut_insert_gap": "0.1",
    "mut_remove_gap": "0.1",
    "parsimony_pressure": "0.0",
    "optimize_weights": "False",
    "enable_gaps": "True",
    "matching_mode": "substring",
    "regex_alphabet": "data/translation/regex_alphabet.csv",
    "max_depth_tree": "4",
    "min_braces": "1",
    "max_braces": "3",
    "init_method": "half",
    "mut_replace_rule": "0.1",
    "mut_replace_subtree": "0.1",
    "mut_add_aa": "0.1",
    "mut_replace_node": "0.1",
    "exp_char_classes": "False",
    "exp_variable_gaps": "False",
    "exp_weighted_positions": "False",
    "exp_pw_threshold": "0.7",
    "exp_composition": "False",
    "exp_match_count": "False",
    "exp_circular": "False",
    "max_variable_gap": "6",
    "mut_char_class": "0.1",
    "mut_variable_gap": "0.1",
    "mut_position_weight": "0.1",
    "mut_composition": "0.05",
    "use_cgp": "False",
    "cgp_generations": "1000",
    "cgp_model_size": "32",
    "cgp_parents": "1",
    "cgp_children": "4",
    "cgp_mutation_type": "full",
    "cgp_mutation_rate": "1.0",
    "cgp_selection_type": "paretoelite",
    "cgp_fitness_function": "correlation",
    "cgp_n_elites": "1",
    "cgp_tournament_size": "4",
    "cgp_step_size": "10",
}


# ═══════════════════════════════════════════════════════════════════════════
# Helper: redirect stdout/stderr to a queue for the console panel
# ═══════════════════════════════════════════════════════════════════════════
class _QueueWriter:
    """File-like object that puts written text onto a queue."""

    def __init__(self, q: queue.Queue):
        self._q = q

    def write(self, text: str):
        if text:
            self._q.put(text)

    def flush(self):
        pass


# ═══════════════════════════════════════════════════════════════════════════
# Main Application
# ═══════════════════════════════════════════════════════════════════════════
class POETApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("POET V2.0b")
        self.geometry("1400x900")
        self.minsize(1100, 700)

        # State
        self._config: dict[str, str] = copy.deepcopy(DEFAULT_CONFIG)
        self._widgets: dict[str, tk.Widget] = {}  # key → entry/combo/checkbutton
        self._vars: dict[str, tk.Variable] = {}  # key → StringVar / BooleanVar
        self._run_thread: threading.Thread | None = None
        self._stop_event = threading.Event()
        self._output_queue: queue.Queue = queue.Queue()
        self._evo_path: str = ""  # path to the evo.csv being written

        self._build_ui()

        # Auto-load config.ini if present
        if os.path.exists("config.ini"):
            try:
                parser = eletility.ConfigParser()
                self._config = parser.read("config.ini")
                self._status_var.set("Loaded config.ini")
            except Exception:
                pass  # fall back to defaults

        self._load_config_into_ui(self._config)
        self._poll_queue()

    # -------------------------------------------------------------------
    # UI construction
    # -------------------------------------------------------------------
    def _build_ui(self):
        # Main paned window: left (config) | right (run + chart)
        self._panes = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        self._panes.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        # ── Left: config editor ──
        left_frame = ttk.Frame(self._panes)
        self._panes.add(left_frame, weight=1)

        # Toolbar
        tb = ttk.Frame(left_frame)
        tb.pack(fill=tk.X, padx=PAD, pady=(PAD, 0))
        ttk.Button(tb, text="Load Config", command=self._on_load_config).pack(
            side=tk.LEFT, padx=2
        )
        ttk.Button(tb, text="Save Config", command=self._on_save_config).pack(
            side=tk.LEFT, padx=2
        )
        ttk.Button(tb, text="Reset Defaults", command=self._on_reset).pack(
            side=tk.LEFT, padx=2
        )
        ttk.Separator(tb, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=6, pady=2
        )
        ttk.Button(
            tb, text="Run Single", command=self._on_run_single, style="Accent.TButton"
        ).pack(side=tk.LEFT, padx=2)
        ttk.Button(tb, text="Stop", command=self._on_stop).pack(side=tk.LEFT, padx=2)

        # Notebook (tabbed config sections)
        nb = ttk.Notebook(left_frame)
        nb.pack(fill=tk.BOTH, expand=True, padx=PAD, pady=PAD)

        for section_name, fields in CONFIG_SECTIONS.items():
            page = ttk.Frame(nb)
            nb.add(page, text=section_name)

            # Scrollable
            canvas = tk.Canvas(page, highlightthickness=0)
            scrollbar = ttk.Scrollbar(page, orient=tk.VERTICAL, command=canvas.yview)
            inner = ttk.Frame(canvas)
            inner.bind(
                "<Configure>",
                lambda e, c=canvas: c.configure(scrollregion=c.bbox("all")),
            )
            canvas.create_window((0, 0), window=inner, anchor=tk.NW)
            canvas.configure(yscrollcommand=scrollbar.set)
            canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

            # Mouse-wheel scrolling — scoped per canvas via Enter/Leave
            def _bind_wheel(event, c=canvas):
                c.bind_all(
                    "<MouseWheel>",
                    lambda e, _c=c: _c.yview_scroll(int(-1 * (e.delta / 120)), "units"),
                )

            def _unbind_wheel(event, c=canvas):
                c.unbind_all("<MouseWheel>")

            canvas.bind("<Enter>", _bind_wheel)
            canvas.bind("<Leave>", _unbind_wheel)

            for row_idx, (key, label, wtype, tooltip) in enumerate(fields):
                ttk.Label(inner, text=label, width=LABEL_W, anchor=tk.W).grid(
                    row=row_idx, column=0, sticky=tk.W, padx=4, pady=2
                )

                if wtype == "bool":
                    var = tk.BooleanVar(value=False)
                    w = ttk.Checkbutton(inner, variable=var)
                    w.grid(row=row_idx, column=1, sticky=tk.W, padx=4, pady=2)
                    self._vars[key] = var
                    self._widgets[key] = w

                elif wtype.startswith("choice:"):
                    choices = wtype.split(":", 1)[1].split(",")
                    var = tk.StringVar(value=choices[0])
                    w = ttk.Combobox(
                        inner,
                        textvariable=var,
                        values=choices,
                        width=ENTRY_W,
                        state="readonly",
                    )
                    w.grid(row=row_idx, column=1, sticky=tk.W, padx=4, pady=2)
                    self._vars[key] = var
                    self._widgets[key] = w

                elif wtype == "file":
                    var = tk.StringVar()
                    frm = ttk.Frame(inner)
                    frm.grid(row=row_idx, column=1, sticky=tk.EW, padx=4, pady=2)
                    e = ttk.Entry(frm, textvariable=var, width=ENTRY_W)
                    e.pack(side=tk.LEFT, fill=tk.X, expand=True)
                    ttk.Button(
                        frm,
                        text="...",
                        width=3,
                        command=lambda v=var: v.set(
                            filedialog.askopenfilename(
                                filetypes=[("CSV", "*.csv"), ("All", "*.*")]
                            )
                            or v.get()
                        ),
                    ).pack(side=tk.LEFT, padx=2)
                    self._vars[key] = var
                    self._widgets[key] = e

                else:  # "entry"
                    var = tk.StringVar()
                    w = ttk.Entry(inner, textvariable=var, width=ENTRY_W)
                    w.grid(row=row_idx, column=1, sticky=tk.W, padx=4, pady=2)
                    self._vars[key] = var
                    self._widgets[key] = w

                if tooltip:
                    target = self._widgets[key]
                    _create_tooltip(target, tooltip)

        # ── Experiment tab ──
        exp_page = ttk.Frame(nb)
        nb.add(exp_page, text="Experiments")
        self._build_experiment_tab(exp_page)

        # ── Right: output + chart ──
        right_frame = ttk.Frame(self._panes)
        self._panes.add(right_frame, weight=2)

        # Chart
        chart_frame = ttk.LabelFrame(right_frame, text="Fitness Over Generations")
        chart_frame.pack(fill=tk.BOTH, expand=True, padx=PAD, pady=(PAD, 0))

        self._fig, self._ax = plt.subplots(figsize=(7, 3.5), dpi=90)
        self._fig.subplots_adjust(left=0.10, right=0.97, top=0.92, bottom=0.15)
        self._canvas_chart = FigureCanvasTkAgg(self._fig, master=chart_frame)
        self._canvas_chart.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Console output
        console_frame = ttk.LabelFrame(right_frame, text="Console Output")
        console_frame.pack(fill=tk.BOTH, expand=True, padx=PAD, pady=PAD)

        self._console = tk.Text(
            console_frame,
            wrap=tk.WORD,
            font=("Consolas", 9),
            bg="#1e1e1e",
            fg="#d4d4d4",
            insertbackground="#d4d4d4",
            state=tk.DISABLED,
            height=12,
        )
        console_sb = ttk.Scrollbar(
            console_frame, orient=tk.VERTICAL, command=self._console.yview
        )
        self._console.configure(yscrollcommand=console_sb.set)
        self._console.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        console_sb.pack(side=tk.RIGHT, fill=tk.Y)

        # Status bar
        self._status_var = tk.StringVar(value="Ready")
        ttk.Label(self, textvariable=self._status_var, anchor=tk.W).pack(
            fill=tk.X, padx=PAD, pady=(0, 2)
        )

    # -------------------------------------------------------------------
    # Experiment tab
    # -------------------------------------------------------------------
    def _build_experiment_tab(self, parent: ttk.Frame):
        top = ttk.Frame(parent)
        top.pack(fill=tk.X, padx=PAD, pady=PAD)

        ttk.Label(top, text="Config files for comparison:").pack(anchor=tk.W)

        list_frame = ttk.Frame(top)
        list_frame.pack(fill=tk.BOTH, expand=True, pady=4)

        self._exp_listbox = tk.Listbox(list_frame, height=8, selectmode=tk.EXTENDED)
        self._exp_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        exp_sb = ttk.Scrollbar(
            list_frame, orient=tk.VERTICAL, command=self._exp_listbox.yview
        )
        exp_sb.pack(side=tk.RIGHT, fill=tk.Y)
        self._exp_listbox.configure(yscrollcommand=exp_sb.set)

        btn_frame = ttk.Frame(top)
        btn_frame.pack(fill=tk.X)
        ttk.Button(btn_frame, text="Add Config(s)", command=self._exp_add).pack(
            side=tk.LEFT, padx=2
        )
        ttk.Button(btn_frame, text="Remove Selected", command=self._exp_remove).pack(
            side=tk.LEFT, padx=2
        )
        ttk.Button(btn_frame, text="Clear All", command=self._exp_clear).pack(
            side=tk.LEFT, padx=2
        )

        params = ttk.Frame(top)
        params.pack(fill=tk.X, pady=6)

        ttk.Label(params, text="Replicates:").grid(row=0, column=0, padx=4)
        self._exp_reps_var = tk.StringVar(value="10")
        ttk.Entry(params, textvariable=self._exp_reps_var, width=6).grid(
            row=0, column=1, padx=4
        )

        ttk.Label(params, text="Generations:").grid(row=0, column=2, padx=4)
        self._exp_gens_var = tk.StringVar(value="50")
        ttk.Entry(params, textvariable=self._exp_gens_var, width=6).grid(
            row=0, column=3, padx=4
        )

        ttk.Label(params, text="Workers:").grid(row=0, column=4, padx=4)
        self._exp_workers_var = tk.StringVar(value="1")
        ttk.Entry(params, textvariable=self._exp_workers_var, width=6).grid(
            row=0, column=5, padx=4
        )

        ttk.Button(top, text="Run Experiment", command=self._on_run_experiment).pack(
            anchor=tk.W, pady=4
        )

        # Preset buttons for common configs
        presets = ttk.LabelFrame(parent, text="Quick Presets")
        presets.pack(fill=tk.X, padx=PAD, pady=PAD)

        ttk.Button(
            presets,
            text="Substring vs Regex",
            command=lambda: self._exp_set_preset(
                [
                    "configs/experiment_substring.ini",
                    "configs/experiment_regex.ini",
                ]
            ),
        ).pack(side=tk.LEFT, padx=2, pady=2)

        ttk.Button(
            presets,
            text="All Experimental Features",
            command=lambda: self._exp_set_preset(
                [
                    "configs/exp_gaps_only.ini",
                    "configs/exp_char_classes.ini",
                    "configs/exp_variable_gaps.ini",
                    "configs/exp_charclass_vargaps.ini",
                    "configs/exp_cc_vargaps_ols.ini",
                    "configs/exp_weighted_pos.ini",
                    "configs/exp_match_count.ini",
                    "configs/exp_circular.ini",
                    "configs/exp_composition.ini",
                    "configs/exp_all_features.ini",
                ]
            ),
        ).pack(side=tk.LEFT, padx=2, pady=2)

        ttk.Button(
            presets,
            text="CC + VarGaps + OLS",
            command=lambda: self._exp_set_preset(
                [
                    "configs/exp_gaps_only.ini",
                    "configs/exp_charclass_vargaps.ini",
                    "configs/exp_cc_vargaps_ols.ini",
                ]
            ),
        ).pack(side=tk.LEFT, padx=2, pady=2)

    def _exp_add(self):
        files = filedialog.askopenfilenames(
            title="Select experiment config files",
            filetypes=[("INI", "*.ini"), ("All", "*.*")],
            initialdir="configs",
        )
        for f in files:
            self._exp_listbox.insert(tk.END, f)

    def _exp_remove(self):
        for idx in reversed(self._exp_listbox.curselection()):
            self._exp_listbox.delete(idx)

    def _exp_clear(self):
        self._exp_listbox.delete(0, tk.END)

    def _exp_set_preset(self, paths: list[str]):
        self._exp_listbox.delete(0, tk.END)
        for p in paths:
            self._exp_listbox.insert(tk.END, p)

    # -------------------------------------------------------------------
    # Config I/O
    # -------------------------------------------------------------------
    def _config_from_ui(self) -> dict[str, str]:
        """Read all UI widgets back into a flat config dict."""
        cfg: dict[str, str] = {}
        for key, var in self._vars.items():
            if isinstance(var, tk.BooleanVar):
                cfg[key] = str(var.get())
            else:
                cfg[key] = var.get()
        return cfg

    def _load_config_into_ui(self, cfg: dict[str, str]):
        """Push a config dict into all UI widgets."""
        for key, var in self._vars.items():
            val = cfg.get(key, DEFAULT_CONFIG.get(key, ""))
            if isinstance(var, tk.BooleanVar):
                var.set(val.strip().lower() in ("true", "1", "yes"))
            else:
                var.set(val)

    def _on_load_config(self):
        path = filedialog.askopenfilename(
            title="Load config file",
            filetypes=[("INI", "*.ini"), ("All", "*.*")],
            initialdir=".",
        )
        if not path:
            return
        try:
            parser = eletility.ConfigParser()
            cfg = parser.read(path)
            self._config = cfg
            self._load_config_into_ui(cfg)
            self._status_var.set(f"Loaded: {path}")
        except Exception as e:
            messagebox.showerror("Load Error", str(e))

    def _on_save_config(self):
        path = filedialog.asksaveasfilename(
            title="Save config file",
            defaultextension=".ini",
            filetypes=[("INI", "*.ini"), ("All", "*.*")],
            initialdir=".",
        )
        if not path:
            return
        cfg = self._config_from_ui()
        try:
            self._write_config_file(path, cfg)
            self._status_var.set(f"Saved: {path}")
        except Exception as e:
            messagebox.showerror("Save Error", str(e))

    @staticmethod
    def _write_config_file(path: str, cfg: dict[str, str]):
        """Write a flat dict to an INI-style file with section comments."""
        section_keys: dict[str, list[str]] = {}
        for section_name, fields in CONFIG_SECTIONS.items():
            section_keys[section_name] = [f[0] for f in fields]

        written_keys: set[str] = set()
        with open(path, "w") as f:
            for section_name, keys in section_keys.items():
                f.write(f"# {section_name}\n")
                for key in keys:
                    val = cfg.get(key, "")
                    if val:
                        f.write(f"{key} = {val}\n")
                    written_keys.add(key)
                f.write("\n")

            # Write any remaining keys not in our sections
            for key, val in cfg.items():
                if key not in written_keys and val:
                    f.write(f"{key} = {val}\n")

    def _on_reset(self):
        self._config = copy.deepcopy(DEFAULT_CONFIG)
        self._load_config_into_ui(self._config)
        self._status_var.set("Defaults restored")

    # -------------------------------------------------------------------
    # Run: single POET run
    # -------------------------------------------------------------------
    def _on_run_single(self):
        if self._run_thread and self._run_thread.is_alive():
            messagebox.showwarning("Busy", "A run is already in progress.")
            return

        cfg = self._config_from_ui()
        self._stop_event.clear()
        self._clear_console()
        self._clear_chart()

        self._evo_path = cfg.get("output_evo", "output/evo.csv")
        self._status_var.set("Running POET...")

        self._run_thread = threading.Thread(
            target=self._run_poet_thread, args=(cfg,), daemon=True
        )
        self._run_thread.start()
        self._poll_chart()

    def _run_poet_thread(self, cfg: dict[str, str]):
        """Run POET in a background thread, capturing stdout."""
        import importlib
        import random as rand

        import archivist
        import optimizer
        import pop as population

        # Reload modules so re-runs start clean
        importlib.reload(archivist)
        importlib.reload(population)
        importlib.reload(optimizer)

        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = _QueueWriter(self._output_queue)
        sys.stderr = _QueueWriter(self._output_queue)

        try:
            print("=" * 50)
            print("POET V2.0b — GUI Run")
            print("=" * 50)
            print()

            rand.seed(int(cfg["seed"]))

            arch = archivist.Archivist(cfg)
            arch.setup()

            pop = population.Population(cfg)
            opt = optimizer.Optimizer(cfg, pop)

            t0 = time.time()
            opt.optimize()
            elapsed = time.time() - t0

            print(
                f"\nPOET evolution completed in {elapsed:.2f}s ({elapsed/60:.2f} min)."
            )

            # CGP post-processing
            use_cgp = cfg.get("use_cgp", "False").strip().lower() == "true"
            if use_cgp:
                try:
                    from cgp_runner import run_cgp_on_model

                    cgp_t0 = time.time()
                    run_cgp_on_model(cfg)
                    cgp_el = time.time() - cgp_t0
                    print(f"CGP completed in {cgp_el:.2f}s ({cgp_el/60:.2f} min).")
                except ImportError as e:
                    print(f"WARNING: CGP skipped — {e}")

            total = time.time() - t0
            print(f"\nTotal: {total:.2f}s ({total/60:.2f} min).")

        except Exception as e:
            print(f"\n*** ERROR: {e}")
            import traceback

            traceback.print_exc()
        finally:
            sys.stdout, sys.stderr = old_stdout, old_stderr
            self._output_queue.put("\n__DONE__")

    # -------------------------------------------------------------------
    # Run: experiment mode
    # -------------------------------------------------------------------
    def _on_run_experiment(self):
        if self._run_thread and self._run_thread.is_alive():
            messagebox.showwarning("Busy", "A run is already in progress.")
            return

        config_paths = list(self._exp_listbox.get(0, tk.END))
        if len(config_paths) < 2:
            messagebox.showwarning(
                "Experiment", "Select at least 2 config files for comparison."
            )
            return

        try:
            replicates = int(self._exp_reps_var.get())
            gens = int(self._exp_gens_var.get())
            workers = int(self._exp_workers_var.get())
        except ValueError:
            messagebox.showerror(
                "Input Error", "Replicates/Gens/Workers must be integers."
            )
            return

        # Read all configs
        parser = eletility.ConfigParser()
        configs = []
        for p in config_paths:
            try:
                configs.append(parser.read(p))
            except Exception as e:
                messagebox.showerror("Config Error", f"Failed to read {p}:\n{e}")
                return

        self._stop_event.clear()
        self._clear_console()
        self._clear_chart()
        self._evo_path = ""  # experiment mode generates its own plots
        self._status_var.set("Running experiment...")

        self._run_thread = threading.Thread(
            target=self._run_experiment_thread,
            args=(configs, replicates, gens, workers),
            daemon=True,
        )
        self._run_thread.start()

    def _run_experiment_thread(
        self, configs: list[dict], replicates: int, gens: int, workers: int
    ):
        """Run experiment.run_experiment in a background thread."""
        import experiment

        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = _QueueWriter(self._output_queue)
        sys.stderr = _QueueWriter(self._output_queue)

        try:
            experiment.run_experiment(configs, replicates, gens, workers)
        except Exception as e:
            print(f"\n*** EXPERIMENT ERROR: {e}")
            import traceback

            traceback.print_exc()
        finally:
            sys.stdout, sys.stderr = old_stdout, old_stderr
            self._output_queue.put("\n__DONE__")

            # Show experiment plot if generated
            plot_path = "output/experiment/experiment_plot.png"
            if os.path.exists(plot_path):
                self._output_queue.put(f"\n__SHOW_PLOT__:{plot_path}")

    # -------------------------------------------------------------------
    # Stop
    # -------------------------------------------------------------------
    def _on_stop(self):
        if self._run_thread and self._run_thread.is_alive():
            self._stop_event.set()
            self._status_var.set("Stopping... (will finish current generation)")
            # We can't truly kill the thread, but we signal it
            messagebox.showinfo(
                "Stop",
                "Stop signal sent.\n\n"
                "The current generation will finish before stopping.\n"
                "If POET doesn't stop, close the window.",
            )
        else:
            self._status_var.set("Nothing to stop")

    # -------------------------------------------------------------------
    # Console output polling
    # -------------------------------------------------------------------
    def _poll_queue(self):
        """Drain the output queue and append to the console widget."""
        try:
            while True:
                text = self._output_queue.get_nowait()

                if text.strip() == "__DONE__":
                    self._status_var.set("Done")
                    self._update_chart_final()
                    continue

                if "__SHOW_PLOT__:" in text:
                    plot_path = text.split("__SHOW_PLOT__:", 1)[1].strip()
                    self._show_experiment_plot(plot_path)
                    continue

                self._console.configure(state=tk.NORMAL)
                self._console.insert(tk.END, text)
                self._console.see(tk.END)
                self._console.configure(state=tk.DISABLED)
        except queue.Empty:
            pass

        self.after(100, self._poll_queue)

    def _clear_console(self):
        self._console.configure(state=tk.NORMAL)
        self._console.delete("1.0", tk.END)
        self._console.configure(state=tk.DISABLED)

    # -------------------------------------------------------------------
    # Live chart
    # -------------------------------------------------------------------
    def _clear_chart(self):
        self._ax.clear()
        self._ax.set_xlabel("Generation", fontsize=10)
        self._ax.set_ylabel("Best Fitness", fontsize=10)
        self._ax.set_title("POET Evolution Progress", fontsize=11)
        self._ax.grid(True, alpha=0.3)
        self._canvas_chart.draw_idle()

    def _poll_chart(self):
        """Periodically re-read the evo CSV and update the chart."""
        if not self._evo_path or not os.path.exists(self._evo_path):
            if self._run_thread and self._run_thread.is_alive():
                self.after(500, self._poll_chart)
            return

        gens, best, avg = self._parse_evo_for_chart(self._evo_path)
        if gens:
            self._ax.clear()
            self._ax.plot(gens, best, "b-", linewidth=1.5, label="Best Fitness")
            self._ax.plot(
                gens, avg, "r--", linewidth=1.0, alpha=0.7, label="Avg Fitness"
            )
            self._ax.set_xlabel("Generation", fontsize=10)
            self._ax.set_ylabel("Fitness (lower = better)", fontsize=10)
            self._ax.set_title("POET Evolution Progress", fontsize=11)
            self._ax.legend(fontsize=9, loc="upper right")
            self._ax.grid(True, alpha=0.3)
            self._canvas_chart.draw_idle()

        # Keep polling while the thread is alive
        if self._run_thread and self._run_thread.is_alive():
            self.after(1000, self._poll_chart)

    def _update_chart_final(self):
        """One last chart update when the run finishes."""
        if self._evo_path and os.path.exists(self._evo_path):
            gens, best, avg = self._parse_evo_for_chart(self._evo_path)
            if gens:
                self._ax.clear()
                self._ax.plot(gens, best, "b-", linewidth=1.5, label="Best Fitness")
                self._ax.plot(
                    gens, avg, "r--", linewidth=1.0, alpha=0.7, label="Avg Fitness"
                )
                self._ax.set_xlabel("Generation", fontsize=10)
                self._ax.set_ylabel("Fitness (lower = better)", fontsize=10)
                self._ax.set_title(
                    f"POET Evolution — {len(gens)} generations", fontsize=11
                )
                self._ax.legend(fontsize=9, loc="upper right")
                self._ax.grid(True, alpha=0.3)
                self._canvas_chart.draw_idle()

    @staticmethod
    def _parse_evo_for_chart(
        path: str,
    ) -> tuple[list[int], list[float], list[float]]:
        """Parse the evo CSV into (gen_indices, best_fitness, avg_fitness)."""
        gens, best, avg = [], [], []
        try:
            with open(path, "r") as f:
                for line_no, line in enumerate(f):
                    if line_no == 0:
                        continue  # header
                    parts = line.strip().split(",")
                    if len(parts) >= 6:
                        try:
                            gens.append(int(parts[0].strip()))
                            best.append(float(parts[1].strip()))
                            avg.append(float(parts[5].strip()))
                        except ValueError:
                            pass
        except (IOError, OSError):
            pass
        return gens, best, avg

    def _show_experiment_plot(self, path: str):
        """Load and display an experiment comparison plot."""
        try:
            from PIL import Image

            img = Image.open(path)
            # Show in the chart area
            self._ax.clear()
            self._ax.imshow(img)
            self._ax.axis("off")
            self._ax.set_title("Experiment Comparison", fontsize=11)
            self._canvas_chart.draw_idle()
        except ImportError:
            # No PIL — just show the path
            self._ax.clear()
            self._ax.text(
                0.5,
                0.5,
                f"Plot saved to:\n{path}",
                ha="center",
                va="center",
                fontsize=12,
                transform=self._ax.transAxes,
            )
            self._ax.axis("off")
            self._canvas_chart.draw_idle()


# ═══════════════════════════════════════════════════════════════════════════
# Tooltip helper
# ═══════════════════════════════════════════════════════════════════════════
class _ToolTip:
    """Hover tooltip for a widget."""

    def __init__(self, widget: tk.Widget, text: str):
        self._widget = widget
        self._text = text
        self._tw: tk.Toplevel | None = None
        widget.bind("<Enter>", self._show)
        widget.bind("<Leave>", self._hide)

    def _show(self, event=None):
        if self._tw:
            return
        x = self._widget.winfo_rootx() + 20
        y = self._widget.winfo_rooty() + self._widget.winfo_height() + 2
        self._tw = tw = tk.Toplevel(self._widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        label = tk.Label(
            tw,
            text=self._text,
            background="#ffffcc",
            relief=tk.SOLID,
            borderwidth=1,
            font=("Segoe UI", 9),
            padx=6,
            pady=2,
        )
        label.pack()

    def _hide(self, event=None):
        if self._tw:
            self._tw.destroy()
            self._tw = None


def _create_tooltip(widget: tk.Widget, text: str):
    if text:
        _ToolTip(widget, text)


# ═══════════════════════════════════════════════════════════════════════════
# Entry point
# ═══════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    app = POETApp()
    app.mainloop()
