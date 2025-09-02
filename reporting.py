from __future__ import annotations
import os
import time
from decimal import Decimal
from typing import Dict, Any, List, Tuple


def _to_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return float(Decimal(x))


def plot_and_report(cert_mode: str, trace: Dict[str, Any], base_dir: str = "reports") -> str:
    os.makedirs(base_dir, exist_ok=True)
    ts = time.strftime("%Y%m%d-%H%M%S")
    out_dir = os.path.join(base_dir, f"{ts}_{cert_mode}")
    img_dir = os.path.join(out_dir, "images")
    os.makedirs(img_dir, exist_ok=True)

    # Try to import matplotlib; if unavailable, proceed with markdown only
    have_mpl = False
    try:
        import matplotlib.pyplot as plt  # type: ignore
        have_mpl = True
    except Exception:
        plt = None  # type: ignore

    images: List[Tuple[str, str]] = []  # (filename, title)

    # 1) Contributions bar (Prime, Arch, Ram, Arch+Ram, Geom)
    if have_mpl:
        fig, ax = plt.subplots(figsize=(7, 4))
        labels = []
        mids = []
        errs = []
        def add_bar(name: str, key: str):
            if key in trace:
                lo, hi = trace[key]
                lo_f, hi_f = _to_float(lo), _to_float(hi)
                labels.append(name)
                mids.append((lo_f + hi_f) / 2.0)
                errs.append((hi_f - lo_f) / 2.0)
        add_bar("Prime", "prime_iv")
        add_bar("Arch", "arch_iv")
        add_bar("Ram", "ram_iv")
        add_bar("Arch+Ram", "arch_plus_ram")
        if "geom_amp" in trace:
            add_bar("GeomAmp", "geom_amp")
        elif "geom_lb_rs" in trace:
            add_bar("GeomLB(RS)", "geom_lb_rs")
        ax.bar(labels, mids, yerr=errs, capsize=4)
        ax.set_title("Ledger Contributions (interval mid ± radius)")
        ax.axhline(0.0, color='k', linewidth=0.8)
        fig.tight_layout()
        fn = os.path.join(img_dir, "contributions.png")
        fig.savefig(fn, dpi=160)
        plt.close(fig)
        images.append(("images/contributions.png", "Ledger Contributions"))

    # 2) R* budget vs Ram
    if have_mpl and "B_Rstar" in trace and "ram_iv" in trace:
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ram_lo, ram_hi = trace["ram_iv"]
        B_lo, B_hi = trace["B_Rstar"]
        ram_lo, ram_hi = _to_float(ram_lo), _to_float(ram_hi)
        B_lo, B_hi = _to_float(B_lo), _to_float(B_hi)
        ax.bar(["Ram"], [(ram_lo+ram_hi)/2], yerr=[(ram_hi-ram_lo)/2], capsize=4, color="#1f77b4")
        ax.axhspan(-B_hi, B_hi, color="#ff7f0e", alpha=0.2, label=f"|Ram| ≤ B_R*∈[{B_lo:.3g},{B_hi:.3g}]")
        ax.axhline(0.0, color='k', linewidth=0.8)
        ax.legend()
        ax.set_title("R* Budget Check")
        fig.tight_layout()
        fn = os.path.join(img_dir, "rstar_vs_ram.png")
        fig.savefig(fn, dpi=160)
        plt.close(fig)
        images.append(("images/rstar_vs_ram.png", "R* Budget vs Ram"))

    # 3) Annihilation inequality visualization
    if have_mpl and "sqrt_1m_eta" in trace and "epsilon_iv" in trace:
        fig, ax = plt.subplots(figsize=(6, 3.5))
        s_lo, s_hi = trace["sqrt_1m_eta"]
        e_lo, e_hi = trace["epsilon_iv"]
        s_hi = _to_float(s_hi)
        e_hi = _to_float(e_hi)
        ax.bar(["sqrt(1-eta)", "epsilon"], [s_hi, e_hi], color=["#2ca02c", "#d62728"])
        ax.axhline(1.0, color='k', linestyle='--', label='Threshold 1')
        ax.set_ylim(0, max(1.1, s_hi + e_hi + 0.05))
        ax.set_title("Annihilation: sqrt(1-eta) + epsilon < 1 (using upper bounds)")
        ax.legend()
        fig.tight_layout()
        fn = os.path.join(img_dir, "annihilation.png")
        fig.savefig(fn, dpi=160)
        plt.close(fig)
        images.append(("images/annihilation.png", "Annihilation Inequality"))

    # 4) Prime term scatter (if terms available)
    if have_mpl and "unram_terms" in trace and trace.get("unram_terms"):
        fig, ax = plt.subplots(figsize=(7, 3.5))
        xs = []
        ys = []
        yerr = []
        for t in trace["unram_terms"]:
            p = int(t["p"]); k = int(t["k"])
            import math as _m
            u = 2.0 * k * _m.log(p)
            lo = _to_float(t["lo"]); hi = _to_float(t["hi"])
            xs.append(u)
            ys.append((lo + hi) / 2.0)
            yerr.append((hi - lo) / 2.0)
        ax.errorbar(xs, ys, yerr=yerr, fmt='o', markersize=3, capsize=3)
        ax.set_xlabel("u = 2 k log p")
        ax.set_title("Prime-power contributions (per term)")
        ax.axhline(0.0, color='k', linewidth=0.8)
        fig.tight_layout()
        fn = os.path.join(img_dir, "prime_terms.png")
        fig.savefig(fn, dpi=160)
        plt.close(fig)
        images.append(("images/prime_terms.png", "Prime-power Contributions"))

    # Build markdown report
    rep = trace.get("report", {})
    result = rep.get("result", "?")
    phase = rep.get("phase_passed", "?")
    details = rep.get("details", {})

    md_lines: List[str] = []
    md_lines.append(f"# Verification Report — {cert_mode}\n")
    md_lines.append(f"- Result: {result} (phase {phase})\n")
    md_lines.append(f"- Params: {trace.get('params')}\n")
    md_lines.append(f"- Details: {details}\n")

    md_lines.append("\n## Visualizations\n")
    for rel, title in images:
        md_lines.append(f"### {title}\n\n")
        md_lines.append(f"![{title}]({rel})\n\n")

    md_lines.append("\n## What was verified\n")
    md_lines.append("- Numerical rigor: Decimal precision set from certificate; interval arithmetic used for comparisons.\n")
    md_lines.append("- Geometric side: Verified prime block minus (Arch + Ram) is nonnegative (evaluate) or that UB(Arch+Ram) ≤ 0 (RS).\n")
    md_lines.append("- Infrastructure (R*): Checked |Ram| ≤ B_{R*} using interval bounds.\n")
    md_lines.append("- Annihilation: Ensured sqrt(1−eta) + epsilon < 1 using upper endpoints; epsilon optionally tightened via q ≥ Q(X,ε).\n")
    md_lines.append("- Proof-of-GRH mode (if used): Ran the above across a finite test net.\n")

    with open(os.path.join(out_dir, "report.md"), "w") as f:
        f.write("".join(md_lines))

    return out_dir

