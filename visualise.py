"""
visualise.py — Interactive Python visualiser for the Relativistic Black-Scholes pricer.

All computation stays in rel_bs.exe; this script calls it with --csv <mode>,
reads the CSV output, and renders matplotlib figures.

Usage:
    python visualise.py                    # render all figures
    python visualise.py iv_surface         # single panel
    python visualise.py density greeks     # two panels

Available panels:
    iv_surface   IV smile heatmap + smile curves by tau_rel
    density      Risk-neutral density: RBS vs Gaussian vs Cauchy
    greeks       Delta and Gamma surface heatmaps
    american     American put premium vs strike and tau_rel
    barrier      Barrier discount vs barrier level
    calibration  RMSE landscape: how well tau_rel is identified

Requires: matplotlib, numpy, pandas (pip install matplotlib numpy pandas)
"""

import subprocess, sys, io, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

# ── locate executable ─────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(SCRIPT_DIR, "rel_bs.exe")
if not os.path.exists(EXE):
    EXE = os.path.join(SCRIPT_DIR, "rel_bs")   # Linux/Mac fallback


def run_csv(mode: str) -> pd.DataFrame:
    """Call rel_bs --csv <mode> and return result as a DataFrame."""
    result = subprocess.run([EXE, "--csv", mode], capture_output=True, text=True, check=True)
    return pd.read_csv(io.StringIO(result.stdout))


# ── styling ───────────────────────────────────────────────────────────────────
PALETTE = plt.cm.plasma
plt.rcParams.update({
    "figure.facecolor": "#0f0f0f",
    "axes.facecolor":   "#1a1a1a",
    "axes.edgecolor":   "#444",
    "axes.labelcolor":  "#ccc",
    "xtick.color":      "#aaa",
    "ytick.color":      "#aaa",
    "text.color":       "#ddd",
    "grid.color":       "#333",
    "grid.linewidth":   0.5,
    "legend.facecolor": "#222",
    "legend.edgecolor": "#555",
})
ACCENT_COLORS = ["#ff6b6b", "#ffd93d", "#6bcb77", "#4d96ff", "#c77dff",
                 "#ff9a3c", "#00b4d8", "#e9c46a", "#f4a261", "#e76f51"]


# ═══════════════════════════════════════════════════════════════════════════════
#  Panel 1 — IV Surface
# ═══════════════════════════════════════════════════════════════════════════════
def plot_iv_surface(df: pd.DataFrame):
    fig = plt.figure(figsize=(14, 6), facecolor="#0f0f0f")
    fig.suptitle("RBS Implied Volatility Surface", fontsize=14, color="#eee", y=0.97)
    gs = GridSpec(1, 2, figure=fig, wspace=0.35)

    # Left: heatmap
    ax1 = fig.add_subplot(gs[0])
    pivot = df.pivot(index="K", columns="tau_rel", values="iv_pct")
    tau_labels = [f"{t:.0e}" for t in pivot.columns]
    im = ax1.imshow(
        pivot.values, aspect="auto", origin="lower",
        cmap="plasma", vmin=pivot.values.min(), vmax=pivot.values.max()
    )
    ax1.set_xticks(range(len(tau_labels)))
    ax1.set_xticklabels(tau_labels, rotation=45, ha="right", fontsize=7)
    ax1.set_yticks(range(len(pivot.index)))
    ax1.set_yticklabels([f"{k:.0f}" for k in pivot.index], fontsize=7)
    ax1.set_xlabel("τ_rel")
    ax1.set_ylabel("Strike K")
    ax1.set_title("IV (%) heatmap", fontsize=10, color="#bbb")
    cb = fig.colorbar(im, ax=ax1, shrink=0.8)
    cb.set_label("IV (%)", color="#ccc")
    cb.ax.yaxis.set_tick_params(color="#ccc")

    # Right: smile curves per tau_rel (select a subset)
    ax2 = fig.add_subplot(gs[1])
    tau_show = sorted(df["tau_rel"].unique())
    step = max(1, len(tau_show) // 6)
    for i, tr in enumerate(tau_show[::step]):
        sub = df[df["tau_rel"] == tr].sort_values("K")
        ax2.plot(sub["K"], sub["iv_pct"],
                 color=ACCENT_COLORS[i % len(ACCENT_COLORS)],
                 label=f"τ={tr:.0e}", linewidth=1.8)
    ax2.axhline(20.0, color="#555", linestyle="--", linewidth=0.8, label="BS flat 20%")
    ax2.set_xlabel("Strike K")
    ax2.set_ylabel("Implied Volatility (%)")
    ax2.set_title("IV Smile by τ_rel", fontsize=10, color="#bbb")
    ax2.legend(fontsize=7, loc="lower center", ncol=2)
    ax2.grid(True)

    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
#  Panel 2 — Risk-neutral Density
# ═══════════════════════════════════════════════════════════════════════════════
def plot_density(df: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), facecolor="#0f0f0f")
    fig.suptitle("Risk-Neutral Density: RBS vs Gaussian vs Cauchy  (τ_rel=0.10)", fontsize=13, color="#eee")

    inside  = df[df["inside_cone"] == 1]
    outside = df[df["inside_cone"] == 0]

    for ax, title, scale in zip(axes, ["Linear scale", "Log scale (tail detail)"], ["linear", "log"]):
        ax.fill_between(inside["x"],  inside["k_rbs"],    alpha=0.25, color="#4d96ff", label="_")
        ax.plot(df["x"], df["k_rbs"],    color="#4d96ff",  linewidth=2.0, label="RBS (Breeden-Litzenberger)")
        ax.plot(df["x"], df["k_gauss"],  color="#6bcb77",  linewidth=1.5, linestyle="--", label="Gaussian (BS)")
        ax.plot(df["x"], df["k_cauchy"], color="#ff6b6b",  linewidth=1.5, linestyle=":",  label="Cauchy (γ_C=σ)")

        # shade outside-cone region
        ax.fill_between(outside["x"], 0, outside["k_cauchy"],
                        alpha=0.12, color="#ff6b6b", label="_")

        # cone boundary lines
        cone_x = df[df["inside_cone"] == 0]["x"].min() if len(outside) else 0.0
        for xb in [-abs(cone_x), abs(cone_x)]:
            ax.axvline(xb, color="#ffd93d", linewidth=0.9, linestyle="--", alpha=0.7,
                       label="Light cone" if xb > 0 else "_")

        ax.set_xlabel("Log-return  x = ln(S_T / K)")
        ax.set_ylabel("Density")
        ax.set_title(title, fontsize=10, color="#bbb")
        ax.legend(fontsize=8)
        ax.set_yscale(scale)
        if scale == "log":
            ax.set_ylim(bottom=1e-5)
        ax.grid(True)

    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
#  Panel 3 — Greeks Surface
# ═══════════════════════════════════════════════════════════════════════════════
def plot_greeks(df: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), facecolor="#0f0f0f")
    fig.suptitle("RBS Greeks Surface  (Call, S=100, T=1yr, r=5%, σ=20%)", fontsize=13, color="#eee")

    for ax, col, title, cmap in zip(
        axes,
        ["delta", "gamma"],
        ["Delta  Δ = ∂C/∂S", "Gamma  Γ = ∂²C/∂S²"],
        ["plasma", "viridis"]
    ):
        pivot = df.pivot(index="K", columns="tau_rel", values=col)
        tau_labels = [f"{t:.0e}" for t in pivot.columns]
        im = ax.imshow(pivot.values, aspect="auto", origin="lower", cmap=cmap)
        ax.set_xticks(range(len(tau_labels)))
        ax.set_xticklabels(tau_labels, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels([f"{k:.0f}" for k in pivot.index], fontsize=7)
        ax.set_xlabel("τ_rel")
        ax.set_ylabel("Strike K")
        ax.set_title(title, fontsize=10, color="#bbb")
        cb = fig.colorbar(im, ax=ax, shrink=0.8)
        cb.ax.yaxis.set_tick_params(color="#ccc")

    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
#  Panel 4 — American Put Premium
# ═══════════════════════════════════════════════════════════════════════════════
def plot_american(df: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), facecolor="#0f0f0f")
    fig.suptitle("American Put Early-Exercise Premium  (RBS)", fontsize=13, color="#eee")

    tau_vals = sorted(df["tau_rel"].unique())

    # Left: premium in $ by strike, curve per tau_rel
    ax = axes[0]
    for i, tr in enumerate(tau_vals):
        sub = df[df["tau_rel"] == tr].sort_values("K")
        ax.plot(sub["K"], sub["premium"], color=ACCENT_COLORS[i], linewidth=2,
                label=f"τ={tr:.0e}")
    ax.set_xlabel("Strike K")
    ax.set_ylabel("Premium  (Am − Eu)")
    ax.set_title("Early-exercise premium ($)", fontsize=10, color="#bbb")
    ax.legend(fontsize=8)
    ax.grid(True)

    # Right: premium as % of European price (heatmap)
    ax2 = axes[1]
    pivot = df.pivot(index="K", columns="tau_rel", values="premium_pct")
    tau_labels = [f"{t:.0e}" for t in pivot.columns]
    im = ax2.imshow(pivot.values, aspect="auto", origin="lower", cmap="hot")
    ax2.set_xticks(range(len(tau_labels)))
    ax2.set_xticklabels(tau_labels, rotation=45, ha="right", fontsize=8)
    ax2.set_yticks(range(len(pivot.index)))
    ax2.set_yticklabels([f"{k:.0f}" for k in pivot.index], fontsize=8)
    ax2.set_xlabel("τ_rel")
    ax2.set_ylabel("Strike K")
    ax2.set_title("Premium as % of European put", fontsize=10, color="#bbb")
    cb = fig.colorbar(im, ax=ax2, shrink=0.8)
    cb.set_label("Premium (%)", color="#ccc")
    cb.ax.yaxis.set_tick_params(color="#ccc")

    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
#  Panel 5 — Barrier Discount
# ═══════════════════════════════════════════════════════════════════════════════
def plot_barrier(df: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), facecolor="#0f0f0f")
    fig.suptitle("Down-and-Out Barrier Call Discount  (ATM, S=K=100)", fontsize=13, color="#eee")

    tau_vals = sorted(df["tau_rel"].unique())

    ax = axes[0]
    for i, tr in enumerate(tau_vals):
        sub = df[df["tau_rel"] == tr].sort_values("H")
        ax.plot(sub["H"], sub["barrier_price"], color=ACCENT_COLORS[i],
                linewidth=2, label=f"τ={tr:.0e}")
    ax.axhline(df["vanilla"].iloc[0], color="#555", linestyle="--", linewidth=1, label="Vanilla")
    ax.set_xlabel("Barrier level H")
    ax.set_ylabel("Option Price")
    ax.set_title("Barrier call price vs H", fontsize=10, color="#bbb")
    ax.legend(fontsize=8)
    ax.grid(True)

    ax2 = axes[1]
    for i, tr in enumerate(tau_vals):
        sub = df[df["tau_rel"] == tr].sort_values("H")
        ax2.plot(sub["H"], sub["discount_pct"], color=ACCENT_COLORS[i],
                 linewidth=2, label=f"τ={tr:.0e}")
    ax2.set_xlabel("Barrier level H")
    ax2.set_ylabel("Discount from Vanilla (%)")
    ax2.set_title("Barrier discount % vs H", fontsize=10, color="#bbb")
    ax2.legend(fontsize=8)
    ax2.grid(True)

    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
#  Panel 6 — Calibration Landscape
# ═══════════════════════════════════════════════════════════════════════════════
def plot_calibration(df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(10, 5), facecolor="#0f0f0f")
    fig.suptitle("Calibration Landscape: RMSE vs τ_rel  (synthetic IV surface, τ_true=0.05)",
                 fontsize=13, color="#eee")

    ax.plot(df["tau_rel"], df["rmse_bp"], color="#4d96ff", linewidth=2)
    ax.axvline(0.05, color="#ffd93d", linestyle="--", linewidth=1.2, label="τ_true = 0.05")
    ax.fill_between(df["tau_rel"], df["rmse_bp"],
                    where=df["rmse_bp"] < 1.0, alpha=0.2, color="#6bcb77",
                    label="RMSE < 1 bp")

    ax.set_xscale("log")
    ax.set_xlabel("τ_rel  (log scale)")
    ax.set_ylabel("RMSE (basis points)")
    ax.set_title("Minimum at τ_true → golden-section search converges to 9.4e-09 error",
                 fontsize=9, color="#bbb")
    ax.legend(fontsize=9)
    ax.grid(True, which="both")

    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
#  Dispatch
# ═══════════════════════════════════════════════════════════════════════════════
PANELS = {
    "iv_surface":  (plot_iv_surface,  "iv_surface"),
    "density":     (plot_density,     "density"),
    "greeks":      (plot_greeks,      "greeks"),
    "american":    (plot_american,    "american"),
    "barrier":     (plot_barrier,     "barrier"),
    "calibration": (plot_calibration, "calibration"),
}

def main():
    requested = sys.argv[1:] if len(sys.argv) > 1 else list(PANELS.keys())
    unknown = [p for p in requested if p not in PANELS]
    if unknown:
        print(f"Unknown panels: {unknown}")
        print(f"Available: {list(PANELS.keys())}")
        sys.exit(1)

    figs = []
    for name in requested:
        plot_fn, csv_mode = PANELS[name]
        print(f"  computing {csv_mode} ...", end=" ", flush=True)
        df = run_csv(csv_mode)
        print(f"({len(df)} rows)")
        fig = plot_fn(df)
        fig.canvas.manager.set_window_title(f"RBS — {name}") if hasattr(fig.canvas, "manager") else None
        figs.append(fig)

    plt.show()


if __name__ == "__main__":
    main()
