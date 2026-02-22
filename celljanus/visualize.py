"""
Visualisation module for CellJanus.

Generates publication-quality plots from Kraken2/Bracken results:
  • Horizontal bar chart of top-N species abundance
  • Pie / donut chart of community composition
  • Stacked bar chart (multi-sample comparison)
  • Taxonomic abundance heatmap
  • QC summary dashboard

Every plot is saved as **both PNG (raster) and PDF (vector)** by default.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")  # non-interactive backend

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from celljanus.config import CellJanusConfig
from celljanus.utils import get_logger


# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------

_STYLE_APPLIED = False


def _apply_style() -> None:
    """Apply publication-quality matplotlib defaults once."""
    global _STYLE_APPLIED
    if _STYLE_APPLIED:
        return
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 11,
            "axes.titlesize": 14,
            "axes.titleweight": "bold",
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
        }
    )
    _STYLE_APPLIED = True


PALETTE = sns.color_palette("husl", 24)


def _save(fig: plt.Figure, path: Path, dpi: int = 300) -> list[Path]:
    """Save figure as both PNG and PDF. Returns list of saved paths."""
    _apply_style()
    path.parent.mkdir(parents=True, exist_ok=True)
    log = get_logger()
    saved: list[Path] = []

    # Always save PNG
    png_path = path.with_suffix(".png")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    saved.append(png_path)
    log.info(f"Saved plot → {png_path}")

    # Always save PDF (vector)
    pdf_path = path.with_suffix(".pdf")
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight", facecolor="white")
    saved.append(pdf_path)
    log.info(f"Saved plot → {pdf_path}")

    plt.close(fig)
    return saved


# ---------------------------------------------------------------------------
# 1. Top-N Abundance Bar Chart
# ---------------------------------------------------------------------------


def plot_abundance_bar(
    df: pd.DataFrame,
    output_path: Path,
    *,
    top_n: int = 20,
    value_col: str = "bracken_estimated",
    name_col: str = "name",
    title: str = "Top Microbial Species by Read Count",
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """Horizontal bar chart of top-N most abundant species."""
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    data = df.nlargest(top_n, value_col).copy()
    data = data.sort_values(value_col)
    n = len(data)

    fig, ax = plt.subplots(figsize=(10, max(4, n * 0.4)))
    colors = sns.color_palette("viridis", n)
    bars = ax.barh(data[name_col], data[value_col], color=colors, edgecolor="white", linewidth=0.5)

    # Add value labels on bars
    max_val = data[value_col].max()
    for bar, val in zip(bars, data[value_col]):
        pct = val / data[value_col].sum() * 100
        label = f" {val:,.0f} ({pct:.1f}%)"
        ax.text(
            bar.get_width() + max_val * 0.01,
            bar.get_y() + bar.get_height() / 2,
            label,
            va="center",
            ha="left",
            fontsize=9,
            color="#333333",
        )

    ax.set_xlabel("Estimated Read Count")
    ax.set_title(title, pad=12)
    ax.set_xlim(0, max_val * 1.25)  # room for labels
    ax.ticklabel_format(axis="x", style="plain")
    fig.tight_layout()

    return _save(fig, output_path, dpi=cfg.dpi)


# ---------------------------------------------------------------------------
# 2. Relative Abundance Pie / Donut
# ---------------------------------------------------------------------------


def plot_abundance_pie(
    df: pd.DataFrame,
    output_path: Path,
    *,
    top_n: int = 12,
    value_col: str = "bracken_estimated",
    name_col: str = "name",
    title: str = "Microbial Community Composition",
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """Donut chart showing relative abundance of top species."""
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    data = df.nlargest(top_n, value_col).copy()
    other = df[~df.index.isin(data.index)][value_col].sum()
    if other > 0:
        other_row = pd.DataFrame([{name_col: "Other", value_col: other}])
        data = pd.concat([data, other_row], ignore_index=True)

    total = data[value_col].sum()
    n = len(data)
    colors = sns.color_palette("husl", n)

    fig, ax = plt.subplots(figsize=(10, 8))
    wedges, texts, autotexts = ax.pie(
        data[value_col],
        autopct=lambda pct: f"{pct:.1f}%" if pct > 2 else "",
        colors=colors,
        pctdistance=0.78,
        startangle=140,
        wedgeprops=dict(width=0.42, edgecolor="white", linewidth=2),
    )
    for t in autotexts:
        t.set_fontsize(9)
        t.set_color("#333333")

    # Legend instead of crowded labels
    legend_labels = [
        f"{row[name_col]}  ({row[value_col] / total * 100:.1f}%)" for _, row in data.iterrows()
    ]
    ax.legend(
        wedges,
        legend_labels,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        fontsize=9,
        frameon=False,
    )

    # Centre text
    ax.text(0, 0, f"{total:,.0f}\nreads", ha="center", va="center", fontsize=13, fontweight="bold")
    ax.set_title(title, pad=16)
    fig.tight_layout()

    return _save(fig, output_path, dpi=cfg.dpi)


# ---------------------------------------------------------------------------
# 3. Stacked Bar Chart (for multi-sample comparison)
# ---------------------------------------------------------------------------


def plot_stacked_bar(
    samples: dict[str, pd.DataFrame],
    output_path: Path,
    *,
    top_n: int = 15,
    value_col: str = "fraction",
    name_col: str = "name",
    title: str = "Taxonomic Composition Across Samples",
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """
    Stacked bar chart comparing community composition across samples.

    Parameters
    ----------
    samples : dict
        {sample_name: bracken_df} mapping.
    """
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    # Collect all top species across samples
    all_species: set[str] = set()
    for df in samples.values():
        all_species.update(df.nlargest(top_n, value_col)[name_col].tolist())

    species_list = sorted(all_species)
    sample_names = list(samples.keys())

    matrix = pd.DataFrame(0.0, index=sample_names, columns=species_list)
    for sname, df in samples.items():
        for _, row in df.iterrows():
            if row[name_col] in species_list:
                matrix.loc[sname, row[name_col]] = row[value_col]

    # Add "Other" column
    matrix["Other"] = 1.0 - matrix.sum(axis=1)
    matrix = matrix.clip(lower=0)

    fig, ax = plt.subplots(figsize=(12, max(4, len(sample_names) * 0.6)))
    colors = sns.color_palette("husl", len(matrix.columns))
    matrix.plot(kind="barh", stacked=True, ax=ax, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Relative Abundance")
    ax.set_title(title, pad=12)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, ncol=1, frameon=False)
    fig.tight_layout()

    return _save(fig, output_path, dpi=cfg.dpi)


# ---------------------------------------------------------------------------
# 4. Taxonomic Tree Heatmap
# ---------------------------------------------------------------------------


def plot_taxonomy_heatmap(
    df: pd.DataFrame,
    output_path: Path,
    *,
    top_n: int = 30,
    value_col: str = "bracken_estimated",
    name_col: str = "name",
    title: str = "Microbial Abundance Heatmap",
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """Single-sample heatmap of species abundance (log10 scale)."""
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    data = df.nlargest(top_n, value_col).copy()
    data = data.set_index(name_col)[[value_col]]
    data.columns = ["Abundance"]
    data["log10"] = np.log10(data["Abundance"].clip(lower=1))

    n = len(data)
    fig, ax = plt.subplots(figsize=(5, max(4, n * 0.35)))
    sns.heatmap(
        data[["log10"]],
        annot=data[["Abundance"]].values,
        fmt=".0f",
        cmap="YlOrRd",
        cbar_kws={"label": "log₁₀(reads)", "shrink": 0.6},
        linewidths=0.8,
        linecolor="white",
        ax=ax,
    )
    ax.set_title(title, pad=12)
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.tick_params(axis="y", rotation=0)
    fig.tight_layout()

    return _save(fig, output_path, dpi=cfg.dpi)


# ---------------------------------------------------------------------------
# 5. QC Summary Dashboard
# ---------------------------------------------------------------------------


def plot_qc_dashboard(
    qc_report: dict,
    align_stats: dict,
    classify_summary: str,
    output_path: Path,
    *,
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """
    Combined dashboard showing QC, alignment, and classification summaries.

    Parameters
    ----------
    qc_report : dict
        Keys: total_reads_before, total_reads_after, q30_rate_before, etc.
    align_stats : dict
        Keys: total_reads, aligned_0_times, overall_alignment_rate, etc.
    classify_summary : str
        Text summary from classify module.
    """
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    fig = plt.figure(figsize=(20, 7))
    fig.suptitle(
        "CellJanus Pipeline Summary",
        fontsize=18,
        fontweight="bold",
        y=0.98,
    )

    # Use gridspec for better control
    gs = fig.add_gridspec(1, 3, wspace=0.25, left=0.05, right=0.95, top=0.88, bottom=0.08)

    # --- Panel 1: QC ---
    ax = fig.add_subplot(gs[0])
    ax.axis("off")
    reads_before = qc_report.get("total_reads_before", 0)
    reads_after = qc_report.get("total_reads_after", 0)
    retained_pct = 100.0 * reads_after / reads_before if reads_before > 0 else 0.0
    qc_lines = [
        ("Reads before QC", f"{reads_before:>10,}"),
        ("Reads after QC", f"{reads_after:>10,}"),
        ("Retained", f"{retained_pct:>9.1f}%"),
        ("Q30 before", f"{qc_report.get('q30_rate_before', 0):>9.2%}"),
        ("Q30 after", f"{qc_report.get('q30_rate_after', 0):>9.2%}"),
        ("Adapters trimmed", f"{qc_report.get('adapter_trimmed_reads', 0):>10,}"),
    ]
    qc_text = "\n".join(f"{k:<18s} {v}" for k, v in qc_lines)
    ax.text(
        0.5,
        0.92,
        "Quality Control",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="center",
        va="top",
        color="#1a73e8",
    )
    ax.text(
        0.5,
        0.78,
        qc_text,
        transform=ax.transAxes,
        fontsize=11,
        va="top",
        ha="center",
        fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#e8f0fe", edgecolor="#a8c7fa", alpha=0.95),
    )

    # --- Panel 2: Alignment ---
    ax = fig.add_subplot(gs[1])
    ax.axis("off")
    total = align_stats.get("total_reads", 0)
    mapped = total - align_stats.get("aligned_0_times", 0)
    unmapped = align_stats.get("aligned_0_times", 0)
    rate = align_stats.get("overall_alignment_rate", "N/A")
    align_lines = [
        ("Total reads", f"{total:>10,}"),
        ("Host-mapped", f"{mapped:>10,}"),
        ("Unmapped (microbe)", f"{unmapped:>10,}"),
        ("Alignment rate", f"{rate:>10}"),
    ]
    align_text = "\n".join(f"{k:<20s} {v}" for k, v in align_lines)
    ax.text(
        0.5,
        0.92,
        "Host Alignment",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="center",
        va="top",
        color="#0d904f",
    )
    ax.text(
        0.5,
        0.78,
        align_text,
        transform=ax.transAxes,
        fontsize=11,
        va="top",
        ha="center",
        fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#e6f4ea", edgecolor="#a8dab5", alpha=0.95),
    )

    # --- Panel 3: Classification ---
    ax = fig.add_subplot(gs[2])
    ax.axis("off")
    ax.text(
        0.5,
        0.92,
        "Microbial Classification",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="center",
        va="top",
        color="#c5221f",
    )
    ax.text(
        0.5,
        0.78,
        classify_summary,
        transform=ax.transAxes,
        fontsize=10,
        va="top",
        ha="center",
        fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#fce8e6", edgecolor="#f5b7b1", alpha=0.95),
    )

    return _save(fig, output_path, dpi=cfg.dpi)


# ---------------------------------------------------------------------------
# Convenience: generate all default plots
# ---------------------------------------------------------------------------


def generate_all_plots(
    bracken_df: pd.DataFrame,
    output_dir: Path,
    *,
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """Generate the standard set of visualisation outputs (PNG + PDF each)."""
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    output_dir = Path(output_dir) / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: list[Path] = []
    paths.extend(plot_abundance_bar(bracken_df, output_dir / "abundance_bar.png", cfg=cfg))
    paths.extend(plot_abundance_pie(bracken_df, output_dir / "abundance_pie.png", cfg=cfg))
    paths.extend(plot_taxonomy_heatmap(bracken_df, output_dir / "abundance_heatmap.png", cfg=cfg))
    return paths
