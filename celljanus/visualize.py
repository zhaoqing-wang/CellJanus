"""
Visualisation module for CellJanus.

Generates publication-quality plots from Kraken2/Bracken results:
  • Stacked bar chart of taxonomic composition
  • Horizontal bar chart of top-N species abundance
  • Pie / donut chart of community composition
  • Krona-style sunburst (interactive HTML via matplotlib)
  • QC summary dashboard
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
# Colour palette
# ---------------------------------------------------------------------------

PALETTE = sns.color_palette("husl", 24)


def _save(fig: plt.Figure, path: Path, dpi: int = 300) -> Path:
    """Save and close a figure."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    get_logger().info(f"Saved plot → {path}")
    return path


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
) -> Path:
    """Horizontal bar chart of top-N most abundant species."""
    if cfg is None:
        cfg = CellJanusConfig()

    data = df.nlargest(top_n, value_col).copy()
    data = data.sort_values(value_col)

    fig, ax = plt.subplots(figsize=(10, max(4, top_n * 0.35)))
    colors = sns.color_palette("viridis", len(data))
    ax.barh(data[name_col], data[value_col], color=colors)
    ax.set_xlabel("Estimated Read Count")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.ticklabel_format(axis="x", style="scientific", scilimits=(0, 3))

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

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
) -> Path:
    """Donut chart showing relative abundance of top species."""
    if cfg is None:
        cfg = CellJanusConfig()

    data = df.nlargest(top_n, value_col).copy()
    other = df[~df.index.isin(data.index)][value_col].sum()
    if other > 0:
        other_row = pd.DataFrame([{name_col: "Other", value_col: other}])
        data = pd.concat([data, other_row], ignore_index=True)

    fig, ax = plt.subplots(figsize=(9, 9))
    colors = sns.color_palette("husl", len(data))
    wedges, texts, autotexts = ax.pie(
        data[value_col],
        labels=data[name_col],
        autopct=lambda pct: f"{pct:.1f}%" if pct > 2 else "",
        colors=colors,
        pctdistance=0.82,
        startangle=140,
        wedgeprops=dict(width=0.45, edgecolor="white", linewidth=1.5),
    )
    for t in texts:
        t.set_fontsize(8)
    for t in autotexts:
        t.set_fontsize(7)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=20)

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
) -> Path:
    """
    Stacked bar chart comparing community composition across samples.

    Parameters
    ----------
    samples : dict
        {sample_name: bracken_df} mapping.
    """
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
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, ncol=1)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

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
) -> Path:
    """Single-sample heatmap of species abundance (log10 scale)."""
    if cfg is None:
        cfg = CellJanusConfig()

    data = df.nlargest(top_n, value_col).copy()
    data = data.set_index(name_col)[[value_col]]
    data.columns = ["Abundance"]
    data["log10"] = np.log10(data["Abundance"].clip(lower=1))

    fig, ax = plt.subplots(figsize=(4, max(4, top_n * 0.3)))
    sns.heatmap(
        data[["log10"]],
        annot=data[["Abundance"]].values,
        fmt=".0f",
        cmap="YlOrRd",
        cbar_kws={"label": "log10(reads)"},
        linewidths=0.5,
        ax=ax,
    )
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_ylabel("")
    ax.set_xlabel("")

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
) -> Path:
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
    if cfg is None:
        cfg = CellJanusConfig()

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("CellJanus Pipeline Summary", fontsize=16, fontweight="bold", y=1.02)

    # --- Panel 1: QC ---
    ax = axes[0]
    ax.axis("off")
    ax.set_title("Quality Control", fontsize=13, fontweight="bold")
    reads_before = qc_report.get("total_reads_before", 0)
    reads_after = qc_report.get("total_reads_after", 0)
    retained_pct = 100.0 * reads_after / reads_before if reads_before > 0 else 0.0
    qc_text = (
        f"Reads before QC:  {reads_before:>12,}\n"
        f"Reads after QC:   {reads_after:>12,}\n"
        f"Retained:         {retained_pct:>11.1f}%\n"
        f"Q30 before:       {qc_report.get('q30_rate_before', 0):>11.2%}\n"
        f"Q30 after:        {qc_report.get('q30_rate_after', 0):>11.2%}\n"
        f"Adapters trimmed: {qc_report.get('adapter_trimmed_reads', 0):>12,}"
    )
    ax.text(
        0.1,
        0.5,
        qc_text,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="center",
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="#e8f4fd", alpha=0.9),
    )

    # --- Panel 2: Alignment ---
    ax = axes[1]
    ax.axis("off")
    ax.set_title("Host Alignment (hg38)", fontsize=13, fontweight="bold")
    total = align_stats.get("total_reads", 0)
    mapped = total - align_stats.get("aligned_0_times", 0)
    unmapped = align_stats.get("aligned_0_times", 0)
    rate = align_stats.get("overall_alignment_rate", "N/A")
    align_text = (
        f"Total reads:      {total:>12,}\n"
        f"Host-mapped:      {mapped:>12,}\n"
        f"Unmapped (microbe):{unmapped:>11,}\n"
        f"Alignment rate:   {rate:>12}"
    )
    ax.text(
        0.1,
        0.5,
        align_text,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="center",
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="#e8fde8", alpha=0.9),
    )

    # --- Panel 3: Classification ---
    ax = axes[2]
    ax.axis("off")
    ax.set_title("Microbial Classification", fontsize=13, fontweight="bold")
    ax.text(
        0.1,
        0.5,
        classify_summary,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="center",
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="#fde8e8", alpha=0.9),
    )

    plt.tight_layout()
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
    """Generate the standard set of visualisation outputs."""
    if cfg is None:
        cfg = CellJanusConfig()

    output_dir = Path(output_dir) / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    fmt = cfg.plot_format

    paths = [
        plot_abundance_bar(bracken_df, output_dir / f"abundance_bar.{fmt}", cfg=cfg),
        plot_abundance_pie(bracken_df, output_dir / f"abundance_pie.{fmt}", cfg=cfg),
        plot_taxonomy_heatmap(bracken_df, output_dir / f"abundance_heatmap.{fmt}", cfg=cfg),
    ]
    return paths
