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
            "font.size": 12,
            "axes.titlesize": 16,
            "axes.titleweight": "bold",
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 11,
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
            fontsize=11,
            color="#333333",
        )

    ax.set_xlabel("Estimated Read Count", fontsize=13)
    ax.set_title(title, fontsize=16, pad=12)
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

    fig, ax = plt.subplots(figsize=(12, 8))
    wedges, texts, autotexts = ax.pie(
        data[value_col],
        autopct=lambda pct: f"{pct:.1f}%" if pct > 2 else "",
        colors=colors,
        pctdistance=0.75,
        startangle=140,
        wedgeprops=dict(width=0.40, edgecolor="white", linewidth=2),
    )
    for t in autotexts:
        t.set_fontsize(13)
        t.set_fontweight("bold")
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
        fontsize=13,
        frameon=False,
    )

    # Centre text
    ax.text(0, 0, f"{total:,.0f}\nreads", ha="center", va="center", fontsize=16, fontweight="bold")
    ax.set_title(title, fontsize=16, pad=16)
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
        annot_kws={"fontsize": 12},
        cmap="YlOrRd",
        cbar_kws={"label": r"$\log_{10}$(reads)", "shrink": 0.6},
        linewidths=0.8,
        linecolor="white",
        ax=ax,
    )
    ax.set_title(title, fontsize=16, pad=12)
    ax.tick_params(axis="both", labelsize=12)
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

    fig = plt.figure(figsize=(18, 5.5))
    fig.suptitle(
        "CellJanus Pipeline Summary",
        fontsize=18,
        fontweight="bold",
        y=1.0,
    )

    # Tight gridspec — minimal vertical whitespace
    gs = fig.add_gridspec(
        1,
        3,
        wspace=0.15,
        left=0.02,
        right=0.98,
        top=0.85,
        bottom=0.05,
    )

    # Helper: draw a panel with title + key-value text filling the box
    def _draw_panel(gs_idx, title, color, bg, border, lines_text, font_sz=12):
        ax = fig.add_subplot(gs[gs_idx])
        ax.axis("off")
        # Title
        ax.text(
            0.5,
            0.97,
            title,
            transform=ax.transAxes,
            fontsize=15,
            fontweight="bold",
            ha="center",
            va="top",
            color=color,
        )
        # Content box — vertically centered in remaining space
        ax.text(
            0.5,
            0.48,
            lines_text,
            transform=ax.transAxes,
            fontsize=font_sz,
            va="center",
            ha="center",
            fontfamily="monospace",
            linespacing=1.6,
            bbox=dict(
                boxstyle="round,pad=0.8",
                facecolor=bg,
                edgecolor=border,
                alpha=0.95,
            ),
        )
        return ax

    # --- Panel 1: QC ---
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
    _draw_panel(0, "Quality Control", "#1a73e8", "#e8f0fe", "#a8c7fa", qc_text)

    # --- Panel 2: Alignment ---
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
    _draw_panel(1, "Host Alignment", "#0d904f", "#e6f4ea", "#a8dab5", align_text)

    # --- Panel 3: Classification ---
    _draw_panel(
        2, "Microbial Classification", "#c5221f", "#fce8e6", "#f5b7b1", classify_summary, font_sz=11
    )

    return _save(fig, output_path, dpi=cfg.dpi)


# ---------------------------------------------------------------------------
# Convenience: generate all default plots
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 6. scRNA-seq Cell × Species Heatmap
# ---------------------------------------------------------------------------


def plot_cell_species_heatmap(
    matrix_df: pd.DataFrame,
    output_path: Path,
    *,
    top_n_species: int = 20,
    top_n_cells: int = 100,
    max_display_cells: int = 200,
    sample_method: str = "stratified",
    title: str = "Per-Cell Microbial Abundance",
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """
    Heatmap showing microbial abundance per cell barcode.

    Optimized for large datasets (10,000+ cells) with smart sampling.

    Parameters
    ----------
    matrix_df : pd.DataFrame
        Cells × species matrix (cell barcodes as index).
    top_n_species : int
        Show only top N species by total abundance.
    top_n_cells : int
        Show only top N cells by total microbial reads. For large datasets,
        this is capped at max_display_cells.
    max_display_cells : int
        Maximum cells to display in heatmap (default: 200). Prevents memory issues.
    sample_method : str
        'top' = show top cells by read count
        'stratified' = sample across read count distribution
        'random' = random sample
    """
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    total_cells = len(matrix_df)
    total_species = len(matrix_df.columns)

    # Efficient species selection using vectorized operations
    species_totals = matrix_df.sum().sort_values(ascending=False)
    top_species = species_totals.head(top_n_species).index.tolist()

    # Smart cell sampling for large datasets
    cell_totals = matrix_df.sum(axis=1).sort_values(ascending=False)
    n_display = min(top_n_cells, max_display_cells, total_cells)

    if sample_method == "stratified" and total_cells > n_display * 2:
        # Stratified sampling: include top, middle, and bottom read-count cells
        n_top = n_display // 2
        n_middle = n_display // 4
        n_low = n_display - n_top - n_middle
        idx_top = cell_totals.head(n_top).index.tolist()
        idx_middle = cell_totals.iloc[total_cells // 4 : total_cells // 4 + n_middle].index.tolist()
        idx_low = cell_totals.iloc[-n_low:].index.tolist()
        sample_cells = idx_top + idx_middle + idx_low
    elif sample_method == "random" and total_cells > n_display:
        import random

        sample_cells = random.sample(list(cell_totals.index), n_display)
    else:
        sample_cells = cell_totals.head(n_display).index.tolist()

    plot_data = matrix_df.loc[sample_cells, top_species].copy()

    # Log transform for better visualization
    plot_data_log = np.log10(plot_data + 1)

    n_cells = len(plot_data)
    n_species = len(top_species)

    # Dynamic figure sizing with max limits
    fig_width = min(20, max(8, n_species * 0.5))
    fig_height = min(15, max(6, n_cells * 0.12))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Optimized heatmap: disable gridlines for large matrices
    sns.heatmap(
        plot_data_log,
        cmap="YlOrRd",
        cbar_kws={"label": r"$\log_{10}$(count + 1)", "shrink": 0.6},
        linewidths=0 if n_cells > 50 else 0.1,
        linecolor="white",
        ax=ax,
        xticklabels=True,
        yticklabels=n_cells <= 50,
    )

    sampling_note = (
        "" if total_cells <= n_display else f" (showing {n_display}/{total_cells} cells)"
    )
    ax.set_title(f"{title}{sampling_note}", fontsize=16, pad=12)
    ax.set_xlabel(f"Species (top {n_species} of {total_species})", fontsize=13)
    ax.set_ylabel("Cell Barcode" if n_cells <= 50 else f"Cells (n={n_cells})", fontsize=13)
    ax.tick_params(axis="x", rotation=45, labelsize=10)
    ax.tick_params(axis="y", rotation=0, labelsize=8)

    fig.tight_layout()
    return _save(fig, output_path, dpi=cfg.dpi)


def plot_cell_microbe_summary(
    matrix_df: pd.DataFrame,
    output_path: Path,
    *,
    title: str = "scRNA-seq Microbial Detection Summary",
    top_n_species: int = 15,
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """
    Summary plots for scRNA-seq microbial detection:
    - Distribution of microbial reads per cell
    - Species richness per cell
    - Top species across all cells

    Optimized for datasets with 10,000+ cells.
    """
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    total_cells = len(matrix_df)
    total_species = len(matrix_df.columns)

    fig = plt.figure(figsize=(16, 5))
    gs = fig.add_gridspec(1, 3, wspace=0.3)

    # Panel 1: Reads per cell distribution
    ax1 = fig.add_subplot(gs[0])
    reads_per_cell = matrix_df.sum(axis=1)

    # Use log scale bins for large dynamic range
    if reads_per_cell.max() > reads_per_cell.median() * 100:
        bins = np.logspace(
            np.log10(max(1, reads_per_cell.min())), np.log10(reads_per_cell.max() + 1), 51
        )
        ax1.hist(reads_per_cell + 1, bins=bins, color="#1a73e8", edgecolor="white", alpha=0.8)
        ax1.set_xscale("log")
        ax1.set_xlabel(r"Microbial Reads per Cell ($\log_{10}$)", fontsize=12)
    else:
        ax1.hist(reads_per_cell, bins=50, color="#1a73e8", edgecolor="white", alpha=0.8)
        ax1.set_xlabel("Microbial Reads per Cell", fontsize=12)

    ax1.set_ylabel("Number of Cells", fontsize=12)
    ax1.set_title(f"Read Depth (n={total_cells:,} cells)", fontsize=14, fontweight="bold")
    ax1.axvline(
        reads_per_cell.median(),
        color="red",
        linestyle="--",
        label=f"Median: {reads_per_cell.median():.0f}",
    )
    ax1.legend(fontsize=10)

    # Panel 2: Species richness per cell
    ax2 = fig.add_subplot(gs[1])
    richness = (matrix_df > 0).sum(axis=1)
    max_richness = int(richness.max())
    bins = range(0, min(max_richness + 2, 101))  # Cap at 100 bins
    ax2.hist(richness.clip(upper=100), bins=bins, color="#0d904f", edgecolor="white", alpha=0.8)
    ax2.set_xlabel("Number of Species Detected", fontsize=12)
    ax2.set_ylabel("Number of Cells", fontsize=12)
    ax2.set_title(f"Species Richness ({total_species:,} total)", fontsize=14, fontweight="bold")

    # Panel 3: Top species bar chart
    ax3 = fig.add_subplot(gs[2])
    species_totals = matrix_df.sum().sort_values(ascending=True).tail(top_n_species)
    bars = ax3.barh(
        range(len(species_totals)), species_totals.values, color="#c5221f", edgecolor="white"
    )
    ax3.set_yticks(range(len(species_totals)))
    ax3.set_yticklabels(species_totals.index, fontsize=10)
    ax3.set_xlabel("Total Reads Across All Cells", fontsize=12)
    ax3.set_title(f"Top {top_n_species} Species", fontsize=14, fontweight="bold")

    # Use log scale for x-axis if range is large
    if species_totals.max() > species_totals.min() * 100:
        ax3.set_xscale("log")

    fig.suptitle(title, fontsize=16, fontweight="bold", y=1.02)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    return _save(fig, output_path, dpi=cfg.dpi)


def plot_cell_bacteria_dotplot(
    matrix_df: pd.DataFrame,
    output_path: Path,
    *,
    top_n_species: int = 15,
    top_n_cells: int = 50,
    max_display_cells: int = 100,
    min_reads_per_cell: int = 1,
    title: str = "Cell–Bacteria Association",
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """
    Dot plot showing cell-bacteria associations.

    Dot size = number of reads, color = fraction of cell's microbial reads.

    Optimized for large datasets with smart filtering and sampling.

    Parameters
    ----------
    matrix_df : pd.DataFrame
        Cells × species matrix.
    top_n_species : int
        Number of top species to display.
    top_n_cells : int
        Number of top cells to display.
    max_display_cells : int
        Maximum cells to show (prevents memory issues).
    min_reads_per_cell : int
        Minimum reads to include a cell.
    """
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    total_cells = len(matrix_df)
    total_species = len(matrix_df.columns)

    # Pre-filter cells with minimum reads
    cell_totals = matrix_df.sum(axis=1)
    valid_cells = cell_totals[cell_totals >= min_reads_per_cell]

    # Filter to top species and cells
    species_totals = matrix_df.sum().sort_values(ascending=False)
    top_species = species_totals.head(top_n_species).index.tolist()

    n_display = min(top_n_cells, max_display_cells, len(valid_cells))
    top_cells = valid_cells.sort_values(ascending=False).head(n_display).index.tolist()

    plot_data = matrix_df.loc[top_cells, top_species].copy()

    # Calculate fraction for each cell
    row_sums = plot_data.sum(axis=1)
    fraction_data = plot_data.div(row_sums, axis=0).fillna(0)

    # Vectorized dot plot data generation for large datasets
    x_coords = []
    y_coords = []
    sizes = []
    colors = []

    for i, cell in enumerate(plot_data.index):
        for j, species in enumerate(plot_data.columns):
            count = plot_data.loc[cell, species]
            frac = fraction_data.loc[cell, species]
            if count > 0:
                x_coords.append(j)
                y_coords.append(i)
                sizes.append(np.sqrt(count) * 10 + 5)
                colors.append(frac)

    n_cells = len(top_cells)
    n_species_display = len(top_species)

    # Dynamic figure sizing
    fig_width = min(18, max(8, n_species_display * 0.6))
    fig_height = min(12, max(6, n_cells * 0.12))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    scatter = ax.scatter(
        x_coords,
        y_coords,
        s=sizes,
        c=colors,
        cmap="YlOrRd",
        alpha=0.8,
        edgecolors="white",
        linewidth=0.5,
    )

    ax.set_xticks(range(n_species_display))
    ax.set_xticklabels(top_species, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(n_cells))
    ax.set_yticklabels(top_cells if n_cells <= 30 else [""] * n_cells, fontsize=8)

    ax.set_xlabel(f"Species (top {n_species_display} of {total_species})", fontsize=13)
    ylabel = f"Cells (top {n_cells} of {total_cells})" if n_cells != total_cells else "Cell Barcode"
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(title, fontsize=16, pad=12)

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6)
    cbar.set_label("Fraction of Cell's Reads", fontsize=11)

    fig.tight_layout()
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


def plot_scrnaseq_dashboard(
    matrix_df: pd.DataFrame,
    output_path: Path,
    *,
    species_summary_df: Optional[pd.DataFrame] = None,
    cell_summary_df: Optional[pd.DataFrame] = None,
    pipeline_summary: Optional[dict] = None,
    title: str = "scRNA-seq Microbial Detection Dashboard",
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """
    Dashboard summarising scRNA-seq microbial detection results.

    Four panels:
    1. Top species abundance bar chart
    2. Per-cell read depth distribution
    3. Species prevalence (fraction of cells)
    4. Key metrics text panel
    """
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    total_cells = len(matrix_df)
    total_species = len(matrix_df.columns)

    # ── Panel 1: Top species bar chart ──
    ax1 = fig.add_subplot(gs[0, 0])
    species_totals = matrix_df.sum().sort_values(ascending=True)
    top_sp = species_totals.tail(min(15, total_species))
    colors_bar = plt.cm.YlOrRd(np.linspace(0.3, 0.9, len(top_sp)))
    ax1.barh(range(len(top_sp)), top_sp.values, color=colors_bar, edgecolor="white")
    ax1.set_yticks(range(len(top_sp)))
    ax1.set_yticklabels([s[:35] for s in top_sp.index], fontsize=9)
    ax1.set_xlabel("Total Reads", fontsize=11)
    ax1.set_title("Species Abundance", fontsize=14, fontweight="bold")
    for i, v in enumerate(top_sp.values):
        ax1.text(v + max(top_sp.values) * 0.01, i, f"{int(v)}", va="center", fontsize=9)

    # ── Panel 2: Reads per cell distribution ──
    ax2 = fig.add_subplot(gs[0, 1])
    reads_per_cell = matrix_df.sum(axis=1)
    if len(reads_per_cell) > 0:
        n_bins = min(50, max(10, int(reads_per_cell.max() - reads_per_cell.min() + 1)))
        ax2.hist(reads_per_cell, bins=n_bins, color="#1a73e8", edgecolor="white", alpha=0.85)
        ax2.axvline(
            reads_per_cell.median(),
            color="red",
            linestyle="--",
            linewidth=1.5,
            label=f"Median: {reads_per_cell.median():.0f}",
        )
        ax2.axvline(
            reads_per_cell.mean(),
            color="orange",
            linestyle=":",
            linewidth=1.5,
            label=f"Mean: {reads_per_cell.mean():.1f}",
        )
        ax2.legend(fontsize=9)
    ax2.set_xlabel("Microbial Reads per Cell", fontsize=11)
    ax2.set_ylabel("Number of Cells", fontsize=11)
    ax2.set_title("Per-Cell Read Depth", fontsize=14, fontweight="bold")

    # ── Panel 3: Species prevalence (% cells with each species) ──
    ax3 = fig.add_subplot(gs[1, 0])
    prevalence = (matrix_df > 0).sum() / total_cells * 100
    prevalence = prevalence.sort_values(ascending=True).tail(min(15, total_species))
    colors_prev = plt.cm.Greens(np.linspace(0.3, 0.9, len(prevalence)))
    ax3.barh(range(len(prevalence)), prevalence.values, color=colors_prev, edgecolor="white")
    ax3.set_yticks(range(len(prevalence)))
    ax3.set_yticklabels([s[:35] for s in prevalence.index], fontsize=9)
    ax3.set_xlabel("Cell Prevalence (%)", fontsize=11)
    ax3.set_title("Species Prevalence", fontsize=14, fontweight="bold")
    ax3.set_xlim(0, 105)
    for i, v in enumerate(prevalence.values):
        ax3.text(v + 1, i, f"{v:.0f}%", va="center", fontsize=9)

    # ── Panel 4: Key metrics ──
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis("off")

    if pipeline_summary is None:
        pipeline_summary = {
            "total_cells": total_cells,
            "species_detected": total_species,
            "total_microbial_reads": int(matrix_df.sum().sum()),
            "mean_reads_per_cell": float(reads_per_cell.mean()) if len(reads_per_cell) > 0 else 0,
        }

    richness = (matrix_df > 0).sum(axis=1)

    metrics = [
        ("Total Cells", f"{total_cells:,}"),
        ("Species Detected", f"{total_species}"),
        (
            "Total Microbial Reads",
            f"{pipeline_summary.get('total_microbial_reads', int(matrix_df.sum().sum())):,}",
        ),
        ("Mean Reads / Cell", f"{pipeline_summary.get('mean_reads_per_cell', 0):.1f}"),
        (
            "Median Reads / Cell",
            f"{reads_per_cell.median():.0f}" if len(reads_per_cell) > 0 else "0",
        ),
        ("Mean Species / Cell", f"{richness.mean():.1f}" if len(richness) > 0 else "0"),
        ("Max Species / Cell", f"{richness.max()}" if len(richness) > 0 else "0"),
    ]

    y_pos = 0.92
    ax4.text(
        0.5,
        y_pos + 0.06,
        "Key Metrics",
        transform=ax4.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="center",
        va="top",
    )
    for label, value in metrics:
        ax4.text(
            0.15, y_pos, label, transform=ax4.transAxes, fontsize=12, va="top", fontweight="medium"
        )
        ax4.text(
            0.85,
            y_pos,
            value,
            transform=ax4.transAxes,
            fontsize=12,
            va="top",
            ha="right",
            color="#1a73e8",
            fontweight="bold",
        )
        y_pos -= 0.12

    fig.suptitle(title, fontsize=18, fontweight="bold", y=1.01)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    return _save(fig, output_path, dpi=cfg.dpi)


def plot_scrnaseq_abundance_pie(
    matrix_df: pd.DataFrame,
    output_path: Path,
    *,
    top_n: int = 10,
    title: str = "Microbial Community Composition (All Cells)",
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """Donut chart of total microbial composition across all cells."""
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    species_totals = matrix_df.sum().sort_values(ascending=False)
    top = species_totals.head(top_n)
    if len(species_totals) > top_n:
        other = species_totals.iloc[top_n:].sum()
        top = pd.concat([top, pd.Series({"Other": other})])

    colors = plt.cm.Set3(np.linspace(0, 1, len(top)))
    fig, ax = plt.subplots(figsize=(8, 8))
    wedges, texts, autotexts = ax.pie(
        top.values,
        labels=None,
        colors=colors,
        autopct=lambda p: f"{p:.1f}%" if p > 3 else "",
        startangle=140,
        pctdistance=0.82,
        wedgeprops=dict(width=0.4, edgecolor="white", linewidth=2),
    )

    for t in autotexts:
        t.set_fontsize(11)
        t.set_fontweight("bold")

    ax.legend(
        wedges,
        [f"{s[:30]} ({int(v)})" for s, v in zip(top.index, top.values)],
        title="Species (reads)",
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        fontsize=10,
    )

    total = int(species_totals.sum())
    ax.text(0, 0, f"{total:,}\nreads", ha="center", va="center", fontsize=16, fontweight="bold")
    ax.set_title(title, fontsize=16, pad=20)

    fig.tight_layout()
    return _save(fig, output_path, dpi=cfg.dpi)


def generate_scrnaseq_plots(
    matrix_df: pd.DataFrame,
    output_dir: Path,
    *,
    species_summary_df: Optional[pd.DataFrame] = None,
    cell_summary_df: Optional[pd.DataFrame] = None,
    pipeline_summary: Optional[dict] = None,
    cfg: Optional[CellJanusConfig] = None,
) -> list[Path]:
    """Generate scRNA-seq specific visualisation outputs."""
    _apply_style()
    if cfg is None:
        cfg = CellJanusConfig()

    output_dir = Path(output_dir) / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: list[Path] = []

    if len(matrix_df) > 0 and len(matrix_df.columns) > 0:
        # Dashboard (like bulk pipeline_dashboard)
        paths.extend(
            plot_scrnaseq_dashboard(
                matrix_df,
                output_dir / "scrnaseq_dashboard.png",
                species_summary_df=species_summary_df,
                cell_summary_df=cell_summary_df,
                pipeline_summary=pipeline_summary,
                cfg=cfg,
            )
        )
        # Donut chart
        paths.extend(
            plot_scrnaseq_abundance_pie(
                matrix_df,
                output_dir / "scrnaseq_abundance_pie.png",
                cfg=cfg,
            )
        )
        # Summary panels (3-panel: read depth, richness, top species)
        paths.extend(
            plot_cell_microbe_summary(matrix_df, output_dir / "cell_microbe_summary.png", cfg=cfg)
        )
        # Note: Heatmap and dotplot removed - not suitable for real scRNA-seq (5000+ cells)

    return paths
