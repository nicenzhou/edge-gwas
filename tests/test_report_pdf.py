"""
Helper to write test run results and figures to a single PDF.
Used by test_all_functions.py and test_all_functions_quick.py.
"""
from __future__ import annotations

import os
import textwrap
from datetime import datetime
from typing import List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure


def write_results_pdf(
    output_path: str,
    title: str,
    log_lines: List[str],
    summary: str,
    plot_paths: Optional[List[str]] = None,
    function_status: Optional[List[tuple]] = None,
) -> None:
    """
    Write a PDF report with title, summary, function status, log output, and optional figures.

    Args:
        output_path: Path for the output PDF file.
        title: Report title (e.g. "EDGE-GWAS Full Function Test").
        log_lines: Lines of captured stdout/log.
        summary: Short summary string (e.g. "All passed" or "2 failures").
        plot_paths: Optional list of paths to PNG/JPEG images to include.
        function_status: Optional list of (function_name, status_str) e.g. ("load_plink_data", "OK").
    """
    plot_paths = plot_paths or []
    function_status = function_status or []
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

    with PdfPages(output_path) as pdf:
        # Page 1: Title and summary
        fig = Figure(figsize=(8.5, 11))
        ax = fig.add_subplot(111)
        ax.axis("off")
        y = 0.95
        ax.text(0.5, y, title, transform=ax.transAxes, fontsize=16, ha="center", fontweight="bold")
        y -= 0.08
        ax.text(0.5, y, datetime.now().strftime("%Y-%m-%d %H:%M"), transform=ax.transAxes, fontsize=10, ha="center")
        y -= 0.12
        ax.text(0.1, y, "Summary", transform=ax.transAxes, fontsize=12, fontweight="bold")
        y -= 0.06
        for line in summary.split("\n"):
            ax.text(0.1, y, line, transform=ax.transAxes, fontsize=10)
            y -= 0.04
        if function_status:
            y -= 0.05
            ax.text(0.1, y, "Functions exercised", transform=ax.transAxes, fontsize=12, fontweight="bold")
            y -= 0.05
            for name, status in function_status:
                ax.text(0.15, y, f"  {name}: {status}", transform=ax.transAxes, fontsize=9, family="monospace")
                y -= 0.035
                if y < 0.05:
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)
                    fig = Figure(figsize=(8.5, 11))
                    ax = fig.add_subplot(111)
                    ax.axis("off")
                    y = 0.95
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Log output (wrapped, multiple pages if needed)
        log_text = "\n".join(log_lines)
        width = 100
        wrapped = []
        for line in log_text.split("\n"):
            wrapped.extend(textwrap.wrap(line, width=width) if line.strip() else [""])
        chunk_size = 80
        for i in range(0, len(wrapped), chunk_size):
            chunk = "\n".join(wrapped[i : i + chunk_size])
            fig = Figure(figsize=(8.5, 11))
            ax = fig.add_subplot(111)
            ax.axis("off")
            ax.text(0.05, 0.98, chunk, transform=ax.transAxes, fontsize=7, verticalalignment="top", family="monospace")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # Plot pages
        for path in plot_paths:
            if not os.path.isfile(path):
                continue
            try:
                img = plt.imread(path)
                fig = Figure(figsize=(8.5, 11))
                ax = fig.add_subplot(111)
                ax.imshow(img)
                ax.axis("off")
                ax.set_title(os.path.basename(path), fontsize=10)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
            except Exception:
                pass

    print(f"\nResults PDF saved to: {output_path}")
