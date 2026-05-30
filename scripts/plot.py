import argparse
import glob
import os
import sys
from typing import List, Optional
import numpy as np
import pandas as pd
import seaborn as sns

#!/usr/bin/env python3
"""
plot.py - Read tabular data files and create plots and summary tables.

Usage examples:
    python plot.py data/*.csv --x time --y value1 value2 --plot line --out results --table summary.csv
    python plot.py results/*.tsv --y all --subplot --plot scatter --out figs --table summary.html

Features:
- Reads CSV/TSV/whitespace-delimited files (pandas sniffing).
- Auto-detects columns; accepts column names or 0-based indices.
- Plot types: line, scatter, bar.
- Overlay all files or create one subplot per file.
- Generates summary statistics (mean, median, std) and saves as CSV or HTML.
"""


import matplotlib.pyplot as plt

# Optional seaborn styling if available
try:
        sns.set(style="whitegrid")
except Exception:
        pass


def read_table(path: str) -> pd.DataFrame:
        """Read a table using pandas with delimiter sniffing."""
        # Try common delimiters first for speed
        for sep in [",", "\t", ";", None]:
                try:
                        df = pd.read_csv(path, sep=sep, engine="python", comment="#")
                        if df.shape[1] > 1:
                                return df
                except Exception:
                        continue
        # Fallback: read as whitespace-delimited
        return pd.read_csv(path, delim_whitespace=True, comment="#")


def select_columns(df: pd.DataFrame, tokens: Optional[List[str]]) -> List[str]:
        """Return a list of column names to plot given tokens (names or indices)."""
        if tokens is None or len(tokens) == 0:
                # default: all numeric columns except index-like
                return [c for c in df.columns if np.issubdtype(df[c].dtype, np.number)]
        if len(tokens) == 1 and tokens[0].lower() == "all":
                return list(df.columns)
        cols = []
        for t in tokens:
                if t.isdigit():
                        idx = int(t)
                        if 0 <= idx < len(df.columns):
                                cols.append(df.columns[idx])
                else:
                        if t in df.columns:
                                cols.append(t)
        return cols


def parse_args():
        p = argparse.ArgumentParser(description="Read files and create plots and summary tables.")
        p.add_argument("files", nargs="+", help="Input file paths or glob patterns.")
        p.add_argument("--x", help="Column name or 0-based index to use for x-axis (default: first column).")
        p.add_argument("--y", nargs="+", help="Column names or indices to plot (default: numeric columns).")
        p.add_argument("--plot", choices=["line", "scatter", "bar"], default="line", help="Plot type.")
        p.add_argument("--subplot", action="store_true", help="Create one subplot per file.")
        p.add_argument("--out", default="plot", help="Output file prefix for figures.")
        p.add_argument("--table", help="Write summary table to this path (CSV or HTML by extension).")
        p.add_argument("--show", action="store_true", help="Show interactive plot window.")
        p.add_argument("--dpi", type=int, default=150, help="Output figure DPI.")
        return p.parse_args()


def ensure_files(patterns: List[str]) -> List[str]:
        files = []
        for pat in patterns:
                expanded = glob.glob(pat)
                if expanded:
                        files.extend(expanded)
                elif os.path.exists(pat):
                        files.append(pat)
        files = sorted(set(files))
        if not files:
                print("No input files found for patterns:", patterns, file=sys.stderr)
                sys.exit(2)
        return files


def make_plot(dfs: List[pd.DataFrame], paths: List[str], xcol: Optional[str], ycols: List[str],
                            plot_type: str, subplot: bool, out_prefix: str, dpi: int):
        n = len(dfs)
        if subplot:
                cols = min(3, n)
                rows = (n + cols - 1) // cols
                fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3 * rows), squeeze=False)
                axes = axes.flatten()
                for i, (df, path) in enumerate(zip(dfs, paths)):
                        ax = axes[i]
                        _plot_one(df, ax, xcol, ycols, plot_type, label=os.path.basename(path))
                        ax.set_title(os.path.basename(path))
                # hide extra axes
                for j in range(n, len(axes)):
                        axes[j].set_visible(False)
                fig.tight_layout()
                out_path = f"{out_prefix}_subplots.png"
                fig.savefig(out_path, dpi=dpi)
                print("Saved:", out_path)
        else:
                fig, ax = plt.subplots(figsize=(8, 5))
                for df, path in zip(dfs, paths):
                        _plot_one(df, ax, xcol, ycols, plot_type, label=os.path.basename(path))
                ax.legend()
                fig.tight_layout()
                out_path = f"{out_prefix}.png"
                fig.savefig(out_path, dpi=dpi)
                print("Saved:", out_path)
        return out_path


def _plot_one(df: pd.DataFrame, ax, xcol: Optional[str], ycols: List[str], plot_type: str, label: str):
        if xcol is None:
                x = df.iloc[:, 0]
        else:
                x = df[xcol]
        for ycol in ycols:
                if ycol not in df.columns:
                        continue
                y = df[ycol]
                if plot_type == "line":
                        ax.plot(x, y, label=f"{label}:{ycol}")
                elif plot_type == "scatter":
                        ax.scatter(x, y, label=f"{label}:{ycol}", s=10)
                elif plot_type == "bar":
                        # bar plot uses index spacing; offset bars slightly if multiple series
                        ax.bar(x, y, label=f"{label}:{ycol}", alpha=0.7)
        ax.set_xlabel(x.name if hasattr(x, "name") else "x")
        ax.set_ylabel(", ".join(ycols))


def summary_table(dfs: List[pd.DataFrame], paths: List[str], ycols: List[str]) -> pd.DataFrame:
        rows = []
        for df, path in zip(dfs, paths):
                name = os.path.basename(path)
                for col in ycols:
                        if col not in df.columns:
                                continue
                        ser = pd.to_numeric(df[col], errors="coerce").dropna()
                        rows.append({
                                "file": name,
                                "column": col,
                                "count": int(ser.count()),
                                "mean": float(ser.mean()) if not ser.empty else np.nan,
                                "median": float(ser.median()) if not ser.empty else np.nan,
                                "std": float(ser.std()) if not ser.empty else np.nan,
                                "min": float(ser.min()) if not ser.empty else np.nan,
                                "max": float(ser.max()) if not ser.empty else np.nan,
                        })
        df_sum = pd.DataFrame(rows)
        if not df_sum.empty:
                # combined rows
                combined = df_sum.groupby("column").agg({
                        "count": "sum",
                        "mean": "mean",
                        "median": "mean",
                        "std": "mean",
                        "min": "min",
                        "max": "max"
                }).reset_index()
                combined.insert(0, "file", "combined")
                df_sum = pd.concat([df_sum, combined], ignore_index=True, sort=False)
        return df_sum


def save_table(df: pd.DataFrame, path: str):
        ext = os.path.splitext(path)[1].lower()
        if ext in [".html", ".htm"]:
                df.to_html(path, index=False)
        else:
                df.to_csv(path, index=False)
        print("Saved summary:", path)


def main():
        args = parse_args()
        files = ensure_files(args.files)
        dfs = []
        for f in files:
                try:
                        df = read_table(f)
                        dfs.append(df)
                except Exception as e:
                        print(f"Failed reading {f}: {e}", file=sys.stderr)
        if not dfs:
                print("No dataframes loaded.", file=sys.stderr)
                sys.exit(1)

        # Determine x column name if provided
        xcol = None
        if args.x:
                if args.x.isdigit():
                        idx = int(args.x)
                        if 0 <= idx < len(dfs[0].columns):
                                xcol = dfs[0].columns[idx]
                else:
                        xcol = args.x if args.x in dfs[0].columns else None

        # Determine y columns
        ycols = select_columns(dfs[0], args.y)

        if not ycols:
                print("No columns selected for plotting.", file=sys.stderr)
                sys.exit(1)

        out_fig = make_plot(dfs, files, xcol, ycols, args.plot, args.subplot, args.out, args.dpi)

        # Summary table
        if args.table:
                tbl = summary_table(dfs, files, ycols)
                save_table(tbl, args.table)
        else:
                # print to stdout brief summary
                tbl = summary_table(dfs, files, ycols)
                if not tbl.empty:
                        print(tbl.to_string(index=False))
                else:
                        print("No numeric data for summary.")

        if args.show:
                plt.show()


if __name__ == "__main__":
        main()