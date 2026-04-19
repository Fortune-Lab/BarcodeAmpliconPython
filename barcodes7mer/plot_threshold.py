from __future__ import annotations

import argparse
import math
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def auto_stride(n, target_points=200_000):
    if n <= target_points:
        return 1
    return int(math.ceil(n / target_points))


def detect_columns(df):
    cols = set(df.columns)

    # Harmonized output (slope/curvature)
    if "slope" in cols or "curvature" in cols:
        return {
            "mode": "harmonized",
            "slope_col": "slope" if "slope" in cols else None,
            "slope_label": "slope",
            "curv_col": "curvature" if "curvature" in cols else None,
            "curv_label": "curvature",
        }

    if "window_slope" in cols:
        info = {
            "mode": "window_slope",
            "slope_col": "window_slope",
            "slope_label": "window regression slope (window_slope)",
            "curv_col": None,
            "curv_label": None,
        }
        if "d_window_slope_dx" in cols:
            info["curv_col"] = "d_window_slope_dx"
            info["curv_label"] = "d/dx(window_slope) (d_window_slope_dx)"
        return info

    info = {"mode": "derivatives", "slope_col": None, "slope_label": None, "curv_col": None, "curv_label": None}
    if "dydx" in cols:
        info["slope_col"] = "dydx"
        info["slope_label"] = "d(norm)/dx (Derivative1) (dydx)"
    if "dy2dx2" in cols:
        info["curv_col"] = "dy2dx2"
        info["curv_label"] = "d²(norm)/dx² (Derivative2) (dy2dx2)"
    return info


def resolve_run_name(available_runs, query):
    q = (query or "").strip()
    if not q:
        raise ValueError("Empty --run value")

    runs = [r.strip() for r in available_runs]
    runs_unique = sorted(set(runs))

    if q in runs_unique:
        return q

    q_low = q.lower()

    exact_ci = [r for r in runs_unique if r.lower() == q_low]
    if len(exact_ci) == 1:
        return exact_ci[0]
    if len(exact_ci) > 1:
        raise ValueError(f"Ambiguous run name (case-insensitive exact matches): {exact_ci}")

    pref = [r for r in runs_unique if r.lower().startswith(q_low)]
    if len(pref) == 1:
        return pref[0]
    if len(pref) > 1:
        raise ValueError(f"Ambiguous run prefix '{q}' matched: {pref}")

    sub = [r for r in runs_unique if q_low in r.lower()]
    if len(sub) == 1:
        return sub[0]
    if len(sub) > 1:
        raise ValueError(f"Ambiguous run substring '{q}' matched: {sub}")

    raise ValueError(f"No run matched '{q}'. Available runs: {runs_unique}")


def plot_run(df_run, run, out_png, stride=1, xmax=None, overlay_curvature=False, left_axis="log10_counts"):
    df = df_run.copy()

    if "run" in df.columns:
        df["run"] = df["run"].astype(str).str.strip()

    df["index"] = pd.to_numeric(df["index"], errors="coerce")
    df["counts"] = pd.to_numeric(df["counts"], errors="coerce")
    if "norm" in df.columns:
        df["norm"] = pd.to_numeric(df["norm"], errors="coerce")

    df = df.dropna(subset=["index", "counts"]).sort_values("index")

    if xmax is not None:
        df = df[df["index"] <= xmax].copy()

    df["log10_counts"] = np.log10(df["counts"].where(df["counts"] > 0, np.nan))
    if left_axis == "norm":
        if "norm" not in df.columns:
            raise ValueError("left_axis=norm requested but 'norm' column not present.")
        y_left = df["norm"]
        y_left_label = "norm (log2(counts/sum))"
    else:
        y_left = df["log10_counts"]
        y_left_label = "log10(counts)"

    thr = df.loc[df["is_threshold"].astype(str).str.upper() == "TRUE", "index"]
    thr_x = int(thr.iloc[0]) if len(thr) else None

    info = detect_columns(df)
    slope_col = info["slope_col"]
    slope_label = info["slope_label"]
    curv_col = info["curv_col"]
    curv_label = info["curv_label"]

    if slope_col is not None and slope_col in df.columns:
        df[slope_col] = pd.to_numeric(df[slope_col], errors="coerce")
    if curv_col is not None and curv_col in df.columns:
        df[curv_col] = pd.to_numeric(df[curv_col], errors="coerce")

    dfp = df.iloc[::stride, :].copy() if stride and stride > 1 else df

    plt.rcParams["figure.figsize"] = (10, 6)
    fig, ax1 = plt.subplots()

    ax1.set_title(run)

    if thr_x is not None:
        ax1.axvline(x=thr_x, color="g", linestyle="--", linewidth=1.5, label=f"threshold={thr_x}")

    ax1.plot(dfp["index"], y_left.loc[dfp.index], color="red", linewidth=1, label=y_left_label)
    ax1.set_xlabel("barcode index (rank)")
    ax1.set_ylabel(y_left_label, color="red")
    ax1.tick_params(axis="y", labelcolor="red")

    ax2 = ax1.twinx()
    if slope_col is not None and slope_col in dfp.columns and dfp[slope_col].notna().any():
        ax2.plot(dfp["index"], dfp[slope_col], color="blue", linewidth=1, label=slope_label)
        ax2.set_ylabel(slope_label, color="blue")
        ax2.tick_params(axis="y", labelcolor="blue")
    else:
        ax2.set_ylabel("")

    ax3 = None
    if overlay_curvature and curv_col is not None and curv_col in dfp.columns and dfp[curv_col].notna().any():
        ax3 = ax1.twinx()
        ax3.spines["right"].set_position(("axes", 1.12))
        ax3.spines["right"].set_visible(True)
        ax3.plot(dfp["index"], dfp[curv_col], color="black", linewidth=1, alpha=0.7, label=curv_label)
        ax3.set_ylabel(curv_label, color="black")
        ax3.tick_params(axis="y", labelcolor="black")

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    if ax3 is not None:
        h3, l3 = ax3.get_legend_handles_labels()
        ax1.legend(h1 + h2 + h3, l1 + l2 + l3, loc="best", fontsize=9)
    else:
        ax1.legend(h1 + h2, l1 + l2, loc="best", fontsize=9)

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("threshold_tsv", help="threshold TSV")
    ap.add_argument("--outdir", default="threshold_plots", help="output directory")
    ap.add_argument("--run", default=None, help="run name (supports prefix/substring matching)")
    ap.add_argument("--xmax", type=int, default=None, help="max index to plot (optional)")
    ap.add_argument("--stride", type=int, default=0, help="plot every Nth point (0=auto)")
    ap.add_argument("--overlay_curvature", action="store_true",
                    help="overlay curvature metric on a third axis")
    ap.add_argument("--left_axis", choices=["log10_counts", "norm"], default="log10_counts",
                    help="what to plot on left axis (default log10_counts)")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.threshold_tsv, sep="\t", dtype=str, comment="#")

    required = {"run", "index", "counts", "is_threshold"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns in TSV: {sorted(missing)}")

    df["run"] = df["run"].astype(str).str.strip()

    if args.run is not None:
        available = sorted(df["run"].unique())
        resolved = resolve_run_name(available, args.run)
        df = df[df["run"] == resolved].copy()
        if df.empty:
            raise SystemExit(f"No rows found for resolved run '{resolved}'. Available runs: {available}")
        args.run = resolved

    for run, df_run in df.groupby("run", sort=True):
        n = len(df_run)
        stride = args.stride if args.stride and args.stride > 0 else auto_stride(n)
        out_png = os.path.join(
            args.outdir,
            f"{run}.threshold.{args.left_axis}"
            + (".curvature" if args.overlay_curvature else "")
            + ".png"
        )
        plot_run(
            df_run,
            run,
            out_png,
            stride=stride,
            xmax=args.xmax,
            overlay_curvature=args.overlay_curvature,
            left_axis=args.left_axis,
        )
        print("Wrote", out_png)

    return 0
