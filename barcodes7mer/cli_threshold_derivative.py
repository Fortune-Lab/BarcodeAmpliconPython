from __future__ import annotations

import argparse
import glob
import os

from .threshold_core import build_chimera_percent_map_from_cwd, threshold_table_for_reads_file, write_threshold_tsv


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="barcodeAmplicon-threshold-derivative")
    ap.add_argument("--in", dest="inputs", action="append", default=[], help="repeatable reads TSV input")
    ap.add_argument("--prefix", default="", help="output prefix (used if --out not set)")
    ap.add_argument("--out", default="", help="explicit output path")

    ap.add_argument("--min-reads", "--min_reads", type=int, default=10000)
    ap.add_argument("--percent-cutoff", "--percent_cutoff", type=int, default=0)
    ap.add_argument("--add-chimeras", "--add_chimeras", choices=["TRUE", "FALSE"], default="FALSE")
    ap.add_argument("--filter-umi", "--filter_umi", choices=["TRUE", "FALSE"], default="FALSE")
    ap.add_argument("--run-mode", "--run_mode", dest="assay_mode", choices=["TRUE", "FALSE"], default="TRUE")

    ap.add_argument("--window", type=int, default=1)
    ap.add_argument("--shoulder", type=int, default=0)
    ap.add_argument("--search-start", "--search_start", type=int, default=0)
    ap.add_argument("--search-stop", "--search_stop", type=int, default=0)

    ap.add_argument("--python-compat", "--python_compat", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--min-index", "--min_index", type=int, default=2)
    ap.add_argument("--float-fmt", "--float_fmt", default="{:.12g}")
    args = ap.parse_args(argv)

    inputs = list(args.inputs)
    if not inputs:
        inputs = sorted(glob.glob("*_reads.tsv"))
        if not inputs:
            raise SystemExit("No input reads TSV provided via --in, and no *_reads.tsv found in current directory")

    for f in inputs:
        if not os.path.exists(f):
            raise SystemExit(f"Input file not found: {f}")

    if args.out:
        out_path = args.out
        prefix_for_oneoff = None
    else:
        if not args.prefix:
            raise SystemExit("Required: either --out <path> or --prefix <output_prefix>")
        prefix = os.path.splitext(args.prefix)[0]
        out_path = f"{prefix}window{args.window}_derivative.tsv"
        prefix_for_oneoff = prefix

    chimera_pct = build_chimera_percent_map_from_cwd()

    all_rows = []
    for reads_file in inputs:
        run, rows, failed = threshold_table_for_reads_file(
            reads_file,
            method="derivative",
            min_reads=args.min_reads,
            percent_cutoff=args.percent_cutoff,
            add_chimeras=(args.add_chimeras == "TRUE"),
            filter_umi=(args.filter_umi == "TRUE"),
            assay_mode=(args.assay_mode == "TRUE"),
            window=args.window,
            shoulder=args.shoulder,
            search_start=args.search_start,
            search_stop=args.search_stop,
            python_compat=args.python_compat,
            min_index_allowed=args.min_index,
            out_prefix_for_chimera_oneoff=prefix_for_oneoff,
            float_fmt=args.float_fmt,
            chimera_pct_map=chimera_pct,
        )
        if failed:
            import sys
            print(f"[WARN] {run}: FAILED: {failed}", file=sys.stderr)
            continue
        all_rows.extend(rows)

    write_threshold_tsv(out_path, all_rows, float_fmt=args.float_fmt)
    print(f"Wrote {out_path}")
    return 0
