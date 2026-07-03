from __future__ import annotations

import argparse
import os

from .reduce_mcounts import reduce_keys_dir_to_mcounts


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="barcodeAmplicon-reduce-mcounts")
    ap.add_argument("--keys-dir", "--keys_dir", required=True, help="Directory containing keys_part_*.tsv")
    ap.add_argument("--run", default="", help="Run name (default: basename of keys-dir)")
    ap.add_argument("--out", default="", help="Output mcounts TSV (default: {run}_mcounts.tsv)")
    ap.add_argument("--tmp-dir", "--tmp_dir", default=os.environ.get("TMPDIR", "/tmp"))
    ap.add_argument("--sort-mem", "--sort_mem", default="8G", help="GNU sort -S value (e.g. 4G, 8G)")
    ap.add_argument("--sort-bin", "--sort_bin", default="sort")

    # NEW:
    ap.add_argument(
        "--part",
        type=int,
        default=None,
        help="Reduce only a single partition index (0..parts-1). Intended for SLURM arrays.",
    )

    args = ap.parse_args(argv)

    run = args.run or os.path.basename(os.path.abspath(args.keys_dir.rstrip("/")))
    out_path = args.out or f"{run}_mcounts.tsv"

    reduce_keys_dir_to_mcounts(
        args.keys_dir,
        run=run,
        out_path=out_path,
        tmp_dir=args.tmp_dir,
        sort_mem=args.sort_mem,
        sort_bin=args.sort_bin,
        only_part=args.part,
    )
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
