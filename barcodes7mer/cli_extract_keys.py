from __future__ import annotations

import argparse

from .barcode_reader import read_sample_list
from .extract_keys import ExtractKeysConfig, extract_keys_for_sample


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="barcodeAmplicon-extract-keys")
    ap.add_argument("sample_list_tsv", help="TSV: run<TAB>fastq(.gz)")
    ap.add_argument("--out-dir", default="keys_out", help="Output directory (creates one subdir per run)")
    ap.add_argument("--parts", type=int, default=128, help="Number of key partitions per run")
    ap.add_argument("--progress-every", "--progress_every", type=int, default=1_000_000)
    ap.add_argument("--log-file", "--log_file", default="ExtractKeys.progress.tsv")
    ap.add_argument("--strict-fastq", "--strict_fastq", action=argparse.BooleanOptionalAction, default=True)
    args = ap.parse_args(argv)

    samples = read_sample_list(args.sample_list_tsv)
    cfg = ExtractKeysConfig(
        out_dir=args.out_dir,
        parts=args.parts,
        progress_every=args.progress_every,
        strict_fastq=args.strict_fastq,
        log_file=args.log_file,
    )

    for s in samples:
        extract_keys_for_sample(s.run, s.r1_fastq_gz, cfg)

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
