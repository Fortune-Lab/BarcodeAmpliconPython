from __future__ import annotations

import argparse

from .barcode_reader import ensure_progress_log_header, process_sample, read_sample_list


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="barcodeAmplicon-barcode-reader")
    ap.add_argument("sample_list_tsv")
    ap.add_argument("--progress-every", "--progress_every", type=int, default=1_000_000)
    ap.add_argument("--log-file", "--log_file", default="BarcodeReader.progress.tsv")
    ap.add_argument("--strict-fastq", "--strict_fastq", action=argparse.BooleanOptionalAction, default=True)
    args = ap.parse_args(argv)

    samples = read_sample_list(args.sample_list_tsv)

    ensure_progress_log_header(args.log_file)
    with open(args.log_file, "at", encoding="utf-8", newline="\n") as log_fh:
        for s in samples:
            process_sample(
                s.run,
                s.r1_fastq_gz,
                progress_every=args.progress_every,
                log_fh=log_fh,
                strict_fastq=args.strict_fastq,
            )
    return 0
