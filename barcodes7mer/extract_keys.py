from __future__ import annotations

import os
import zlib
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, TextIO, Tuple

from .barcode_reader import FULL_RE, QTAGS
from .fastq import iter_fastq_sequences


def _timestamp_local() -> str:
    import datetime
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def ensure_progress_log_header(path: str) -> None:
    need_header = (not os.path.exists(path)) or (os.path.getsize(path) == 0)
    if not need_header:
        return
    with open(path, "at", encoding="utf-8", newline="\n") as fh:
        fh.write(
            "\t".join(
                [
                    "timestamp",
                    "run",
                    "event",
                    "reads_processed",
                    "matched",
                    "percent_matched",
                    "reads_per_sec",
                    "elapsed_sec",
                    "input_files",
                ]
            )
            + "\n"
        )


def _neighbors_within_hamming2(seq: str) -> Iterable[str]:
    """
    Generate all sequences within Hamming distance <=2 of seq (including seq itself).
    seq is assumed length 12 (qtag).
    """
    bases = "ACGT"
    L = len(seq)

    yield seq  # distance 0

    # distance 1
    for i in range(L):
        orig = seq[i]
        for b in bases:
            if b == orig:
                continue
            yield seq[:i] + b + seq[i + 1 :]

    # distance 2
    for i in range(L - 1):
        orig_i = seq[i]
        for j in range(i + 1, L):
            orig_j = seq[j]
            for bi in bases:
                if bi == orig_i:
                    continue
                for bj in bases:
                    if bj == orig_j:
                        continue
                    # replace i and j
                    yield seq[:i] + bi + seq[i + 1 : j] + bj + seq[j + 1 :]


def build_qtag_lookup(qtags: Dict[str, str]) -> Dict[str, Optional[str]]:
    """
    Returns lookup[obs_12mer] -> qtag_id if uniquely assignable within <=2 mismatches,
    else None (ambiguous) or missing (no hit).
    """
    lut: Dict[str, Optional[str]] = {}
    for canon_seq, qtag_id in qtags.items():
        for neigh in _neighbors_within_hamming2(canon_seq):
            if neigh in lut:
                # if conflict, mark ambiguous
                if lut[neigh] != qtag_id:
                    lut[neigh] = None
            else:
                lut[neigh] = qtag_id
    return lut


def _stable_partition_index(qtag_id: str, barcode: str, parts: int) -> int:
    key = (qtag_id + "\t" + barcode).encode("utf-8")
    return zlib.crc32(key) % parts


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


@dataclass(frozen=True)
class ExtractKeysConfig:
    out_dir: str
    parts: int = 128
    progress_every: int = 1_000_000
    strict_fastq: bool = True
    log_file: str = "ExtractKeys.progress.tsv"


def extract_keys_for_sample(run: str, fastq_path: str, cfg: ExtractKeysConfig) -> Tuple[int, int]:
    """
    Writes partitioned key files:

        {cfg.out_dir}/{run}/keys_part_XXX.tsv

    Each line:
        qtag_id<TAB>barcode<TAB>umi

    Returns: (n_total_reads_processed, n_matched_reads_emitted)
    """
    if cfg.parts <= 0:
        raise ValueError("--parts must be >= 1")

    qtag_lut = build_qtag_lookup(QTAGS)

    run_dir = os.path.join(cfg.out_dir, run)
    _ensure_dir(run_dir)

    ndigits = max(3, len(str(cfg.parts - 1)))
    part_paths = [
        os.path.join(run_dir, f"keys_part_{i:0{ndigits}d}.tsv") for i in range(cfg.parts)
    ]

    fhs: List[TextIO] = []
    try:
        for p in part_paths:
            fhs.append(open(p, "wt", encoding="utf-8", newline="\n"))

        n_total = 0
        n_matched = 0
        import time

        t0 = time.time()

        def _msg(log_fh: Optional[TextIO], text: str) -> None:
            if log_fh is not None:
                log_fh.write(f"# {_timestamp_local()} [{run}] {text}\n")
                log_fh.flush()

        def _log_tsv(log_fh: Optional[TextIO], event: str) -> None:
            if log_fh is None:
                return
            elapsed = time.time() - t0
            rate = int(n_total / elapsed) if elapsed > 0 else 0
            pct = (100.0 * n_matched / n_total) if n_total else 0.0
            log_fh.write(
                "\t".join(
                    [
                        _timestamp_local(),
                        run,
                        event,
                        str(n_total),
                        str(n_matched),
                        f"{pct:.2f}",
                        str(rate),
                        f"{elapsed:.3f}",
                        "1",
                    ]
                )
                + "\n"
            )
            log_fh.flush()

        ensure_progress_log_header(cfg.log_file)
        with open(cfg.log_file, "at", encoding="utf-8", newline="\n") as log_fh:
            _msg(log_fh, f"Extracting keys for run={run} from {fastq_path} into {run_dir} (parts={cfg.parts})")
            _log_tsv(log_fh, "start")

            for seq in iter_fastq_sequences(
                fastq_path,
                strict_4line=cfg.strict_fastq,
                normalize_whitespace=True,
            ):
                n_total += 1
                m = FULL_RE.search(seq)
                if not m:
                    pass
                else:
                    umi = m.group(1)
                    barcode_motif = m.group(2)
                    qtag_obs = m.group(3)

                    # barcode is 8nt immediately after "CGA"
                    barcode = barcode_motif[3:11]

                    qtag_id = qtag_lut.get(qtag_obs, None)
                    if qtag_id is None:
                        # either no hit or ambiguous hit
                        pass
                    else:
                        pi = _stable_partition_index(qtag_id, barcode, cfg.parts)
                        fhs[pi].write(f"{qtag_id}\t{barcode}\t{umi}\n")
                        n_matched += 1

                if cfg.progress_every and (n_total % cfg.progress_every == 0):
                    elapsed = time.time() - t0
                    rate = int(n_total / elapsed) if elapsed > 0 else 0
                    pct = (100.0 * n_matched / n_total) if n_total else 0.0
                    _msg(
                        log_fh,
                        f"[{run}] {n_total} reads; {n_matched} matched; {pct:.2f}% matched; ~{rate} reads/sec",
                    )
                    _log_tsv(log_fh, "progress")

            _log_tsv(log_fh, "done")
        return n_total, n_matched

    finally:
        for fh in fhs:
            try:
                fh.close()
            except Exception:
                pass
