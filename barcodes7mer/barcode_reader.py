from __future__ import annotations

import os
import re
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .fastq import iter_fastq_sequences
from .util import hamming


# IMPORTANT: corrected barcode motif to match your validated historical outputs
BARCODE_RE_STR = r"CGA[ACTG]{3}C[ACTG]{4}AATTCGATGG"
MCOUNT_RE_STR  = r"[ATCG]{0,3}C[ACTG]{3}C[ACTG]{3}C[ACTG]{3}GCGCAACGCG"

FULL_RE = re.compile(
    rf"({MCOUNT_RE_STR})[ATCG]+({BARCODE_RE_STR})[ATCG]+TGGTGTTCAAGCTT([ATCG]{{12}})"
)

QTAGS: Dict[str, str] = {
    "TCGGCTAGATGT": "19",
    "AGGAACACCAAG": "23",
    "TCGCCGAGCAGT": "22",
    "CGAGCGCGAGGA": "24",
    "TGGCGAATATGG": "25",
    "TCTTCTACAACA": "26",
    "AGCACGCCTTGT": "27",
    "GCAACTTCTTCA": "26_2",
    "AAGAAGTCCAAC": "17",
}
QTAG_SEQS: List[str] = list(QTAGS.keys())


@dataclass(frozen=True)
class Sample:
    run: str
    r1_fastq_gz: str


def read_sample_list(sample_list_tsv: str) -> List[Sample]:
    """
    Accept either:
      run<TAB>path
    or
      run<whitespace>path
    """
    out: List[Sample] = []
    with open(sample_list_tsv, "rt", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                continue

            run: Optional[str] = None
            path: Optional[str] = None

            if "\t" in line:
                parts = line.split("\t", 1)
                if len(parts) == 2:
                    run = parts[0].strip()
                    path = parts[1].strip()
            else:
                parts = line.split()
                if len(parts) >= 2:
                    run = parts[0].strip()
                    path = parts[1].strip()

            if not run or not path:
                continue

            out.append(Sample(run=run, r1_fastq_gz=path))

    if not out:
        raise ValueError(f"No valid sample entries found in {sample_list_tsv}")
    return out


def bin_qtags(obs: str, candidates: Sequence[str], *, max_mm: int = 2) -> List[str]:
    hits: List[str] = []
    for cand in candidates:
        if len(cand) != len(obs):
            continue
        if hamming(cand, obs) <= max_mm:
            hits.append(cand)
    return hits


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


def process_sample(
    run: str,
    input_fastq: str,
    *,
    progress_every: int = 1_000_000,
    log_fh=None,
    strict_fastq: bool = True,
) -> Tuple[int, int]:
    reads_out = f"{run}_reads.tsv"
    chimera_out = f"{run}_chimera_data.txt"
    unmatched_out = f"{run}_unmatched_regex.tsv"

    chimera_counts: Dict[str, Dict[str, int]] = {}

    n_total = 0
    n_matched = 0
    t0 = time.time()

    def msg(text: str) -> None:
        if log_fh is not None:
            log_fh.write(f"# {_timestamp_local()} [{run}] {text}\n")
            log_fh.flush()

    def log_tsv(event: str) -> None:
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

    msg(f"Processing {run} from input file {input_fastq}...")
    log_tsv("start")

    with open(reads_out, "wt", encoding="utf-8", newline="\n") as of_reads, \
         open(chimera_out, "wt", encoding="utf-8", newline="\n") as of_chim, \
         open(unmatched_out, "wt", encoding="utf-8", newline="\n") as of_unm:

        of_reads.write("\t".join(["run", "qtag_id", "barcode", "mcount_tag", "read_sequence"]) + "\n")
        of_chim.write("\t".join(["run", "barcode", "qtag_id", "read_count"]) + "\n")
        of_unm.write("\t".join(["run", "read_sequence"]) + "\n")

        for seq in iter_fastq_sequences(
            input_fastq,
            strict_4line=strict_fastq,
            normalize_whitespace=True,
        ):
            n_total += 1

            m = FULL_RE.search(seq)
            if not m:
                of_unm.write(f"{run}\t{seq}\n")
            else:
                mcount_tag = m.group(1)
                barcode_motif = m.group(2)
                qtag_obs = m.group(3)

                barcode = barcode_motif[3:11]  # Perl @c[3..10]

                hits = bin_qtags(qtag_obs, QTAG_SEQS, max_mm=2)
                if len(hits) != 1:
                    pass
                else:
                    qtag_id = QTAGS.get(hits[0])
                    if qtag_id is None:
                        pass
                    else:
                        of_reads.write("\t".join([run, qtag_id, barcode, mcount_tag, seq]) + "\n")
                        chimera_counts.setdefault(barcode, {}).setdefault(qtag_id, 0)
                        chimera_counts[barcode][qtag_id] += 1
                        n_matched += 1

            if progress_every and (n_total % progress_every == 0):
                elapsed = time.time() - t0
                rate = int(n_total / elapsed) if elapsed > 0 else 0
                pct = (100.0 * n_matched / n_total) if n_total else 0.0
                msg(f"  [{run}] {n_total} reads processed; {n_matched} matched; {pct:.2f}% matched; ~{rate} reads/sec")
                log_tsv("progress")

        for bc in sorted(chimera_counts.keys()):
            qids = sorted(chimera_counts[bc].keys())
            if len(qids) <= 1:
                continue
            for qid in qids:
                of_chim.write(f"{run}\t{bc}\t{qid}\t{chimera_counts[bc][qid]}\n")

    elapsed = time.time() - t0
    pct = (100.0 * n_matched / n_total) if n_total else 0.0
    msg(f"Done {run}: {n_total} reads processed; {n_matched} matched; {pct:.2f}% matched; elapsed {elapsed:.1f}s")
    log_tsv("done")
    return n_total, n_matched
