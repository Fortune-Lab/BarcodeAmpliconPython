from __future__ import annotations

import gzip
import re
from typing import Iterator, TextIO

_WS_RE = re.compile(r"\s+")


def _open_text_maybe_gz(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


def iter_fastq_sequences(
    path: str,
    *,
    strict_4line: bool = True,
    normalize_whitespace: bool = True,
) -> Iterator[str]:
    """
    FASTQ iterator yielding sequences.

    strict_4line:
      Enforce standard 4-line FASTQ records (@, seq, +, qual).

    normalize_whitespace:
      Remove ANY whitespace inside sequence line (helps parity with seqtk-style normalization).
    """
    with _open_text_maybe_gz(path) as fh:
        rec = 0
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline()
            p = fh.readline()
            q = fh.readline()

            rec += 1
            if not (s and p and q):
                raise ValueError(f"{path}: truncated FASTQ record at record #{rec}")

            if strict_4line:
                if not h.startswith("@"):
                    raise ValueError(f"{path}: record #{rec}: header does not start with '@': {h.rstrip()}")
                if not p.startswith("+"):
                    raise ValueError(f"{path}: record #{rec}: plus-line does not start with '+': {p.rstrip()}")

            seq = s.rstrip("\n")
            if normalize_whitespace:
                seq = _WS_RE.sub("", seq)

            yield seq.upper()
