from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional


def log2(x: float) -> float:
    return math.log(x, 2)


def hamming(a: str, b: str) -> int:
    if a == b:
        return 0
    if len(a) != len(b):
        raise ValueError(f"hamming() unequal lengths: {len(a)} vs {len(b)}")
    return sum(1 for ca, cb in zip(a, b) if ca != cb)


def fmt_float(x: Optional[float], *, float_fmt: str = "{:.12g}") -> str:
    if x is None:
        return ""
    return float_fmt.format(x)


@dataclass(frozen=True)
class ThresholdResult:
    threshold_index0: Optional[int]
    note: str
