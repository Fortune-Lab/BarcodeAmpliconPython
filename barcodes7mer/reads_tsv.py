from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import DefaultDict, Dict, Set


@dataclass
class ReadsData:
    data: Dict[str, Dict[str, Set[str]]]          # data[qtag][barcode] = set(umi)
    umi_counts: Dict[str, Dict[str, int]]         # umi_counts[umi][barcode] = count


def load_reads_tsv(path: str) -> ReadsData:
    data: DefaultDict[str, DefaultDict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    umi_counts: DefaultDict[str, DefaultDict[str, int]] = defaultdict(lambda: defaultdict(int))

    with open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            a = line.split("\t")
            if a and a[0] == "run":
                continue
            if len(a) < 5:
                continue
            _run, qtag, barcode, umi = a[0], a[1], a[2], a[3]
            umi_counts[umi][barcode] += 1
            data[qtag][barcode].add(umi)

    return ReadsData(data=dict(data), umi_counts={k: dict(v) for k, v in umi_counts.items()})


def filter_umi_inplace(
    data: Dict[str, Dict[str, Set[str]]],
    umi_counts: Dict[str, Dict[str, int]],
    *,
    enabled: bool,
) -> None:
    if not enabled:
        return

    for qtag in list(data.keys()):
        for barcode in list(data[qtag].keys()):
            keep = set()
            for umi in data[qtag][barcode]:
                if umi_counts.get(umi, {}).get(barcode, 0) <= 1:
                    keep.add(umi)
            data[qtag][barcode] = keep


def mcounts(data: Dict[str, Dict[str, Set[str]]]) -> Dict[str, Dict[str, int]]:
    out: Dict[str, Dict[str, int]] = {}
    for qtag, by_bc in data.items():
        for bc, umis in by_bc.items():
            out.setdefault(qtag, {})[bc] = len(umis)
    return out
