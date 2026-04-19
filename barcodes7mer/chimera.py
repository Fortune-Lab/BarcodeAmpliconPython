from __future__ import annotations

import glob
from typing import Dict, List


def load_chimera_counts_from_cwd(pattern: str = "*_chimera_data.txt") -> Dict[str, Dict[str, List[int]]]:
    out: Dict[str, Dict[str, List[int]]] = {}
    for path in glob.glob(pattern):
        with open(path, "rt", encoding="utf-8") as fh:
            first = True
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if first:
                    first = False
                    continue
                a = line.split("\t")
                if len(a) < 4:
                    continue
                run, barcode, _qtag_id, count_s = a[0], a[1], a[2], a[3]
                try:
                    c = int(count_s)
                except ValueError:
                    continue
                out.setdefault(run, {}).setdefault(barcode, []).append(c)
    return out


def chimera_percent(chimera_counts: Dict[str, Dict[str, List[int]]]) -> Dict[str, Dict[str, float]]:
    out: Dict[str, Dict[str, float]] = {}
    for run, by_bc in chimera_counts.items():
        for bc, counts in by_bc.items():
            if not counts:
                continue
            s = sum(counts)
            if s == 0:
                continue
            m = max(counts)
            out.setdefault(run, {})[bc] = 100.0 * ((s - m) / s)
    return out
