from __future__ import annotations

import glob
import os
import subprocess
from typing import Iterator, List, Optional, Tuple


def _iter_partition_files(keys_dir: str) -> List[str]:
    paths = sorted(glob.glob(os.path.join(keys_dir, "keys_part_*.tsv")))
    if not paths:
        raise FileNotFoundError(f"No keys_part_*.tsv found in {keys_dir}")
    return paths


def _sorted_lines_via_sort(
    in_path: str,
    *,
    tmp_dir: str,
    sort_mem: str,
    sort_bin: str = "sort",
) -> Iterator[str]:
    """
    Yields sorted lines from GNU sort:
      sort -T tmp_dir -S sort_mem -t '\t' -k1,1 -k2,2 -k3,3 in_path
    """
    env = dict(os.environ)
    env["LC_ALL"] = "C"

    cmd = [
        sort_bin,
        "-T",
        tmp_dir,
        "-S",
        sort_mem,
        "-t",
        "\t",
        "-k1,1",
        "-k2,2",
        "-k3,3",
        in_path,
    ]

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        env=env,
    )
    assert proc.stdout is not None
    assert proc.stderr is not None

    try:
        for line in proc.stdout:
            yield line
    finally:
        err = proc.stderr.read()
        rc = proc.wait()
        if rc != 0:
            raise RuntimeError(f"sort failed (rc={rc}) on {in_path}\nSTDERR:\n{err}")


def reduce_partition_to_mcounts(
    part_path: str,
    *,
    run: str,
    tmp_dir: str,
    sort_mem: str,
    sort_bin: str = "sort",
) -> Iterator[Tuple[str, str, str, int]]:
    """
    Input partition line: qtag_id<TAB>barcode<TAB>umi
    Output rows: (run, qtag_id, barcode, mcount)
    """
    cur_qtag: Optional[str] = None
    cur_bc: Optional[str] = None
    last_umi: Optional[str] = None
    mcount = 0

    for line in _sorted_lines_via_sort(part_path, tmp_dir=tmp_dir, sort_mem=sort_mem, sort_bin=sort_bin):
        line = line.rstrip("\n")
        if not line:
            continue
        a = line.split("\t")
        if len(a) < 3:
            continue

        qtag_id, barcode, umi = a[0], a[1], a[2]

        if (cur_qtag is None) or (qtag_id != cur_qtag) or (barcode != cur_bc):
            if cur_qtag is not None:
                yield (run, cur_qtag, cur_bc, mcount)  # type: ignore[arg-type]
            cur_qtag, cur_bc = qtag_id, barcode
            last_umi = None
            mcount = 0

        if umi != last_umi:
            mcount += 1
            last_umi = umi

    if cur_qtag is not None:
        yield (run, cur_qtag, cur_bc, mcount)  # type: ignore[arg-type]


def reduce_keys_dir_to_mcounts(
    keys_dir: str,
    *,
    run: str,
    out_path: str,
    tmp_dir: str,
    sort_mem: str = "8G",
    sort_bin: str = "sort",
    only_part: Optional[int] = None,
) -> None:
    """
    Reduce all partitions in keys_dir into out_path.
    If only_part is provided, reduce just that partition index (0-based in sorted list of keys_part_*.tsv).
    """
    part_files = _iter_partition_files(keys_dir)

    if only_part is not None:
        if only_part < 0 or only_part >= len(part_files):
            raise ValueError(f"only_part={only_part} out of range (0..{len(part_files)-1})")
        part_files = [part_files[only_part]]

    out_dir = os.path.dirname(os.path.abspath(out_path)) or "."
    os.makedirs(out_dir, exist_ok=True)

    with open(out_path, "wt", encoding="utf-8", newline="\n") as of:
        of.write("run\tqtag_id\tbarcode\tmcount\n")
        for part_path in part_files:
            for (run_s, qtag_id, barcode, mcount) in reduce_partition_to_mcounts(
                part_path,
                run=run,
                tmp_dir=tmp_dir,
                sort_mem=sort_mem,
                sort_bin=sort_bin,
            ):
                of.write(f"{run_s}\t{qtag_id}\t{barcode}\t{mcount}\n")
