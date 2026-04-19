from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .chimera import chimera_percent, load_chimera_counts_from_cwd
from .reads_tsv import filter_umi_inplace, load_reads_tsv, mcounts as compute_mcounts
from .util import ThresholdResult, fmt_float, hamming, log2


@dataclass(frozen=True)
class ThresholdRow:
    run: str
    index1: int
    qbid: str
    counts: int
    norm: float
    slope: Optional[float]
    curvature: Optional[float]
    percent: str
    percent_chimera: str  # numeric "x.yyyy" OR "NA"
    is_threshold: str
    note: str


def remove_barcodes_one_off(
    values: List[Tuple[int, str, str]],  # (count, barcode, qtag)
    *,
    add_chimeras: bool,
    run: str,
    out_prefix: Optional[str] = None,
) -> List[Tuple[int, str, str]]:
    cache = set()
    out_path = f"{run}.chimera_one_off.txt" if out_prefix is None else f"{out_prefix}_{run}.chimera_one_off.txt"

    with open(out_path, "wt", encoding="utf-8", newline="\n") as df:
        for i in range(len(values)):
            if i in cache:
                continue
            for j in range(i + 1, len(values)):
                if j in cache:
                    continue
                mm = hamming(values[i][1], values[j][1])
                if mm <= 1:
                    df.write(
                        "\t".join(
                            [run, values[i][1], str(values[i][0]), values[j][1], str(values[j][0]), str(mm)]
                        )
                        + "\n"
                    )
                    cache.add(j)
                    if add_chimeras:
                        values[i] = (values[i][0] + values[j][0], values[i][1], values[i][2])

    return [values[i] for i in range(len(values)) if i not in cache]


def derivative1_uniform(y: Sequence[float]) -> List[float]:
    n = len(y)
    if n == 0:
        return []
    if n == 1:
        return [0.0]
    if n == 2:
        d = y[1] - y[0]
        return [d, d]

    d = [0.0] * n
    d[0] = (-3.0 * y[0] + 4.0 * y[1] - y[2]) / 2.0
    for i in range(1, n - 1):
        d[i] = (y[i + 1] - y[i - 1]) / 2.0
    d[n - 1] = (3.0 * y[n - 1] - 4.0 * y[n - 2] + y[n - 3]) / 2.0
    return d


def derivative2_uniform(y: Sequence[float]) -> List[float]:
    n = len(y)
    if n == 0:
        return []
    if n < 3:
        return [0.0] * n
    if n == 3:
        dd = (y[2] - 2.0 * y[1] + y[0])
        return [dd, dd, dd]

    dd = [0.0] * n
    dd[0] = (2.0 * y[0] - 5.0 * y[1] + 4.0 * y[2] - y[3])
    for i in range(1, n - 1):
        dd[i] = (y[i + 1] - 2.0 * y[i] + y[i - 1])
    dd[n - 1] = (2.0 * y[n - 1] - 5.0 * y[n - 2] + 4.0 * y[n - 3] - y[n - 4])
    return dd


def endpoint_slopes(y: Sequence[float], window: int, *, shoulder: int = 0) -> List[Optional[float]]:
    n = len(y)
    slope: List[Optional[float]] = [None] * n
    if window <= 1:
        return slope

    for i in range(n):
        idx1 = i + 1
        if shoulder and idx1 > shoulder:
            continue
        if idx1 < window:
            continue
        left_idx1 = idx1 - window + 1
        slope[i] = (y[idx1 - 1] - y[left_idx1 - 1]) / float(window - 1)
    return slope


def derivative_of_slopes(s: Sequence[Optional[float]]) -> List[Optional[float]]:
    n = len(s)
    d: List[Optional[float]] = [None] * n
    for i in range(n):
        if s[i] is None:
            continue
        if i > 0 and i < n - 1 and s[i - 1] is not None and s[i + 1] is not None:
            d[i] = (s[i + 1] - s[i - 1]) / 2.0
        elif i > 0 and s[i - 1] is not None:
            d[i] = s[i] - s[i - 1]
        elif i < n - 1 and s[i + 1] is not None:
            d[i] = s[i + 1] - s[i]
    return d


def rolling_slope(y: Sequence[float], window: int, *, shoulder: int = 0, python_compat: bool = True) -> List[Optional[float]]:
    n = len(y)
    out: List[Optional[float]] = [None] * n
    if window <= 1 or n < 3:
        return out

    pref_y = [0.0]
    pref_xy = [0.0]
    for i in range(n):
        x = i + 1
        pref_y.append(pref_y[-1] + y[i])
        pref_xy.append(pref_xy[-1] + x * y[i])

    def sum_int(a: int, b: int) -> float:
        return (b * (b + 1) / 2.0) - ((a - 1) * a / 2.0)

    def sum_sq(a: int, b: int) -> float:
        def ss(k: int) -> float:
            return k * (k + 1) * (2 * k + 1) / 6.0
        return ss(b) - ss(a - 1)

    for i in range(n):
        idx = i + 1
        if shoulder and idx > shoulder:
            continue
        if idx < window:
            continue

        if python_compat:
            a = idx - window + 2
            b = idx
        else:
            a = idx - window + 1
            b = idx

        if a < 1:
            continue
        nwin = b - a + 1
        if nwin < 2:
            continue

        start0, end0 = a - 1, b - 1
        sum_y = pref_y[end0 + 1] - pref_y[start0]
        sum_xy = pref_xy[end0 + 1] - pref_xy[start0]
        sum_x = sum_int(a, b)
        sum_x2 = sum_sq(a, b)

        den = nwin * sum_x2 - sum_x * sum_x
        if den == 0:
            continue

        m = (nwin * sum_xy - sum_x * sum_y) / den
        out[i] = m

    return out


def fill_undef_with_last(vals: Sequence[Optional[float]]) -> List[float]:
    first = None
    for v in vals:
        if v is not None:
            first = v
            break
    if first is None:
        return [0.0] * len(vals)

    out: List[float] = []
    last = first
    for v in vals:
        if v is None:
            out.append(last)
        else:
            out.append(v)
            last = v
    return out


def argmin_defined_in_range(arr: Sequence[Optional[float]], *, start_idx1: int, stop_idx1: int, min_idx1: int) -> Optional[int]:
    best_i = None
    best_v = None
    for i, v in enumerate(arr):
        idx1 = i + 1
        if idx1 < min_idx1 or idx1 < start_idx1 or idx1 > stop_idx1:
            continue
        if v is None:
            continue
        if best_v is None or v < best_v:
            best_v = v
            best_i = i
    return best_i


def compute_threshold(
    *,
    series_for_argmin: Sequence[Optional[float]],
    window: int,
    n_points: int,
    shoulder: int,
    search_start: int,
    search_stop: int,
    min_index_allowed: int,
    disable_when_window_le_1: bool,
    min_start_idx1: int = 1,
) -> ThresholdResult:
    stop_default = shoulder if (shoulder and shoulder < n_points) else n_points
    start_default = min_start_idx1

    sstart = search_start if search_start else start_default
    sstop = search_stop if search_stop else stop_default

    sstart = max(1, sstart, min_index_allowed, min_start_idx1)
    sstop = min(n_points, sstop)

    if disable_when_window_le_1 and window <= 1:
        return ThresholdResult(None, "window_disabled")
    if sstop <= sstart:
        return ThresholdResult(None, "search_interval_empty")

    i0 = argmin_defined_in_range(series_for_argmin, start_idx1=sstart, stop_idx1=sstop, min_idx1=min_index_allowed)
    if i0 is None:
        return ThresholdResult(None, "no_valid_threshold_candidates")
    return ThresholdResult(i0, "")


def build_chimera_percent_map_from_cwd() -> Dict[str, Dict[str, float]]:
    chim_counts = load_chimera_counts_from_cwd()
    return chimera_percent(chim_counts)


def threshold_table_for_reads_file(
    reads_file: str,
    *,
    method: str,  # "derivative" or "regression"
    min_reads: int,
    percent_cutoff: int,
    add_chimeras: bool,
    filter_umi: bool,
    assay_mode: bool,
    window: int,
    shoulder: int,
    search_start: int,
    search_stop: int,
    python_compat: bool,
    min_index_allowed: int,
    out_prefix_for_chimera_oneoff: Optional[str],
    float_fmt: str,
    chimera_pct_map: Dict[str, Dict[str, float]],
) -> Tuple[str, List[ThresholdRow], Optional[str]]:
    import os

    base = os.path.basename(reads_file)
    run = base[: -len("_reads.tsv")] if base.endswith("_reads.tsv") else os.path.splitext(base)[0]

    rd = load_reads_tsv(reads_file)
    filter_umi_inplace(rd.data, rd.umi_counts, enabled=filter_umi)
    mc = compute_mcounts(rd.data)

    values: List[Tuple[int, str, str]] = []
    for qtag, by_bc in mc.items():
        for bc, n in by_bc.items():
            if n and n > 0:
                values.append((n, bc, qtag))
    values.sort(key=lambda t: t[0], reverse=True)

    if assay_mode:
        values = remove_barcodes_one_off(values, add_chimeras=add_chimeras, run=run, out_prefix=out_prefix_for_chimera_oneoff)
        values.sort(key=lambda t: t[0], reverse=True)

    total = sum(v[0] for v in values)
    if total <= min_reads:
        return run, [], f"sum_reads<={min_reads}"
    if len(values) < 3:
        return run, [], "too_few_points"

    counts = [v[0] for v in values]
    barcode = [v[1] for v in values]
    qtag = [v[2] for v in values]

    percent = [f"{(100.0 * c / total):.4f}" for c in counts]
    norm = [log2(c / total) for c in counts]

    if method == "derivative":
        if window <= 1:
            slope = [float(x) for x in derivative1_uniform(norm)]
            curvature = [float(x) for x in derivative2_uniform(norm)]
        else:
            slope = endpoint_slopes(norm, window, shoulder=shoulder)
            curvature = derivative_of_slopes(slope)

        thr = compute_threshold(
            series_for_argmin=curvature,
            window=window,
            n_points=len(values),
            shoulder=shoulder,
            search_start=search_start,
            search_stop=search_stop,
            min_index_allowed=min_index_allowed,
            disable_when_window_le_1=False,
            min_start_idx1=min_index_allowed,
        )

    elif method == "regression":
        slope = rolling_slope(norm, window, shoulder=shoulder, python_compat=python_compat)

        slope_fill = fill_undef_with_last(slope)
        d_slope_all = derivative1_uniform(slope_fill)
        curvature = [d_slope_all[i] if slope[i] is not None else None for i in range(len(slope))]

        thr = compute_threshold(
            series_for_argmin=slope,  # argmin(rolling slope), Perl-compatible
            window=window,
            n_points=len(values),
            shoulder=shoulder,
            search_start=search_start,
            search_stop=search_stop,
            min_index_allowed=min_index_allowed,
            disable_when_window_le_1=True,
            min_start_idx1=(window if window > 1 else 1),
        )
    else:
        raise ValueError("method must be 'derivative' or 'regression'")

    chimera_is_na = (not assay_mode) and (int(percent_cutoff) == 0)

    rows: List[ThresholdRow] = []
    for i in range(len(values)):
        idx1 = i + 1
        qbid = f"{qtag[i]}{barcode[i]}"

        if chimera_is_na:
            chim_s = "NA"
        else:
            chim = chimera_pct_map.get(run, {}).get(barcode[i], 0.0)
            chim_s = f"{chim:.4f}"

        is_thr = "TRUE" if (thr.threshold_index0 is not None and i == thr.threshold_index0) else "FALSE"
        note = thr.note if is_thr == "TRUE" else ""

        rows.append(
            ThresholdRow(
                run=run,
                index1=idx1,
                qbid=qbid,
                counts=counts[i],
                norm=norm[i],
                slope=slope[i],
                curvature=curvature[i],
                percent=percent[i],
                percent_chimera=chim_s,
                is_threshold=is_thr,
                note=note,
            )
        )

    return run, rows, None


def write_threshold_tsv(out_path: str, all_rows: List[ThresholdRow], *, float_fmt: str) -> None:
    header = ["run", "index", "qbid", "counts", "norm", "slope", "curvature", "percent", "percent_chimera", "is_threshold", "note"]
    with open(out_path, "wt", encoding="utf-8", newline="\n") as of:
        of.write("\t".join(header) + "\n")
        for r in all_rows:
            of.write(
                "\t".join(
                    [
                        r.run,
                        str(r.index1),
                        r.qbid,
                        str(r.counts),
                        float_fmt.format(r.norm),
                        fmt_float(r.slope, float_fmt=float_fmt),
                        fmt_float(r.curvature, float_fmt=float_fmt),
                        r.percent,
                        r.percent_chimera,
                        r.is_threshold,
                        r.note,
                    ]
                )
                + "\n"
            )
