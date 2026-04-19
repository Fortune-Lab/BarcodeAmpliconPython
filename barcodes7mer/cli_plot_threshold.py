from __future__ import annotations

import os


def main(argv=None) -> int:
    # Make plotting safe on headless HPC nodes by default
    os.environ.setdefault("MPLBACKEND", "Agg")

    from .plot_threshold import main as plot_main
    return int(plot_main(argv))
