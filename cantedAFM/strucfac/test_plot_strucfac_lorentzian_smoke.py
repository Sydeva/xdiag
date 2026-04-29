#!/usr/bin/env python3
"""Headless smoke test for plot_strucfac_lorentzian.py."""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent
    script = root / "plot_strucfac_lorentzian.py"
    data = root / "teststrucfac.outfile.txt"

    with tempfile.TemporaryDirectory() as tmp:
        out_png = Path(tmp) / "spectrum.png"
        env = os.environ.copy()
        env["MPLBACKEND"] = "Agg"

        cmd = [
            "python3",
            str(script),
            str(data),
            "--momentum",
            "Gamma",
            "--eps",
            "0.04",
            "--output",
            str(out_png),
            "--no-show",
        ]
        subprocess.run(cmd, check=True, cwd=root, env=env)

        assert out_png.exists(), "Output image was not created"
        assert out_png.stat().st_size > 0, "Output image is empty"

    print("Smoke test passed.")


if __name__ == "__main__":
    main()

