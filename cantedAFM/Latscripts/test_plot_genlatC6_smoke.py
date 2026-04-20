#!/usr/bin/env python3
"""Smoke test for plot_genlatC6.py in headless mode."""

import os
import subprocess
import tempfile
from pathlib import Path


def run_case(
    space,
    a,
    b,
    out_file,
    show_bonds=True,
    show_rotation_labels=False,
    display_coordinates=None,
):
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    root = Path(__file__).resolve().parent
    cmd = [
        "python3",
        str(root / "plot_genlatC6.py"),
        space,
        str(a),
        str(b),
        "--save",
        str(out_file),
        "--no-show",
    ]
    if space == "real" and not show_bonds:
        cmd.append("--no-show-bonds")
    if space == "real" and show_rotation_labels:
        cmd.append("--show-rotation-labels")
    if space == "real" and display_coordinates is not None:
        cmd.extend(["--display-coordinates", display_coordinates])
    subprocess.run(cmd, check=True, cwd=root, env=env)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        real_png = tmpdir / "real.png"
        real_no_bonds_png = tmpdir / "real_no_bonds.png"
        real_rot_labels_png = tmpdir / "real_rotation_labels.png"
        real_ws_png = tmpdir / "real_wigner_seitz.png"
        mom_png = tmpdir / "momentum.png"

        run_case("real", 3, 0, real_png)
        run_case("real", 3, 0, real_no_bonds_png, show_bonds=False)
        run_case("real", 3, 0, real_rot_labels_png, show_rotation_labels=True)
        run_case("real", 3, 0, real_ws_png, display_coordinates="wigner-seitz")
        run_case("momentum", 4, 2, mom_png)

        assert real_png.exists() and real_png.stat().st_size > 0
        assert real_no_bonds_png.exists() and real_no_bonds_png.stat().st_size > 0
        assert real_rot_labels_png.exists() and real_rot_labels_png.stat().st_size > 0
        assert real_ws_png.exists() and real_ws_png.stat().st_size > 0
        assert mom_png.exists() and mom_png.stat().st_size > 0

    print("Smoke test passed.")


if __name__ == "__main__":
    main()

