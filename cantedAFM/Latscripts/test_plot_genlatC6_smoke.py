#!/usr/bin/env python3
"""Smoke test for plot_genlatC6.py in headless mode."""

import os
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch

from genlatC6 import build_torus, compute_nnn_bonds_chiral, compute_nn_bonds, site_to_cart
from plot_genlatC6 import _bond_segment_coords, _display_coords, plot_real_space


def run_case(
    space,
    a,
    b,
    out_file,
    show_bonds=True,
    show_nnn_bonds=True,
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
    if space == "real" and show_bonds and not show_nnn_bonds:
        cmd.append("--no-show-nnn-bonds")
    if space == "real" and show_rotation_labels:
        cmd.append("--show-rotation-labels")
    if space == "real" and display_coordinates is not None:
        cmd.extend(["--display-coordinates", display_coordinates])
    subprocess.run(cmd, check=True, cwd=root, env=env)


def assert_has_wrapped_bond(a, b, display_coordinates, bond_fn, label):
    _, sites, T, t_inv_n, ncell, _ = build_torus(a, b)
    coords_canonical = np.array([site_to_cart(site) for site in sites])
    coords = _display_coords(coords_canonical, T, display_coordinates)

    found_wrapped = False
    for i, j in bond_fn(sites, t_inv_n, ncell):
        segment = _bond_segment_coords(coords, i, j, T)
        direct = np.linalg.norm(coords[j] - coords[i])
        wrapped = np.linalg.norm(segment[1] - segment[0])
        if not np.allclose(segment[1], coords[j]) and wrapped + 1e-10 < direct:
            found_wrapped = True
            break

    assert found_wrapped, (
        f"Expected at least one wrapped {label} bond for (a={a}, b={b}, "
        f"display={display_coordinates})."
    )


def assert_nnn_arrow_count(a, b, display_coordinates, show_nnn_bonds, expected_count):
    fig = plot_real_space(
        a,
        b,
        show_nnn_bonds=show_nnn_bonds,
        display_coordinates=display_coordinates,
    )
    try:
        ax = fig.axes[0]
        arrows = [patch for patch in ax.patches if isinstance(patch, FancyArrowPatch)]
        assert len(arrows) == expected_count, (
            f"Expected {expected_count} NNN arrows for (a={a}, b={b}, "
            f"display={display_coordinates}, show_nnn_bonds={show_nnn_bonds}), "
            f"got {len(arrows)}."
        )
    finally:
        plt.close(fig)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        real_png = tmpdir / "real.png"
        real_no_bonds_png = tmpdir / "real_no_bonds.png"
        real_no_nnn_png = tmpdir / "real_no_nnn.png"
        real_rot_labels_png = tmpdir / "real_rotation_labels.png"
        real_ws_png = tmpdir / "real_wigner_seitz.png"
        mom_png = tmpdir / "momentum.png"

        run_case("real", 3, 0, real_png)
        run_case("real", 3, 0, real_no_bonds_png, show_bonds=False)
        run_case("real", 3, 0, real_no_nnn_png, show_nnn_bonds=False)
        run_case("real", 3, 0, real_rot_labels_png, show_rotation_labels=True)
        run_case("real", 3, 0, real_ws_png, display_coordinates="wigner-seitz")
        run_case("momentum", 4, 2, mom_png)

        assert real_png.exists() and real_png.stat().st_size > 0
        assert real_no_bonds_png.exists() and real_no_bonds_png.stat().st_size > 0
        assert real_no_nnn_png.exists() and real_no_nnn_png.stat().st_size > 0
        assert real_rot_labels_png.exists() and real_rot_labels_png.stat().st_size > 0
        assert real_ws_png.exists() and real_ws_png.stat().st_size > 0
        assert mom_png.exists() and mom_png.stat().st_size > 0

        _, sites, _, t_inv_n, ncell, _ = build_torus(3, 0)
        expected_nnn = len(compute_nnn_bonds_chiral(sites, t_inv_n, ncell))

        for display_coordinates in ("canonical", "wigner-seitz"):
            assert_has_wrapped_bond(3, 0, display_coordinates, compute_nn_bonds, "NN")
            assert_has_wrapped_bond(3, 0, display_coordinates, compute_nnn_bonds_chiral, "NNN")
            assert_nnn_arrow_count(3, 0, display_coordinates, show_nnn_bonds=True, expected_count=expected_nnn)
            assert_nnn_arrow_count(3, 0, display_coordinates, show_nnn_bonds=False, expected_count=0)

    print("Smoke test passed.")


if __name__ == "__main__":
    main()

