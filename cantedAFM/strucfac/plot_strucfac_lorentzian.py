#!/usr/bin/env python3
"""Plot Lorentzian-broadened structure factor for a selected momentum sector."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np


MomentumArray = np.ndarray
FloatArray = np.ndarray


def read_long_form(path: Path) -> Tuple[MomentumArray, FloatArray, FloatArray]:
    """Read momentum, energy, weight columns from long-form CSV text."""
    momenta: List[str] = []
    energies: List[float] = []
    weights: List[float] = []

    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if header is None:
            raise ValueError(f"Input file is empty: {path}")

        for row in reader:
            if len(row) < 3:
                continue
            momenta.append(row[0].strip())
            energies.append(float(row[1]))
            weights.append(float(row[2]))

    if not energies:
        raise ValueError(f"No data rows found in {path}")

    return np.array(momenta, dtype=object), np.array(energies), np.array(weights)


def filter_by_momentum(
    momenta: MomentumArray,
    energies: FloatArray,
    weights: FloatArray,
    momentum: str,
) -> Tuple[FloatArray, FloatArray]:
    """Select energy/weight rows for the chosen momentum.

    If momentum has no dot (e.g. "Gamma"), it matches labels by prefix before
    the first dot (e.g. "Gamma.C1.A"). If it includes a dot, exact matching is used.
    """
    if "." in momentum:
        mask = np.array([label == momentum for label in momenta], dtype=bool)
    else:
        mask = np.array([str(label).split(".", 1)[0] == momentum for label in momenta], dtype=bool)

    return energies[mask], weights[mask]


def lorentzian_sum(
    omegas: FloatArray,
    energies: FloatArray,
    weights: FloatArray,
    eps: float,
) -> FloatArray:
    """Compute sum_i w_i * eps/pi / ((omega - E_i)^2 + eps^2)."""
    if eps <= 0.0:
        raise ValueError("eps must be > 0")

    diff = omegas[:, None] - energies[None, :]
    kernel = (eps / np.pi) / (diff * diff + eps * eps)
    return kernel @ weights


def build_omega_grid(
    energies: FloatArray,
    eps: float,
    omega_min: float | None,
    omega_max: float | None,
    npoints: int,
) -> FloatArray:
    """Construct omega grid around pole energies unless explicit bounds are provided."""
    if npoints < 2:
        raise ValueError("npoints must be >= 2")

    if omega_min is None:
        omega_min = float(np.min(energies) - 10.0 * eps)
    if omega_max is None:
        omega_max = float(np.max(energies) + 10.0 * eps)
    if omega_max <= omega_min:
        raise ValueError("omega_max must be larger than omega_min")

    return np.linspace(omega_min, omega_max, npoints)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read long-form structure-factor data, keep one momentum sector, "
            "sum Lorentzians, and plot the broadened spectrum."
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=Path("teststrucfac.outfile.txt"),
        help="Path to long-form data file (default: teststrucfac.outfile.txt)",
    )
    parser.add_argument(
        "--momentum",
        type=str,
        default="Gamma",
        help=(
            "Momentum selector (default: Gamma). "
            "Use prefix (Gamma) or full label (Gamma.C1.A)."
        ),
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=0.04,
        help="Lorentzian broadening epsilon (default: 0.04)",
    )
    parser.add_argument(
        "--omega-min",
        type=float,
        default=None,
        help="Optional omega minimum for plotting grid",
    )
    parser.add_argument(
        "--max-frequency",
        type=float,
        default=2.0,
        help="Maximum frequency (default upper omega bound) for plotting grid (default: 2.0)",
    )
    parser.add_argument(
        "--omega-max",
        type=float,
        default=None,
        help="Optional explicit omega maximum for plotting grid (overrides --max-frequency)",
    )
    parser.add_argument(
        "--npoints",
        type=int,
        default=2000,
        help="Number of omega grid points (default: 2000)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output image path (e.g. spectrum.png)",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open an interactive plot window",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    momenta, energies, weights = read_long_form(args.input)
    sel_energies, sel_weights = filter_by_momentum(momenta, energies, weights, args.momentum)

    if sel_energies.size == 0:
        unique = sorted(set(str(m) for m in momenta))
        raise ValueError(
            f"No rows found for momentum '{args.momentum}'. "
            f"Available labels include: {', '.join(unique[:12])}"
        )

    omega_max = args.omega_max if args.omega_max is not None else args.max_frequency

    omegas = build_omega_grid(
        sel_energies,
        args.eps,
        args.omega_min,
        omega_max,
        args.npoints,
    )
    spectrum = lorentzian_sum(omegas, sel_energies, sel_weights, args.eps)

    plt.figure(figsize=(8.0, 4.8))
    plt.plot(omegas, spectrum, lw=1.5)
    plt.xlabel("omega")
    plt.ylabel("S(omega)")
    plt.title(f"Lorentzian-broadened structure factor ({args.momentum}, eps={args.eps:g})")
    plt.grid(alpha=0.3)
    plt.tight_layout()

    if args.output is not None:
        plt.savefig(args.output, dpi=300)

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()

