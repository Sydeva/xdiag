#!/usr/bin/env python3
"""Plot the lowest Heisenberg eigenvalues by momentum and irrep."""

from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import cycle
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt


def parse_spectrum_file(path: Path) -> Tuple[List[float], List[str], List[str], List[str]]:
    """Return arrays of eigenvalues, momentum labels, irrep labels, and Sz labels."""
    eigenvalues: List[float] = []
    sz_labels: List[str] = []
    momenta: List[str] = []
    irreps: List[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            fields = [item.strip() for item in line.split(",")]
            if len(fields) == 3:
                value_text, sz_text, label_text = fields
            elif len(fields) == 2:
                # Backward compatibility with the older format without explicit Sz.
                value_text, label_text = fields
                sz_text = "0"
            else:
                raise ValueError(f"Invalid data at line {line_number}: '{line}'")

            try:
                value = float(value_text)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid data at line {line_number}: '{line}'"
                ) from exc

            parts = [token for token in label_text.split(".") if token]
            if len(parts) < 2:
                raise ValueError(
                    f"Invalid label at line {line_number}: '{label_text}'"
                )

            momentum = parts[0]
            irrep = ".".join(parts[1:])

            eigenvalues.append(value)
            sz_labels.append(sz_text)
            momenta.append(momentum)
            irreps.append(irrep)

    return eigenvalues, momenta, irreps, sz_labels


def select_lowest_indices(
    eigenvalues: Sequence[float],
    momenta: Sequence[str],
    allowed_momenta: Sequence[str],
    n_lowest: int,
) -> List[int]:
    """Keep the indices of the lowest n eigenvalues for each momentum sector."""
    grouped_indices: Dict[str, List[int]] = defaultdict(list)
    allowed = set(allowed_momenta)

    for idx, momentum in enumerate(momenta):
        if momentum in allowed:
            grouped_indices[momentum].append(idx)

    selected: List[int] = []
    for momentum in allowed_momenta:
        indices = grouped_indices.get(momentum, [])
        indices = sorted(indices, key=lambda i: eigenvalues[i])[:n_lowest]
        selected.extend(indices)

    return selected


def build_irrep_styles(irreps: Sequence[str]) -> Dict[str, Tuple[str, str]]:
    """Assign a distinct marker and color pair to each irrep."""
    unique_irreps = sorted(set(irreps))
    markers = cycle(["o", "s", "^", "D", "v", "P", "X", "<", ">", "*"])
    cmap = plt.get_cmap("tab20")

    styles: Dict[str, Tuple[str, str]] = {}
    for idx, irrep in enumerate(unique_irreps):
        color = cmap(idx % cmap.N)
        marker = next(markers)
        styles[irrep] = (color, marker)

    return styles


def gamma_reference_energy(
    eigenvalues: Sequence[float], momenta: Sequence[str], gamma_label: str = "Gamma"
) -> float:
    """Return the minimum energy in the Gamma sector used for vertical offset."""
    gamma_values = [eigenvalues[i] for i, momentum in enumerate(momenta) if momentum == gamma_label]
    if not gamma_values:
        raise ValueError("Could not find any Gamma entries to define the energy offset.")
    return min(gamma_values)


def compute_horizontal_offsets(
    eigenvalues: Sequence[float],
    momenta: Sequence[str],
    sz_labels: Sequence[str],
    keep_indices: Sequence[int],
    momentum_order: Sequence[str],
    degeneracy_tol: float,
    jitter_width: float,
) -> Dict[int, float]:
    """Spread near-degenerate points at each momentum/Sz subcluster by x-offset."""
    grouped_indices: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    for idx in keep_indices:
        grouped_indices[(momenta[idx], sz_labels[idx])].append(idx)

    offsets = {idx: 0.0 for idx in keep_indices}
    for momentum in momentum_order:
        sz_values = sorted({sz_labels[i] for i in keep_indices if momenta[i] == momentum}, key=sz_sort_key)
        for sz in sz_values:
            indices = sorted(grouped_indices.get((momentum, sz), []), key=lambda i: eigenvalues[i])
            if not indices:
                continue

            cluster: List[int] = [indices[0]]
            clusters: List[List[int]] = []
            for idx in indices[1:]:
                if abs(eigenvalues[idx] - eigenvalues[cluster[-1]]) <= degeneracy_tol:
                    cluster.append(idx)
                else:
                    clusters.append(cluster)
                    cluster = [idx]
            clusters.append(cluster)

            for group in clusters:
                if len(group) == 1:
                    continue
                ordered_group = sorted(group)
                half_span = min(max(jitter_width, 0.0), 0.2)
                step = (2.0 * half_span) / (len(ordered_group) - 1)
                for position, idx in enumerate(ordered_group):
                    offsets[idx] = -half_span + position * step

    return offsets


def sz_sort_key(value: str) -> Tuple[int, float | str]:
    """Sort Sz labels numerically when possible, otherwise lexicographically."""
    try:
        return (0, float(value))
    except ValueError:
        return (1, value)


def sz_distance_label(sz_value: str, nsites: int) -> str:
    """Format |Sz - Nsites/2| for axis labels, with compact integer output."""
    try:
        distance = (float(sz_value) - 0.5 * nsites)
    except ValueError:
        return sz_value

    rounded = round(distance)
    if abs(distance - rounded) < 1e-12:
        return str(int(rounded))
    return f"{distance:.6g}"


def compute_sz_cluster_positions(
    keep_indices: Sequence[int],
    momenta: Sequence[str],
    sz_labels: Sequence[str],
    momentum_order: Sequence[str],
    nsites: int,
    momentum_span: float = 0.34,
) -> Tuple[Dict[Tuple[str, str], float], List[float], List[str]]:
    """Return x positions for momentum/Sz clusters and corresponding tick labels."""
    cluster_positions: Dict[Tuple[str, str], float] = {}
    tick_positions: List[float] = []
    tick_labels: List[str] = []

    for momentum_index, momentum in enumerate(momentum_order):
        momentum_sz_values = sorted(
            {sz_labels[idx] for idx in keep_indices if momenta[idx] == momentum},
            key=sz_sort_key,
        )
        if not momentum_sz_values:
            continue

        if len(momentum_sz_values) == 1:
            offsets = [0.0]
        else:
            half_span = min(max(momentum_span, 0.0), 0.45)
            step = (2.0 * half_span) / (len(momentum_sz_values) - 1)
            offsets = [-half_span + step * i for i in range(len(momentum_sz_values))]

        for sz, offset in zip(momentum_sz_values, offsets):
            xpos = momentum_index + offset
            cluster_positions[(momentum, sz)] = xpos
            tick_positions.append(xpos)
            tick_labels.append(f"{momentum}\n Sz ={sz_distance_label(sz, nsites)}")

    return cluster_positions, tick_positions, tick_labels


def plot_lowest_spectrum(
    eigenvalues: Sequence[float],
    momenta: Sequence[str],
    irreps: Sequence[str],
    sz_labels: Sequence[str],
    keep_indices: Sequence[int],
    momentum_order: Sequence[str],
    output: Path | None,
    show: bool,
    degeneracy_tol: float,
    jitter_width: float,
    nsites: int,
) -> None:
    """Create a scatter plot with marker style/color determined by irrep."""
    if not keep_indices:
        raise ValueError("No points selected for plotting.")

    x_positions = {momentum: i for i, momentum in enumerate(momentum_order)}
    cluster_positions, tick_positions, tick_labels = compute_sz_cluster_positions(
        keep_indices=keep_indices,
        momenta=momenta,
        sz_labels=sz_labels,
        momentum_order=momentum_order,
        nsites=nsites,
    )
    energy_shift = gamma_reference_energy(eigenvalues=eigenvalues, momenta=momenta)
    horizontal_offsets = compute_horizontal_offsets(
        eigenvalues=eigenvalues,
        momenta=momenta,
        sz_labels=sz_labels,
        keep_indices=keep_indices,
        momentum_order=momentum_order,
        degeneracy_tol=degeneracy_tol,
        jitter_width=jitter_width,
    )
    selected_irreps = [irreps[i] for i in keep_indices]
    styles = build_irrep_styles(selected_irreps)

    plt.figure(figsize=(8, 5))

    for irrep in sorted(set(selected_irreps)):
        color, marker = styles[irrep]
        irrep_indices = [i for i in keep_indices if irreps[i] == irrep]
        x_vals = [
            cluster_positions.get((momenta[i], sz_labels[i]), x_positions[momenta[i]])
            + horizontal_offsets[i]
            for i in irrep_indices
        ]
        y_vals = [eigenvalues[i] - energy_shift for i in irrep_indices]
        plt.scatter(
            x_vals,
            y_vals,
            label=irrep,
            color=color,
            marker=marker,
            s=60,
            edgecolors="black",
            linewidths=0.4,
        )

    plt.xticks(tick_positions, tick_labels, rotation=45, ha="right")
    plt.xlabel(r"Momentum with $|Sz - N/2|$ subclusters")
    plt.ylabel(r"Shifted eigenvalue  $E - E_{0,\Gamma}$")
    plt.title("Lowest eigenvalues per momentum sector")
    plt.grid(axis="y", linestyle="--", alpha=0.3)
    for xpos in x_positions.values():
        plt.axvline(x=xpos, color="gray", linestyle=":", linewidth=0.6, alpha=0.35)
    plt.legend(title="Irrep", fontsize="small", ncol=2)
    plt.tight_layout()

    if output is not None:
        plt.savefig(output, dpi=300)

    if show:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Plot the lowest eigenvalues from a Heisenberg spectrum file "
            "as a function of momentum with optional Sz subclusters."
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=Path("HeisenbergZeemanSzvar24.outfile.txt"),
        help="Path to the source spectrum .txt file.",
    )
    parser.add_argument(
        "--n-lowest",
        type=int,
        default=5,
        help="Number of lowest eigenvalues kept per momentum sector.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional path where the plot image is written.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open an interactive plot window.",
    )
    parser.add_argument(
        "--degeneracy-tol",
        type=float,
        default=1e-3,
        help="Tolerance used to identify close/degenerate eigenvalues.",
    )
    parser.add_argument(
        "--jitter-width",
        type=float,
        default=0.06,
        help="Maximum symmetric horizontal offset used to separate degenerate points.",
    )
    parser.add_argument(
        "--nsites",
        type=int,
        default=24,
        help="Number of lattice sites N used for |Sz - N/2| axis labels.",
    )
    args = parser.parse_args()

    if args.n_lowest <= 0:
        raise ValueError("--n-lowest must be a positive integer.")
    if args.degeneracy_tol < 0.0:
        raise ValueError("--degeneracy-tol must be non-negative.")
    if args.jitter_width < 0.0:
        raise ValueError("--jitter-width must be non-negative.")
    if args.nsites <= 0:
        raise ValueError("--nsites must be a positive integer.")

    momentum_order = ["Gamma", "K", "M"]

    eigenvalues, momenta, irreps, sz_labels = parse_spectrum_file(args.input)
    keep_indices = select_lowest_indices(
        eigenvalues=eigenvalues,
        momenta=momenta,
        allowed_momenta=momentum_order,
        n_lowest=args.n_lowest,
    )

    plot_lowest_spectrum(
        eigenvalues=eigenvalues,
        momenta=momenta,
        irreps=irreps,
        sz_labels=sz_labels,
        keep_indices=keep_indices,
        momentum_order=momentum_order,
        output=args.output,
        show=not args.no_show,
        degeneracy_tol=args.degeneracy_tol,
        jitter_width=args.jitter_width,
        nsites=args.nsites,
    )


if __name__ == "__main__":
    main()

