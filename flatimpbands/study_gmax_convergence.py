"""
Study convergence of projected Haldane interaction matrix elements as a function of gmax.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from genprojectedinteraction import make_haldane_interaction_from_givebloch

_HERE = os.path.dirname(os.path.abspath(__file__))
BLOCH_FILE = os.path.join(_HERE, "bloch_nk24.npy")

# Test parameters
NK = 24
PSEUDO_COEFFS = {1: 1.0, 3: 0.2}
GMAX_VALUES = [0, 1, 2, 3, 4, 5]
BAND = 0

# Test momentum points
TEST_POINTS = [
    ((0, 0), (1, 2), (3, 1)),
    ((5, 5), (7, 8), (10, 10)),
    ((12, 12), (13, 14), (14, 15)),
]


def main():
    print(f"Loading bloch data from {BLOCH_FILE}...")
    
    # Dictionary to store results: {test_point: [gmax_values, matrix_element_values]}
    results = {tp: {"gmax": [], "values": [], "real": [], "imag": []} for tp in TEST_POINTS}
    
    print(f"\nComputing matrix elements for gmax = {GMAX_VALUES}...")
    for gmax in GMAX_VALUES:
        print(f"  gmax = {gmax}...", end=" ", flush=True)
        
        hint = make_haldane_interaction_from_givebloch(
            nkpoint=NK,
            pseudo_coeffs=PSEUDO_COEFFS,
            band=BAND,
            gmax=gmax,
            antisymmetrized=True,
            bloch_file=BLOCH_FILE,
        )
        
        for k1, k2, k3 in TEST_POINTS:
            val = hint(k1, k2, k3)
            results[(k1, k2, k3)]["gmax"].append(gmax)
            results[(k1, k2, k3)]["values"].append(val)
            results[(k1, k2, k3)]["real"].append(val.real)
            results[(k1, k2, k3)]["imag"].append(val.imag)
        
        print("done")
    
    # Print results in tabular form
    print("\n" + "="*80)
    print("Projected Haldane Interaction Matrix Elements vs gmax")
    print("="*80)
    for k1, k2, k3 in TEST_POINTS:
        print(f"\nH_AS[{k1}, {k2}, {k3}]:")
        print(f"{'gmax':<6} {'Real Part':<18} {'Imag Part':<18} {'Magnitude':<18}")
        print("-" * 64)
        for gmax, val in zip(results[(k1, k2, k3)]["gmax"], results[(k1, k2, k3)]["values"]):
            mag = abs(val)
            print(f"{gmax:<6} {val.real:>17.12f} {val.imag:>17.12f} {mag:>17.12f}")
    
    # Plot convergence
    fig, axes = plt.subplots(1, len(TEST_POINTS), figsize=(15, 5))
    if len(TEST_POINTS) == 1:
        axes = [axes]
    
    for ax, (k1, k2, k3) in zip(axes, TEST_POINTS):
        gmax_vals = results[(k1, k2, k3)]["gmax"]
        real_vals = results[(k1, k2, k3)]["real"]
        imag_vals = results[(k1, k2, k3)]["imag"]
        
        ax.plot(gmax_vals, real_vals, 'o-', label='Real part', linewidth=2, markersize=8)
        ax.plot(gmax_vals, imag_vals, 's-', label='Imag part', linewidth=2, markersize=8)
        ax.set_xlabel('gmax', fontsize=12)
        ax.set_ylabel('Matrix element', fontsize=12)
        ax.set_title(f'H_AS{k1, k2, k3}', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_xticks(gmax_vals)
    
    plt.tight_layout()
    outfile = os.path.join(_HERE, "gmax_convergence.png")
    plt.savefig(outfile, dpi=150)
    print(f"\nPlot saved to {outfile}")
    
    # Summary
    print("\n" + "="*80)
    print("Convergence Summary:")
    print("="*80)
    for k1, k2, k3 in TEST_POINTS:
        vals = results[(k1, k2, k3)]["values"]
        if len(vals) > 1:
            last_val = vals[-1]
            prev_val = vals[-2]
            diff = abs(last_val - prev_val)
            print(f"H_AS{k1, k2, k3}:")
            print(f"  gmax={GMAX_VALUES[-2]} -> gmax={GMAX_VALUES[-1]}: Δ = {diff:.2e}")
            print(f"  Final value: {last_val:.12f}")


if __name__ == "__main__":
    main()

