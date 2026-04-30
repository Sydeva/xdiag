"""Script to generate and cache Bloch wavefunctions for Nkpoint=24.

Run once:
    python generate_bloch_nk24.py

The result is written to bloch_nk24.npy next to this script.
Afterwards, make_haldane_interaction_from_givebloch can load it via the
bloch_file parameter instead of recomputing:

    hint = make_haldane_interaction_from_givebloch(
        nkpoint=24,
        pseudo_coeffs={1: 1.0},
        bloch_file="bloch_nk24.npy",
    )
"""

import os
import time
import numpy as np

# Allow running from any directory.
_HERE = os.path.dirname(os.path.abspath(__file__))

import sys
sys.path.insert(0, _HERE)

from genprojectedinteraction import save_bloch_data

NKPOINT = 24
OUT_FILE = os.path.join(_HERE, "bloch_nk24.npy")


def main():
    print(f"Generating Bloch data for Nkpoint={NKPOINT}...")
    t0 = time.time()
    save_bloch_data(NKPOINT, OUT_FILE)
    elapsed = time.time() - t0
    size_mb = os.path.getsize(OUT_FILE) / 1024**2
    print(f"Saved to {OUT_FILE}  ({size_mb:.1f} MB, {elapsed:.1f} s)")

    # Quick sanity check: reload and print shape.
    data = np.load(OUT_FILE, allow_pickle=True).item()
    print(f"bloch shape : {data['bloch'].shape}")
    print(f"trans_vecs  :\n{data['trans_vecs']}")
    print(f"dxs shape   : {data['dxs'].shape}")


if __name__ == "__main__":
    main()

