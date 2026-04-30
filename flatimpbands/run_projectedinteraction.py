import os
import numpy as np

from genprojectedinteraction import make_haldane_interaction_from_givebloch

_HERE = os.path.dirname(os.path.abspath(__file__))
BLOCH_FILE = os.path.join(_HERE, "bloch_nk24.npy")


def main():
    nk = 24
    pseudo_coeffs = {1: 1.0, 3: 0.2}
    hint = make_haldane_interaction_from_givebloch(
        nkpoint=nk,
        pseudo_coeffs=pseudo_coeffs,
        band=0,
        gmax=0,
        antisymmetrized=True,
        bloch_file=BLOCH_FILE,
    )

    k1 = (0, 0)
    k2 = (1, 2)
    k3 = (3, 1)
    val = hint(k1, k2, k3)

    print(f"H_AS{k1, k2, k3} = {val}")
    print(f"l_B^2 = {hint.lb2:.6f}")


if __name__ == "__main__":
    main()

