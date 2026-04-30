"""Check whether projected Haldane interaction depends only on pair momentum sums.

The tested hypothesis is:
  H(k1, k2; k3, k4) depends only on (k1 + k3) and (k2 + k4) modulo Nk.

For each random momentum-conserving tuple (k1, k2, k3, k4), we generate shifted tuples
that keep these pair sums fixed:
  k1' = k1 + q,  k3' = k3 - q,  k2' = k2 - q,  k4' = k4 + q  (mod Nk)

If the hypothesis is true, H should be unchanged under these shifts.
"""

import argparse
import os

import numpy as np

from genprojectedinteraction import make_haldane_interaction_from_givebloch


def _add_pair(a, b, nk):
    return ((a[0] + b[0]) % nk, (a[1] + b[1]) % nk)


def _sub_pair(a, b, nk):
    return ((a[0] - b[0]) % nk, (a[1] - b[1]) % nk)


def _infer_k4(k1, k2, k3, nk):
    return ((k1[0] + k2[0] - k3[0]) % nk, (k1[1] + k2[1] - k3[1]) % nk)


def _random_pair(rng, nk):
    return (int(rng.integers(0, nk)), int(rng.integers(0, nk)))


def run_check(hint, nk, nsamples, nshifts, tol):
    rng = np.random.default_rng(12345)

    max_diff = 0.0
    worst = None

    for _ in range(nsamples):
        k1 = _random_pair(rng, nk)
        k2 = _random_pair(rng, nk)
        k3 = _random_pair(rng, nk)
        k4 = _infer_k4(k1, k2, k3, nk)

        base_val = hint(k1, k2, k3, k4)

        for _ in range(nshifts):
            q = _random_pair(rng, nk)
            if q == (0, 0):
                continue

            k1p = _add_pair(k1, q, nk)
            k3p = _sub_pair(k3, q, nk)
            k2p = _sub_pair(k2, q, nk)
            k4p = _add_pair(k4, q, nk)

            # Sanity: transformed tuple keeps both tested pair sums fixed.
            s13 = _add_pair(k1, k3, nk)
            s24 = _add_pair(k2, k4, nk)
            s13p = _add_pair(k1p, k3p, nk)
            s24p = _add_pair(k2p, k4p, nk)
            if s13 != s13p or s24 != s24p:
                raise RuntimeError("Internal error: pair sums were not preserved")

            shifted_val = hint(k1p, k2p, k3p, k4p)
            diff = abs(shifted_val - base_val)
            if diff > max_diff:
                max_diff = diff
                worst = {
                    "k": (k1, k2, k3, k4),
                    "kp": (k1p, k2p, k3p, k4p),
                    "base": base_val,
                    "shifted": shifted_val,
                    "diff": diff,
                }

    passed = max_diff <= tol
    return passed, max_diff, worst


def main():
    parser = argparse.ArgumentParser(description="Test pair-sum dependence of projected Haldane interaction")
    parser.add_argument("--nk", type=int, default=24, help="Nk grid size")
    parser.add_argument("--gmax", type=int, default=3, help="gmax for interaction")
    parser.add_argument("--band", type=int, default=0, help="band index")
    parser.add_argument("--nsamples", type=int, default=30, help="number of base tuples")
    parser.add_argument("--nshifts", type=int, default=8, help="number of shifts per base tuple")
    parser.add_argument("--tol", type=float, default=1e-10, help="absolute tolerance for invariance")
    parser.add_argument(
        "--bloch-file",
        type=str,
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "bloch_nk24.npy"),
        help="path to cached Bloch .npy file",
    )
    args = parser.parse_args()

    hint_direct = make_haldane_interaction_from_givebloch(
        nkpoint=args.nk,
        pseudo_coeffs={1: 1.0, 3: 0.2},
        band=args.band,
        gmax=args.gmax,
        antisymmetrized=False,
        bloch_file=args.bloch_file,
    )
    hint_antisym = make_haldane_interaction_from_givebloch(
        nkpoint=args.nk,
        pseudo_coeffs={1: 1.0, 3: 0.2},
        band=args.band,
        gmax=args.gmax,
        antisymmetrized=True,
        bloch_file=args.bloch_file,
    )

    print(f"Testing pair-sum invariance with nk={args.nk}, gmax={args.gmax}, tol={args.tol:.1e}")

    mode_pass = {}
    for label, hint in (("direct", hint_direct), ("antisymmetrized", hint_antisym)):
        passed, max_diff, worst = run_check(
            hint=hint,
            nk=args.nk,
            nsamples=args.nsamples,
            nshifts=args.nshifts,
            tol=args.tol,
        )
        mode_pass[label] = passed
        print(f"  mode={label:14s} max|Delta| = {max_diff:.6e} -> {'PASS' if passed else 'FAIL'}")
        if not passed and worst is not None:
            print("    worst original tuple:", worst["k"])
            print("    worst shifted  tuple:", worst["kp"])
            print("    base value    =", worst["base"])
            print("    shifted value =", worst["shifted"])

    # Mark script as failed if either mode fails.
    if not (mode_pass.get("direct", False) and mode_pass.get("antisymmetrized", False)):
        raise SystemExit(1)


if __name__ == "__main__":
    main()


