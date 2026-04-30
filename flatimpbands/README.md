# Flat Impurity Bands: Projected Haldane Interaction

This folder contains tools to generate Bloch wavefunctions and evaluate projected two-body interaction matrix elements on a momentum grid. On 64 sites max 6 electrons

## What is implemented

- `givebloch(Nkpoint, return_metadata=False)` returns Bloch states on the Kwant wraparound momentum grid.
- `make_projected_interaction_callable(...)` builds a lazy callable for matrix elements
  `H(k1, k2, k3[, k4])`, with momentum inputs as grid index pairs `(ix, iy)`.
- `make_haldane_interaction_from_givebloch(...)` is a convenience constructor using `givebloch`.

By default, returned values are antisymmetrized for fermions in both pairs:

`H_AS(k1,k2;k3,k4) = H(k1,k2;k3,k4) - H(k2,k1;k3,k4) - H(k1,k2;k4,k3) + H(k2,k1;k4,k3)`.

## Notes on conventions

- Momentum inputs are integer indices on the `Nkpoint x Nkpoint` grid.
- If `k4` is omitted, momentum conservation is used:
  `k4 = k1 + k2 - k3 (mod Nkpoint)`.
- Orbital positions are taken from `dxs`.
- The Haldane pseudopotential kernel is configured by `pseudo_coeffs = {m: V_m}`.

## Quick run

```bash
python -u xdiag/flatimpbands/run_projectedinteraction.py
```

## Dependencies

- `numpy`
- `scipy`
- `kwant`
- `opt_einsum`

## Square Lattice TOML Generator

- `genlatsquare.py` generates square-lattice TOML files with **translation symmetry only**.
- Inputs are rectangular: `Lx`, `Ly`.
- Momentum sections are labeled numerically: `k1`, `k2`, ...
- Irrep sections include `characters` and `momentum` fields (no `allowed_symmetries`).
- The generated TOML intentionally omits `Interactions`.

Quick run:

```bash
python -u flatimpbands/genlatsquare.py 4 3 --toml flatimpbands/square.12.T.toml --verify
```

Smoke test:

```bash
python -u flatimpbands/test_genlatsquare_smoke.py
```
