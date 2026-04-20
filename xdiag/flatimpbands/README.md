# Flat Impurity Bands: Projected Haldane Interaction

This folder contains tools to generate Bloch wavefunctions and evaluate projected two-body interaction matrix elements on a momentum grid.

## What is implemented

- `givebloch(Nkpoint, return_metadata=False)` returns Bloch states on the Kwant wraparound momentum grid.
- `make_projected_interaction_callable(...)` builds a lazy callable for matrix elements
  `H(k1, k2, k3[, k4])`, with momentum inputs as grid index pairs `(ix, iy)`.
- `make_haldane_interaction_from_givebloch(...)` is a convenience constructor using `givebloch`.

By default, returned values are antisymmetrized for fermions:

`H_AS(k1,k2;k3,k4) = H(k1,k2;k3,k4) - H(k1,k2;k4,k3)`.

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

