# Heisenberg Spectrum Plotter

This helper script reads `HeisenbergZeemanSzvar24.outfile.txt`, stores the data in four arrays,
and plots only the lowest eigenvalues per momentum sector.

Expected line format:

`eigenvalue,Sz,momentum.irrep`

Example:

`-15.1737,8,Gamma.C6.A`

The plotted energies are shifted by subtracting the minimum `Gamma` energy,
so the lowest `Gamma` level is at `0`. Near-degenerate points at the same momentum
are separated by a small horizontal jitter for readability.

## Parsed arrays

- `eigenvalues`: floating-point eigenvalues
- `sz_labels`: Sz labels (for example `8`, `6`, `4`, ...)
- `momenta`: momentum labels (`Gamma`, `K`, `M`)
- `irreps`: irrep labels (for example `C6.A`, `C3.Ea`)

The horizontal axis is grouped by momentum, with subclusters labeled by
`|Sz - Nsites/2|` around each momentum point. The `Nsites` value is controlled
with `--nsites` (default: `24`).

## Run

```bash
cd /home/t30/all/ge45hub/CLionProjects/xdiagAFM/cantedAFM
python3 plot_heisenberg_spectrum.py HeisenbergZeemanSzvar24.outfile.txt --n-lowest 5
```

## Headless run (save to file)

```bash
cd /home/t30/all/ge45hub/CLionProjects/xdiagAFM/cantedAFM
python3 plot_heisenberg_spectrum.py HeisenbergZeemanSzvar24.outfile.txt --n-lowest 5 --output lowest5.png --no-show
```

Optional tuning:

```bash
cd /home/t30/all/ge45hub/CLionProjects/xdiagAFM/cantedAFM
python3 plot_heisenberg_spectrum.py HeisenbergZeemanSzvar24.outfile.txt --n-lowest 5 --degeneracy-tol 1e-3 --jitter-width 0.08
```

Example with explicit system size:

```bash
cd /home/t30/all/ge45hub/CLionProjects/xdiagAFM/cantedAFM
python3 plot_heisenberg_spectrum.py HeisenbergZeemanSzvar24.outfile.txt --n-lowest 5 --nsites 24
```

