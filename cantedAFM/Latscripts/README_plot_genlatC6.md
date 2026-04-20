# plot_genlatC6.py

Small visualizer for the honeycomb C6 torus produced by `genlatC6.py`.

## Requirements

Install plotting dependencies:

```bash
python3 -m pip install -r /home/sideva/PycharmProjects/xdiag/cantedAFM/requirements-plot.txt
```

## Usage

From `/home/sideva/PycharmProjects/xdiag/cantedAFM`:

```bash
python3 plot_genlatC6.py real 3 0
python3 plot_genlatC6.py momentum 4 2
python3 plot_genlatC6.py real 3 0 --no-show-bonds
python3 plot_genlatC6.py real 3 0 --show-rotation-labels
python3 plot_genlatC6.py real 3 0 --display-coordinates wigner-seitz
```

Save figures without opening a GUI window:

```bash
python3 plot_genlatC6.py real 3 0 --save real.png --no-show
python3 plot_genlatC6.py momentum 4 2 --save momentum.png --no-show
```

## Notes

- `space` must be `real` or `momentum`.
- Real-space plots show NN bonds by default; use `--no-show-bonds` to hide them.
- Use `--show-rotation-labels` to overlay the lookup index and its rotated index in different colors.
- Use `--display-coordinates wigner-seitz` for a plotting-only representative choice adapted to the torus Wigner-Seitz cell.
- Momentum plots highlight and label high-symmetry points `Gamma`, `K`, `Kp`, and `M`.
