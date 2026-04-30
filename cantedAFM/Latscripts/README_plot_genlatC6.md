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
python3 plot_genlatC6.py real 3 0 --no-show-nnn-bonds
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
- Real-space plots show NN and NNN bonds by default; use `--no-show-bonds` to hide all bonds or `--no-show-nnn-bonds` to keep only NN bonds.
- NNN bonds are drawn as directed dashed arrows, so the displayed arrow direction shows the order of each index pair `(i, j)`.
- Bonds that cross the torus fundamental domain are drawn to a translated periodic image of the target site, so wrapped bonds remain short and local.
- Use `--show-rotation-labels` to overlay the lookup index and its rotated index in different colors.
- Use `--display-coordinates wigner-seitz` for a plotting-only representative choice adapted to the torus Wigner-Seitz cell; site markers stay in those representatives while bonds may extend to translated periodic images.
- Momentum plots highlight and label high-symmetry points `Gamma`, `K`, `Kp`, and `M`.
