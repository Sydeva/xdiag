# Plot Structure Factor (Lorentzian)

This helper plots a Lorentzian-broadened spectrum from long-form rows in:

- `momentum,energy,weight`

## Script

- `cantedAFM/strucfac/plot_strucfac_lorentzian.py`

## Defaults

- Momentum selector: `Gamma`
- Broadening: `eps = 0.04`
- Maximum frequency: `2.0`
- Input path (if omitted): `teststrucfac.outfile.txt`

## Example

```bash
python3 cantedAFM/strucfac/plot_strucfac_lorentzian.py \
  cantedAFM/strucfac/teststrucfac.outfile.txt \
  --momentum Gamma \
  --eps 0.04 \
  --max-frequency 2.0
```

Use an exact label if needed:

```bash
python3 cantedAFM/strucfac/plot_strucfac_lorentzian.py \
  cantedAFM/strucfac/teststrucfac.outfile.txt \
  --momentum Gamma.C1.A \
  --eps 0.04
```

Headless save:

```bash
python3 cantedAFM/strucfac/plot_strucfac_lorentzian.py \
  cantedAFM/strucfac/teststrucfac.outfile.txt \
  --momentum M1.C1.A \
  --eps 0.06 \
  --output /tmp/m1_spectrum.png \
  --no-show
```

## Smoke test

```bash
python3 cantedAFM/strucfac/test_plot_strucfac_lorentzian_smoke.py
```
