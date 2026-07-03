# resolved-photometry (resphot)

Resolved aperture and per-pixel multi-band photometry for galaxy imaging.

`resphot` measures fluxes from a set of co-registered, multi-band FITS maps of a
galaxy — both in circular apertures and pixel-by-pixel — with unit conversion
and two error estimators. The instrument setup is described in a YAML config, so
the package is not tied to any particular dataset or object.

## Features

- Aperture photometry with two error estimators:
  - **per-pixel quadrature** (for uncorrelated, RMS noise maps), and
  - the **empty-aperture** method (for small apertures on a wide field).
- Per-pixel extraction to a table (one row per pixel, flux + error per band).
- Unit conversion to mJy/pixel from `mjy_per_pixel`, `mjy_per_sr`, or HST `hst_eps`
  (via an AB zeropoint).
- A synthetic-galaxy generator so the demo runs with no external data.

## Install

```bash
git clone https://github.com/<your-username>/resolved-photometry.git
cd resolved-photometry
pip install -r requirements.txt
```

## Quick start

```bash
python examples/run_demo.py
```

This generates a synthetic multi-band galaxy, runs aperture and per-pixel
photometry, writes a per-pixel table, and saves a diagnostic figure to
`examples/demo_output/`.

## Using your own data

Copy `examples/config_example.yaml`, point each band at its flux and noise FITS
files, set the pixel scale and (optionally) redshift, and choose the unit per
band. All maps must share the same pixel grid. Then:

```python
from resphot.config import load_config, resolve_paths
from resphot.pixels import extract_pixel_table

cfg = load_config("my_config.yaml")
cfg = resolve_paths(cfg, base_dir=".")
table = extract_pixel_table(cfg)
table.write("pixel_fluxes.fits", overwrite=True)
```

## Notes

- Supply the **original** (non-PSF-matched) noise maps for error estimation;
  PSF matching correlates pixel noise and biases naive per-pixel errors.
- The synthetic data are illustrative only and do not model any real object.

## License

MIT — see `LICENSE`.

## Citation

If this software is useful in your work, please cite it via the archived
release (see the repository's Zenodo DOI).
