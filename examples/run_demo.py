"""End-to-end demo: generate a synthetic galaxy, then run aperture and
per-pixel photometry on it. Produces a table and a diagnostic figure using
only data created at runtime — no external files needed.

Run from the repo root:
    python examples/run_demo.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from matplotlib.colors import LogNorm
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from astropy.io import fits

from resphot.mockdata import make_mock_galaxy
from resphot.photometry import aperture_sum, error_quadrature, error_empty_aperture
from resphot.pixels import extract_pixel_table
from resphot.config import load_config, resolve_paths

HERE = os.path.dirname(__file__)
MOCK_DIR = os.path.join(HERE, "mock_data")
OUT_DIR = os.path.join(HERE, "demo_output")
os.makedirs(OUT_DIR, exist_ok=True)

bands = ["b1", "b2", "b3", "b4", "b5"]
pixscale = 0.06
r_core_arcsec, r_int_arcsec = 0.14, 1.40

# 1. Synthetic data
print("Generating synthetic galaxy...")
paths = make_mock_galaxy(MOCK_DIR, bands=bands, shape=(61, 61), seed=0)

# 2. Aperture photometry on one band
band0 = bands[0]
with fits.open(paths[band0]["flux"]) as h:
    flux = np.array(h[0].data, dtype=float)
with fits.open(paths[band0]["noise"]) as h:
    noise = np.array(h[0].data, dtype=float)

iy, ix = np.unravel_index(np.nanargmax(flux), flux.shape)
r_core = r_core_arcsec / pixscale
r_int = r_int_arcsec / pixscale

f_core, m_core, _ = aperture_sum(flux, ix, iy, r_core)
f_int, m_int, _ = aperture_sum(flux, ix, iy, r_int)
e_core = error_quadrature(noise, m_core)
e_int = error_quadrature(noise, m_int)
e_core_ea, _, _, _ = error_empty_aperture(noise, ix, iy, r_core)

print(f"  core: {f_core:.3f} +/- {e_core:.3f} mJy "
      f"(empty-aperture err {e_core_ea:.3f})")
print(f"  integrated: {f_int:.3f} +/- {e_int:.3f} mJy")

# 3. Per-pixel table via config
cfg = load_config(os.path.join(HERE, "config_example.yaml"))
cfg = resolve_paths(cfg, base_dir=HERE)
table = extract_pixel_table(cfg, annulus_pix=(25, 30), peak_band=band0)
table_path = os.path.join(OUT_DIR, "pixel_fluxes.fits")
table.write(table_path, overwrite=True)
print(f"  per-pixel table: {table_path} ({len(table)} pixels)")

# 4. Diagnostic figure
fig, ax = plt.subplots(figsize=(6, 5))
pos = flux[flux > 0]
vmin = np.percentile(pos, 50) if pos.size else 1e-3
vmax = np.percentile(flux, 99.5)
im = ax.imshow(flux, origin="lower", cmap="inferno",
               norm=LogNorm(vmin=max(vmin, 1e-3), vmax=vmax))
ax.add_patch(Circle((ix, iy), r_core, fill=False, ec="cyan", lw=2, label="core"))
ax.add_patch(Circle((ix, iy), r_int, fill=False, ec="lime", lw=2, ls="--",
                    label="integrated"))
ax.plot(ix, iy, "w+", ms=12, mew=2)
ax.legend(fontsize=8)
ax.set_title("Synthetic galaxy with apertures")
plt.colorbar(im, ax=ax, label="mJy/pixel")
fig.savefig(os.path.join(OUT_DIR, "demo_apertures.png"), dpi=130,
            bbox_inches="tight")
print(f"  figure: {os.path.join(OUT_DIR, 'demo_apertures.png')}")
print("Done.")
