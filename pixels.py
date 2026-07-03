"""Per-pixel multi-band photometry.

For every pixel in a common grid, assemble a flux value per band (converted to
mJy/pixel) and a per-band per-pixel error taken as the sigma-clipped RMS of the
noise map in a source-free annulus. The RMS is estimated on the noise map you
supply; if that is the ORIGINAL (non-PSF-matched) map, the resulting per-pixel
error is not inflated by PSF-induced correlation.

All maps must already share the same shape and pixel grid.
"""
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clip
from astropy.table import Table

from .units import to_mjy_per_pixel


def _annulus_mask(shape, cx, cy, r_inner, r_outer):
    yy, xx = np.indices(shape)
    r2 = (xx - cx) ** 2 + (yy - cy) ** 2
    return (r2 >= r_inner ** 2) & (r2 <= r_outer ** 2)


def _band_rms(noise_mjy, mask):
    vals = noise_mjy[mask]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan
    return float(np.nanstd(sigma_clip(vals, sigma=3.0, maxiters=5).filled(np.nan)))


def extract_pixel_table(cfg, annulus_pix=(25, 30), peak_band=None):
    """Build a per-pixel photometry table from a validated config.

    Parameters
    ----------
    cfg : dict            validated config (see config.load_config)
    annulus_pix : (int, int)  inner, outer radii (pixels) for the noise annulus
    peak_band : str       label of the band whose peak defines the annulus centre;
                          defaults to the first band.

    Returns
    -------
    astropy.table.Table   one row per pixel, columns <label> and <label>_err.
    """
    bands = cfg["bands"]
    pixscale = cfg["pixscale_arcsec"]
    redshift = cfg.get("redshift")

    # Reference frame from the first band
    with fits.open(bands[0]["flux"]) as h:
        ref = np.array(h[0].data, dtype=float)
        ref_hdr = h[0].header
    ny, nx = ref.shape
    npix = ny * nx
    yy, xx = np.indices((ny, nx))

    # Annulus centre = peak of chosen band
    peak = bands[0]
    if peak_band is not None:
        peak = next((b for b in bands if b["label"] == peak_band), bands[0])
    with fits.open(peak["flux"]) as h:
        peak_data = np.array(h[0].data, dtype=float)
    iy_pk, ix_pk = np.unravel_index(np.nanargmax(peak_data), peak_data.shape)
    bg_mask = _annulus_mask((ny, nx), ix_pk, iy_pk, *annulus_pix)

    # Optional WCS
    try:
        w2d = WCS(ref_hdr).celestial
        world = w2d.pixel_to_world(xx.ravel(), yy.ravel())
        ra = np.array(world.ra.deg).ravel()
        dec = np.array(world.dec.deg).ravel()
    except Exception:
        ra = np.full(npix, np.nan)
        dec = np.full(npix, np.nan)

    df = pd.DataFrame({
        "id": [f"px_{ix:03d}_{iy:03d}" for iy in range(ny) for ix in range(nx)],
        "x": xx.ravel().astype(int),
        "y": yy.ravel().astype(int),
        "ra_deg": ra,
        "dec_deg": dec,
    })
    if redshift is not None:
        df["redshift"] = np.full(npix, redshift)

    for b in bands:
        with fits.open(b["flux"]) as h:
            flux = np.array(h[0].data, dtype=float)
        with fits.open(b["noise"]) as h:
            noise = np.array(h[0].data, dtype=float)
        if flux.shape != (ny, nx) or noise.shape != (ny, nx):
            raise RuntimeError(f"Band '{b['label']}' shape mismatch with reference grid.")

        flux_mjy = to_mjy_per_pixel(flux, b["unit"], pixscale, b.get("zeropoint_ab"))
        noise_mjy = to_mjy_per_pixel(noise, b["unit"], pixscale, b.get("zeropoint_ab"))

        rms = _band_rms(noise_mjy, bg_mask)
        df[b["label"]] = flux_mjy.ravel()
        df[f"{b['label']}_err"] = np.full(npix, rms)

    t = Table.from_pandas(df)
    t.meta["NX"], t.meta["NY"] = nx, ny
    t.meta["PIXSCALE"] = pixscale
    if redshift is not None:
        t.meta["REDSHIFT"] = redshift
    return t
