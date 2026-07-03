"""Synthetic multi-band galaxy generator for demos and tests.

Creates a set of matched flux maps (a smooth Sersic-like disk plus a compact
nuclear component) and matching per-pixel noise maps, written as FITS pairs.
Lets anyone run the pipeline with no external data. Purely illustrative — the
numbers are arbitrary, not a physical model of any real object.
"""
import os
import numpy as np
from astropy.io import fits


def _sersic(ny, nx, cx, cy, amp, r_eff, n=1.0, ellip=0.0, theta=0.0):
    yy, xx = np.indices((ny, nx))
    dx, dy = xx - cx, yy - cy
    ct, st = np.cos(theta), np.sin(theta)
    xr = dx * ct + dy * st
    yr = -dx * st + dy * ct
    q = 1.0 - ellip
    r = np.sqrt(xr ** 2 + (yr / max(q, 1e-3)) ** 2)
    b_n = 2.0 * n - 1.0 / 3.0
    return amp * np.exp(-b_n * ((r / r_eff) ** (1.0 / n) - 1.0))


def make_mock_galaxy(outdir, bands=None, shape=(61, 61), seed=0):
    """Write matched flux maps and original-style noise maps for several bands.

    Parameters
    ----------
    outdir : str          directory to create the FITS files in
    bands : list[str]     band labels; defaults to a small generic set
    shape : (ny, nx)
    seed : int

    Returns
    -------
    dict : {band: {"flux": path, "noise": path}}
    """
    if bands is None:
        bands = ["b1", "b2", "b3", "b4", "b5"]

    os.makedirs(outdir, exist_ok=True)
    ny, nx = shape
    cx, cy = nx // 2, ny // 2
    rng = np.random.default_rng(seed)

    paths = {}
    for i, band in enumerate(bands):
        # Disk brightness rises then falls across bands; nucleus strongest in redder bands
        disk_amp = 1.0 + 0.5 * np.sin(i)
        nuc_amp = 0.3 * (i + 1)

        disk = _sersic(ny, nx, cx, cy, disk_amp, r_eff=8.0, n=1.0, ellip=0.3, theta=0.6)
        nucleus = _sersic(ny, nx, cx, cy, nuc_amp, r_eff=1.5, n=4.0)
        truth = disk + nucleus

        noise_sigma = 0.05 * (1.0 + 0.3 * i)
        flux = truth + rng.normal(0.0, noise_sigma, size=shape)
        noise = np.full(shape, noise_sigma)

        fpath = os.path.join(outdir, f"mock_{band}_flux.fits")
        npath = os.path.join(outdir, f"mock_{band}_noise.fits")
        hdr = fits.Header()
        hdr["BUNIT"] = "mJy/pixel"
        hdr["BAND"] = band
        fits.writeto(fpath, flux.astype("float32"), hdr, overwrite=True)
        fits.writeto(npath, noise.astype("float32"), hdr, overwrite=True)
        paths[band] = {"flux": fpath, "noise": npath}

    return paths
