"""Aperture photometry and error estimation.

Two error estimators, matching common practice for resolved photometry:

  * error_quadrature   : sigma_ap = sqrt(sum sigma_i^2) over aperture pixels.
                         Valid when the noise map is a per-pixel 1-sigma (RMS)
                         map and pixel noise is uncorrelated (i.e. the ORIGINAL,
                         non-PSF-matched noise map).

  * error_empty_aperture : scatter of aperture sums placed at random source-free
                         positions on the noise map. Reliable only when the
                         aperture is small relative to the field of view.
"""
import numpy as np
from astropy.stats import sigma_clip


def aperture_sum(data, cx, cy, r_pix):
    """Sum pixels within a circular aperture.

    Returns
    -------
    (flux, mask, npix) : float, ndarray(bool), int
    """
    ny, nx = data.shape
    yy, xx = np.indices((ny, nx))
    mask = (xx - cx) ** 2 + (yy - cy) ** 2 <= r_pix ** 2
    return float(np.nansum(data[mask])), mask, int(np.sum(mask))


def error_quadrature(noise_map, mask):
    """Quadrature sum of per-pixel 1-sigma noise within a mask."""
    vals = noise_map[mask]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return np.nan
    return float(np.sqrt(np.sum(vals ** 2)))


def error_empty_aperture(noise_map, cx_src, cy_src, r_ap_pix,
                         exclusion_r_pix=None, n_target=2000, seed=42):
    """Empty-aperture error: scatter of random source-free aperture sums.

    Parameters
    ----------
    noise_map : ndarray
    cx_src, cy_src : float   source position to avoid
    r_ap_pix : float         aperture radius (pixels)
    exclusion_r_pix : float  radius around the source to exclude; defaults to 3x aperture
    n_target : int           number of empty apertures to place
    seed : int

    Returns
    -------
    (rms, median, n_placed, sums) : float, float, int, ndarray
    """
    ny, nx = noise_map.shape
    if exclusion_r_pix is None:
        exclusion_r_pix = 3.0 * r_ap_pix

    yy, xx = np.indices((ny, nx))
    rng = np.random.default_rng(seed)

    sums = []
    attempts = 0
    max_att = n_target * 500
    while len(sums) < n_target and attempts < max_att:
        attempts += 1
        cx = rng.uniform(r_ap_pix, nx - r_ap_pix - 1)
        cy = rng.uniform(r_ap_pix, ny - r_ap_pix - 1)
        if exclusion_r_pix > 0 and (cx - cx_src) ** 2 + (cy - cy_src) ** 2 < exclusion_r_pix ** 2:
            continue
        m = (xx - cx) ** 2 + (yy - cy) ** 2 <= r_ap_pix ** 2
        vals = noise_map[m]
        if vals.size == 0 or not np.all(np.isfinite(vals)):
            continue
        sums.append(float(np.sum(vals)))

    sums = np.array(sums)
    if sums.size < 10:
        return np.nan, np.nan, len(sums), sums

    clean = sigma_clip(sums, sigma=3.0, maxiters=5).filled(np.nan)
    return float(np.nanstd(clean)), float(np.nanmedian(clean)), len(sums), sums
