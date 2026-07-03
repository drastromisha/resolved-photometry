"""Unit conversions to mJy per pixel.

Each band declares a `unit` in the config. Supported units:
  - "mjy_per_pixel" : already mJy/pixel, no conversion.
  - "mjy_per_sr"    : surface brightness in MJy/sr -> mJy/pixel (needs pixscale).
  - "hst_eps"       : HST e-/s -> mJy/pixel via an AB zeropoint (needs zeropoint).
"""
import numpy as np

# arcsec^2 per steradian conversion constant
_ARCSEC_PER_RAD = 206265.0


def mjy_sr_to_mjy_pixel_factor(pixscale_arcsec):
    """Multiplicative factor: (MJy/sr) * factor -> mJy/pixel."""
    return 1e9 * (pixscale_arcsec ** 2) / (_ARCSEC_PER_RAD ** 2)


def hst_eps_to_mjy_factor(zeropoint_ab):
    """Multiplicative factor: (HST e-/s) * factor -> mJy/pixel, via AB zeropoint."""
    return 10.0 ** ((23.9 - zeropoint_ab) / 2.5) / 1000.0


def to_mjy_per_pixel(data, unit, pixscale_arcsec=None, zeropoint_ab=None):
    """Convert a map to mJy/pixel according to its declared unit.

    Parameters
    ----------
    data : ndarray
    unit : {"mjy_per_pixel", "mjy_per_sr", "hst_eps"}
    pixscale_arcsec : float, required for "mjy_per_sr"
    zeropoint_ab : float, required for "hst_eps"

    Returns
    -------
    ndarray : converted copy, same shape.
    """
    unit = (unit or "").strip().lower()

    if unit == "mjy_per_pixel":
        return np.array(data, dtype=float)

    if unit == "mjy_per_sr":
        if pixscale_arcsec is None:
            raise ValueError("unit 'mjy_per_sr' requires pixscale_arcsec")
        return np.array(data, dtype=float) * mjy_sr_to_mjy_pixel_factor(pixscale_arcsec)

    if unit == "hst_eps":
        if zeropoint_ab is None:
            raise ValueError("unit 'hst_eps' requires zeropoint_ab")
        return np.array(data, dtype=float) * hst_eps_to_mjy_factor(zeropoint_ab)

    raise ValueError(f"Unknown unit '{unit}'. "
                     "Use 'mjy_per_pixel', 'mjy_per_sr', or 'hst_eps'.")
