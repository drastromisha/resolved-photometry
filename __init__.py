"""resphot: resolved multi-band photometry for galaxy FITS maps.

Aperture and per-pixel photometry across an arbitrary set of imaging bands,
with unit conversion and two error estimators (per-pixel quadrature and the
empty-aperture method). Instrument setup is driven entirely by a YAML config,
so the package is not tied to any particular dataset.
"""
__version__ = "0.1.0"

from .photometry import aperture_sum, error_quadrature, error_empty_aperture
from .units import to_mjy_per_pixel
