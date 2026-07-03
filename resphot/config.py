"""YAML configuration loading and validation.

A config describes the dataset generically: global settings (pixel scale,
redshift, aperture radii, output paths) and a list of bands, each with a flux
file, a noise file, a unit, and any unit-specific parameter (e.g. an AB
zeropoint). Nothing here is tied to a particular instrument or object.
"""
import os
import yaml

REQUIRED_TOP = ["pixscale_arcsec", "bands"]
VALID_UNITS = {"mjy_per_pixel", "mjy_per_sr", "hst_eps"}


def load_config(path):
    """Load and validate a YAML config file. Returns a dict."""
    with open(path, "r") as f:
        cfg = yaml.safe_load(f)

    for key in REQUIRED_TOP:
        if key not in cfg:
            raise ValueError(f"Config missing required key: '{key}'")

    if not cfg["bands"]:
        raise ValueError("Config 'bands' list is empty.")

    for i, b in enumerate(cfg["bands"]):
        for key in ("label", "flux", "noise", "unit"):
            if key not in b:
                raise ValueError(f"Band #{i} missing '{key}'")
        if b["unit"] not in VALID_UNITS:
            raise ValueError(
                f"Band '{b['label']}' has invalid unit '{b['unit']}'. "
                f"Valid: {sorted(VALID_UNITS)}")
        if b["unit"] == "hst_eps" and "zeropoint_ab" not in b:
            raise ValueError(f"Band '{b['label']}' (hst_eps) needs 'zeropoint_ab'")

    # Sensible defaults
    cfg.setdefault("redshift", None)
    cfg.setdefault("aperture_radii_arcsec", {})
    cfg.setdefault("output_dir", "resphot_output")
    return cfg


def resolve_paths(cfg, base_dir=""):
    """Prefix relative band file paths with base_dir (e.g. config's own folder)."""
    for b in cfg["bands"]:
        if not os.path.isabs(b["flux"]):
            b["flux"] = os.path.join(base_dir, b["flux"])
        if not os.path.isabs(b["noise"]):
            b["noise"] = os.path.join(base_dir, b["noise"])
    return cfg
