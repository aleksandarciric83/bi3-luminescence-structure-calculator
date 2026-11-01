"""
bi3_validate.py
===============
Checks for uploaded structure files. The calculation engine accepts almost any text without an
error and fills in defaults (cell a = 5 Angstrom, Z = 1, an empty site list), which would show
wrong numbers without a warning. These functions stop such input before it reaches the engine.
The engine itself is not changed.
"""
import re

CELL_LENGTH_TAGS = ("_cell_length_a", "_cell_length_b", "_cell_length_c")
CELL_ANGLE_TAGS = ("_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")


def _has_tag(text, tag):
    return re.search(re.escape(tag) + r"\s+[^\s#]+", text) is not None


def check_cif_text(text):
    """Return (errors, warnings) for the text of a CIF file, before parsing."""
    errors, warnings = [], []
    for tag in CELL_LENGTH_TAGS:
        if not _has_tag(text, tag):
            errors.append("the cell length %s is missing" % tag)
    if "_atom_site_fract_x" not in text:
        errors.append("no atom sites with fractional coordinates (_atom_site_fract_x) were found")
    if not _has_tag(text, "_cell_formula_units_Z"):
        warnings.append("_cell_formula_units_Z is missing, so Z = 1 was assumed. Z changes the "
                        "stoichiometric coefficients; check that Z = 1 is correct for this structure.")
    missing_angles = [t for t in CELL_ANGLE_TAGS if not _has_tag(text, t)]
    if missing_angles:
        warnings.append("cell angles missing (%s); 90 degrees were assumed." % ", ".join(missing_angles))
    return errors, warnings


def check_json(raw):
    """Return a list of errors for a parsed JSON structure file."""
    if not isinstance(raw, dict):
        return ["the JSON file must contain one object, not a list or a single value"]
    cs = raw.get("crystal_structure", raw)
    positions = cs.get("atomic_positions") if isinstance(cs, dict) else None
    pymatgen_like = bool(raw.get("lattice")) and bool(raw.get("sites"))
    if not positions and not pymatgen_like:
        return ["no atom positions were found (expected 'atomic_positions' inside 'crystal_structure', "
                "or the pymatgen layout with 'lattice' and 'sites')"]
    return []


def check_crystal(crystal):
    """Return a list of errors for a structure that the engine has built."""
    errors = []
    if not crystal.sites_cations:
        errors.append("no cation sites could be built from the file")
    if not crystal.sites_anions:
        errors.append("no anion sites could be built from the file")
    if not (crystal.V and crystal.V > 0):
        errors.append("the unit cell volume is not positive")
    return errors
