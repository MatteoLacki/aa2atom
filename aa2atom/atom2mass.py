from .isotopes import monoisotopes, meanisotopes
from .formulas import chem2atom


def atom2mass(formula, which_mass='monoisotopic'):
    """Get the molecular mass of a chemical formula.

    For a full isotopic envelope, use the IsoSpecPy module.

    Args:
        formula (dict or str): A mapping between atoms and their number of occurences.
        which (str): Either 'monoisotopic' or 'average'.
    Returns:
        float: the calculated mass.
    """
    assert which_mass in ('monoisotopic', 'average')
    res = 0.0
    iso = monoisotopes if which_mass == 'monoisotopic' else meanisotopes
    try:
        formula = chem2atom(formula)
    except TypeError:
        pass
    mass = sum(iso[elem]*atomcnt for elem, atomcnt in formula.items())
    return mass

