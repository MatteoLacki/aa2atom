import re


chem_formula = re.compile('([A-Z][a-z]?)([-]?[0-9]*)')


def chem2atom(chem):
    """Transform a string with chemical formula to a dictionary with atom counts.

    Args:
        chem (str): A string with a chemical formula, or a chemical diff, e.g. 'C10HNa', 'C-1H-2O2'.
    Return:
        dict: counts of atoms.
    """
    global chem_formula
    return {el: 1 if cnt is '' else int(cnt) for el, cnt in chem_formula.findall(chem)}


def test_chem2atom():
	assert chem2atom('CH2Cu-1') == {'C': 1,'H':2,'Cu':-1}