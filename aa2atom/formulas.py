import re
from collections import Counter

chem_formula = re.compile('([A-Z][a-z]?)([-]?[0-9]*)')


def chem2atom(chem):
    """Transform a string with chemical formula to a dictionary with atom counts.

    Args:
        chem (str): A string with a chemical formula, or a chemical diff, e.g. 'C10HNa', 'C-1H-2O2'.
    Return:
        dict: counts of atoms.
    """
    global chem_formula
    return Counter({el: 1 if cnt is '' else int(cnt)
    				for el, cnt in chem_formula.findall(chem)})

def test_chem2atom():
	assert chem2atom('CH2Cu-1') == {'C':1,'H':2,'Cu':-1}

def add_formulas(f, g, check_positivity=True):
	res = Counter()
	for a in set(f.keys()) | set(g.keys()):
		cnt = f.get(a, 0) + g.get(a, 0)
		assert cnt >= 0, "Adding formulas resulted in a negative atom count. This might explain the dark matter!"
		if cnt:
			res[a] = cnt
	return res

def test_add_formulas():
	assert add_formulas({'C':10, 'H':20}, {'C':2, 'H':2}) == {'H':22, 'C':12}
	assert add_formulas({'C':10, 'H':20}, {'C':-10, 'H':-18}) == {'H':2}