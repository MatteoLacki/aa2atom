import re
from collections import Counter

ptms = {
    "phosphorylation":              {"P":1, "H":1,  "O":3},
    "acetylation":                  {"C":2, "H":2,  "O":1},
    "glycation_by_hexose":          {"C":6, "H":10, "O":5},
    "amidation":                    {"H":1, "N":1,  "O":-1},
    "hydroxylation":                {"O":1},
    "methylation":                  {"C":1, "H":2},
    "ubiquitylation":               {"C":378, "H":627, "N":105, "O":117, "S":1},
    "pyroglutamic_acid":            {"H":-2, "O":-1},
    "sulfation":                    {"S":1, "O":3},
    "gamma_carboxyglutamic_acid":   {"C":1, "O":2},
    "dehydralation":                {"H":-2, "O": -1},
    "oxidation":                    {"O":1},
    "carboxymethylation":           {"C":2, "H":2, "O":2},
    "palmitic_acylation":           {"C":16,"H":30,"O":1},
    "oleic_acylation":              {"C":18,"H":32,"O":1},
    "arachidonic_acylation":        {"C":20,"H":30,"O":1},
    "docosahexanoic_acylation":     {"C":22,"H":30,"O":1},
    "carbamidomethylation":         {"C":2, "H":3, "N":1, "O":1}
}

typical_ptm_sites = {
    "phosphorylation": ("S", "T", "Y")
}

plgs_mod_pattern = re.compile("(\w+)\+(\w)\((\d+)\)")
# this should match: Oxidation+M(24), Carbamidomethyl+C(7), ...

plgs_notation = {
    'Carbamidomethyl':'carbamidomethylation',
    'Oxidation':'oxidation'
}
# this should be prolonged

def compare_seq_and_mod(seq, mod):
    assert all(seq[int(i)-1] == aa for _, aa, i in plgs_mod_pattern.findall(mod))


def test_parsing():
    """This tests that Stefan is right about +C(17) meaning that the 17th amino acid is C."""
    for seq, mod in zip(sm.sequence, sm.modification):
        compare_seq_and_mod(seq, mod)


def PLGSptms2atom(mod, seq=''):
    """For a given PLGS modification, calculate the atom_diff.

    Args:
        mod (str): A string with PLGS modification (see 'plgs_notation').
        seq (str): The amino acidic sequence that must be checked.
    Returns:
        collections.Counter: An atom_diff.
    """
    global plgs_mod_pattern, plgs_notation
    atom_diff = Counter()
    for name, aa, i in plgs_mod_pattern.findall(mod):
        atom_diff.update(ptms[plgs_notation[name]])
        if seq:
            assert seq[int(i)-1] == aa, "The modification is poorly specified!"
    return atom_diff
