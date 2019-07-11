from collections import Counter

amino_acids = {  'A': {'H': 5, 'C': 3, 'O': 1, 'N': 1},
                 'C': {'H': 5, 'C': 3, 'S': 1, 'O': 1, 'N': 1},
                 'E': {'H': 7, 'C': 5, 'O': 3, 'N': 1},
                 'D': {'H': 5, 'C': 4, 'O': 3, 'N': 1},
                 'G': {'H': 3, 'C': 2, 'O': 1, 'N': 1},
                 'F': {'H': 9, 'C': 9, 'O': 1, 'N': 1},
                 'I': {'H': 11, 'C': 6, 'O': 1, 'N': 1},
                 'J': {'H': 11, 'C': 6, 'O': 1, 'N': 1},
                 'H': {'H': 7, 'C': 6, 'N': 3, 'O': 1},
                 'K': {'H': 12, 'C': 6, 'N': 2, 'O': 1},
                 'M': {'H': 9, 'C': 5, 'S': 1, 'O': 1, 'N': 1},
                 'L': {'H': 11, 'C': 6, 'O': 1, 'N': 1},
                 'N': {'H': 6, 'C': 4, 'O': 2, 'N': 2},
                 'O': {'H': 19,'C':12, 'O': 2, 'N': 3},
                 'Q': {'H': 8, 'C': 5, 'O': 2, 'N': 2},
                 'P': {'C': 5, 'O': 1, 'H': 7, 'N': 1},
                 'S': {'H': 5, 'C': 3, 'O': 2, 'N': 1},
                 'R': {'H': 12, 'C': 6, 'N': 4, 'O': 1},
                 'T': {'H': 7, 'C': 4, 'O': 2, 'N': 1},
                 'U': {'H': 5, 'C': 3, 'Se': 1, 'O': 1, 'N': 1},
                 'W': {'H': 10, 'C': 11, 'N': 2, 'O': 1},
                 'V': {'H': 9, 'C': 5, 'O': 1, 'N': 1},
                 'Y': {'H': 9, 'C': 9, 'O': 2, 'N': 1}}

aa2info = {  'A': ('Ala', 'Alanine', 'GCT', 'GCG', 'GCA', 'GCG'),
             'C': ('Cys', 'Cysteine', 'TGT', 'TGC'),
             'E': ('Glu', 'Glutamic Acid', 'GAA', 'GAG'),
             'D': ('Asp', 'Aspartic Acid', 'GAT', 'GAC'),
             'G': ('Gly', 'Glycine', 'GGT', 'GGC', 'GGA', 'GGG'),
             'F': ('Phe', 'Phenylalanine', 'TTT', 'TTC'),
             'I': ('Ile', 'Isoleucine', 'ATT', 'ATC', 'ATA'),
             'J': ('Ile/Leu', 'Isoleucine or Leucine'),
             'H': ('His', 'Histidine', 'CAT, CAC'),
             'K': ('Lys', 'Lysine', 'AAA', 'AAG'),
             'M': ('Met', 'Methionine', 'ATG'),
             'L': ('Leu', 'Leucine', 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'),
             'N': ('Asn', 'Asparagine', 'AAT', 'AAC'),
             'O': ('Pyl', 'Pyrrolysine'),
             'Q': ('Gln', 'Glutamine', 'CAA', 'CAG'),
             'P': ('Pro', 'Proline', 'CCT', 'CCC', 'CCA', 'CCG'),
             'S': ('Ser', 'Serine', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
             'R': ('Arg', 'Arginine', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
             'T': ('Thr', 'Threonine', 'ACT', 'ACC', 'ACA', 'ACG'),
             'U': ('Sec', 'Selenocysteine', 'TGT', 'TGC'),
             'W': ('Trp', 'Tryptophan', 'TGG'),
             'V': ('Val', 'Valine', 'GTT', 'GTC', 'GTA', 'GTG'),
             'Y': ('Tyr', 'Tyrosine', 'TAT', 'TAC')}



aa2shortName = {aa: v[0] for aa, v in aa2info.items()}
aa2name = {aa: v[1] for aa, v in aa2info.items()}
aa2encodingSequence = {aa: v[2:] for aa, v in aa2info.items()}


class UnknownAminoAcid(Exception):
    """Exception raised for errors in the input.

    Attributes:
        res: a partial result before the exception was raised. It is stored at property "res" of the exception and can be accessed later on.
        unknown_aas (dict): unspecified amino acids and their counts.  
    """

    def __init__(self, res, unknown_aas, *args, **kwrg):
        super().__init__(*args, **kwrg)
        self.res = res
        self.unknown_aas = unknown_aas


def aa2atom(aaseq, no_water=False):
    """Get the atom counts for an amino acidic sequence.

    Arguments
    =========
    aaseq : str
        A sequence of amino acids, eg. "AAQPQQAA".
    Returns
    =======
    collections.Counter : counts of atoms.
    """
    aacnts = Counter(aaseq)
    res = Counter()
    error = False
    unknown_aas = {}
    for a, aacnt in aacnts.items():
        try:
            for atom, cnt in amino_acids[a].items():
                res[atom] += aacnt*cnt
        except KeyError:
            unknown_aas[a] = unknown_aas.get(a, 0) + aacnt
            error = True
    if not no_water:
        res['H'] += 2
        res['O'] += 1
    if error:
        raise UnknownAminoAcid(res, unknown_aas, "Sequence contained unspecified entries: {}.".format(str(unknown_aas)))
    else:
        return res


def atom2str(atoms):
    """Present atoms in the string form.

    Arguments
    =========
    atoms : dict-like
        Keys correspond to elements and values to atom counts.

    Returns
    =======
    str : The chemical formula (compatible with IsoSpec).
    """
    return "".join(k + str(atoms[k]) for k in sorted(list(atoms.keys())))


def atom2tex(atoms, ce=False):
    f = "".join("{element}_{{{count}}}".format(element=element, count=count) if count > 1 else element
                for element, count in sorted(atoms.items()))
    if ce:
        return "\ce{" + f + "}"
    else:
        return "$"+f+"$"

