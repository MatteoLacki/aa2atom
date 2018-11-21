from linearCounter import LC

amino_acids = {  'A': {'H': 5, 'C': 3, 'O': 1, 'N': 1},
                 'C': {'H': 5, 'C': 3, 'S': 1, 'O': 1, 'N': 1},
                 'E': {'H': 7, 'C': 5, 'O': 3, 'N': 1},
                 'D': {'H': 5, 'C': 4, 'O': 3, 'N': 1},
                 'G': {'H': 3, 'C': 2, 'O': 1, 'N': 1},
                 'F': {'H': 9, 'C': 9, 'O': 1, 'N': 1},
                 'I': {'H': 11, 'C': 6, 'O': 1, 'N': 1},
                 'H': {'H': 7, 'C': 6, 'N': 3, 'O': 1},
                 'K': {'H': 12, 'C': 6, 'N': 2, 'O': 1},
                 'M': {'H': 9, 'C': 5, 'S': 1, 'O': 1, 'N': 1},
                 'L': {'H': 11, 'C': 6, 'O': 1, 'N': 1},
                 'N': {'H': 6, 'C': 4, 'O': 2, 'N': 2},
                 'Q': {'H': 8, 'C': 5, 'O': 2, 'N': 2},
                 'P': {'C': 5, 'O': 1, 'H': 7, 'N': 1},
                 'S': {'H': 5, 'C': 3, 'O': 2, 'N': 1},
                 'R': {'H': 12, 'C': 6, 'N': 4, 'O': 1},
                 'T': {'H': 7, 'C': 4, 'O': 2, 'N': 1},
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
             'H': ('His', 'Histidine', 'CAT, CAC'),
             'K': ('Lys', 'Lysine', 'AAA', 'AAG'),
             'M': ('Met', 'Methionine', 'ATG'),
             'L': ('Leu', 'Leucine', 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'),
             'N': ('Asn', 'Asparagine', 'AAT', 'AAC'),
             'Q': ('Gln', 'Glutamine', 'CAA', 'CAG'),
             'P': ('Pro', 'Proline', 'CCT', 'CCC', 'CCA', 'CCG'),
             'S': ('Ser', 'Serine', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
             'R': ('Arg', 'Arginine', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
             'T': ('Thr', 'Threonine', 'ACT', 'ACC', 'ACA', 'ACG'),
             'W': ('Trp', 'Tryptophan', 'TGG'),
             'V': ('Val', 'Valine', 'GTT', 'GTC', 'GTA', 'GTG'),
             'Y': ('Tyr', 'Tyrosine', 'TAT', 'TAC')}

def aa2atom(aaseq):
    """Get the atom counts for an amino acidic sequence.

    Arguments
    =========
    aaseq : str
        A sequence of amino acids, eg. "AAQPQQAA".
    Returns
    =======
    dict : counts of atoms.
    """
    res = Counter()
    for aa, count in Counter(aaseq).items():
        for atom, acount in amino_acids[aa].items():
            res[atom] += count * acount
    return dict(res)

aa2atom
