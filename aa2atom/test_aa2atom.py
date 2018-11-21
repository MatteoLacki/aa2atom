from aa2atom import aa2atom


def test_if_aa2atom_works_like_before():
    assert aa2atom("ACEF") == {'H': 11, 'C': 9, 'O': 2, 'N': 1}

def test_if_aa2atom_adds_water_properly():
    x = aa2atom("ACEF")
    y = aa2atom("ACEF", add_water=False)
    assert y['H'] + 2 == x['H']
    assert y['O'] + 1 == x['O']