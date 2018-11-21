Simply, map a sequence of amino acids to the corresponding counts of atoms.

Installation
------------

```batch
pip3 install aa2atom
```

Usage
-----
```python
from aa2atom import aa2atom, atom2str, atom2tex

atoms = aa2atom("MDYSNFGNSASKKFQDDTLNRVRKEHEEA")
print(atoms)
# {'S': 1, 'N': 44, 'O': 51, 'C': 144, 'H': 222}

atoms_no_water = aa2atom("MDYSNFGNSASKKFQDDTLNRVRKEHEEA", add_water=False)
print(atoms_no_water)
# {'S': 1, 'N': 44, 'O': 50, 'C': 144, 'H': 220}

print(atom2str(atoms))
# C144H222N44O51S1

print(atom2tex(atoms))
# $C_{144}H_{222}N_{44}O_{51}S$

print(atom2tex(atoms, ce=True))
# \ce{C_{144}H_{222}N_{44}O_{51}S}
```

Terminal usage
```batch
aa2atom MDYSNFGNSASKKFQDDTLNRVRKEHEEA
# C144H222N44O51S1

aa2atom MDYSNFGNSASKKFQDDTLNRVRKEHEEA --nowater
# C144H220N44O50S1

aa2atom MDYSNFGNSASKKFQDDTLNRVRKEHEEA --TeX
# $C_{144}H_{222}N_{44}O_{51}S$

aa2atom MDYSNFGNSASKKFQDDTLNRVRKEHEEA --TeX --ce
# \ce{C_{144}H_{222}N_{44}O_{51}S}
```