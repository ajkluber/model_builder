""" Holding properties of residues


Currently not being used for anything
"""

def residue_mass(code):
    """Masses in atomic mass units."""
    residue_mass = {'ALA':   89.0935, 'ARG':  174.2017, 'ASN':  132.1184,
                    'ASP':  133.1032, 'CYS':  121.1590, 'GLN':  146.1451,
                    'GLU':  147.1299, 'GLY':   75.0669, 'HIS':  155.1552,
                    'ILE':  131.1736, 'LEU':  131.1736, 'LYS':  146.1882,
                    'MET':  149.2124, 'PHE':  165.1900, 'PRO':  115.1310,
                    'SER':  105.0930, 'THR':  119.1197, 'TRP':  204.2262,
                    'TYR':  181.1894, 'VAL':  117.1469, 'SOL':   18.0150}[code]
    return residue_mass

def residue_three_to_one_letter_code(code):
    """Converting from three letter code to one letter FASTA code."""
    residue_code = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N',
                    'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
                    'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
                    'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                    'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                    'SER': 'S', 'THR': 'T', 'TRP': 'W',
                    'TYR': 'Y', 'VAL': 'V'}[code]
    return residue_code
