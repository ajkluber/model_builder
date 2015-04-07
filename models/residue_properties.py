# coding=utf-8
""" Holding properties of residues"""

"""Masses in atomic mass units"""
residue_mass = {'ALA':   89.0935, 'ARG':  174.2017, 'ASN':  132.1184,
                'ASP':  133.1032, 'CYS':  121.1590, 'GLN':  146.1451,
                'GLU':  147.1299, 'GLY':   75.0669, 'HIS':  155.1552,
                'ILE':  131.1736, 'LEU':  131.1736, 'LYS':  146.1882,
                'MET':  149.2124, 'PHE':  165.1900, 'PRO':  115.1310,
                'SER':  105.0930, 'THR':  119.1197, 'TRP':  204.2262,
                'TYR':  181.1894, 'VAL':  117.1469, 'SOL':   18.0150}

"""Converting from three letter code to one letter FASTA code."""
residue_code = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N',
                'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
                'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
                'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W',
                'TYR': 'Y', 'VAL': 'V'}

"""Source: Partial molar volumes of proteins: amino acid side-chain
contributions derived from the partial molar volumes of some tripeptide
Atomic radii in nanometers"""
residue_radii = {'ALA': 0.1844827, 'ARG': 0.3134491, 'ASN': 0.2477519,
                'ASP': 0.2334602, 'CYS': 0.2276212, 'GLN': 0.2733978,
                'GLU': 0.2639170, 'GLY': 0.0000000, 'HIS': 0.2835556,
                'ILE': 0.2889931, 'LEU': 0.2887070, 'LYS': 0.2937731,
                'MET': 0.2916368, 'PHE': 0.3140150, 'PRO': 0.2419109,
                'SER': 0.1936102, 'THR': 0.2376198, 'TRP': 0.3422321,
                'TYR': 0.3168939, 'VAL': 0.2619603, 'AVERAGE': 0.2683678}
