# finds smallest radii that can be used with gaussian contacts

import model_builder as mb
import numpy as np
import csv
from model_builder.models.residue_properties import *

#option = "average"
option = "flavored"
cushion = 0.01 # cushion to separate distances by (angstroms)
pairs = np.loadtxt('1E0Gcacb_conts')
pdb = 'clean_1E0G.pdb'

pairwise_distances = mb.models.pdb_parser.get_pairwise_distances(pdb,pairs)
min_radii = 10 * np.min(pairwise_distances) # convert from nm to angstroms
    
r_flavored = residue_radii.values()
r_keys = residue_radii.keys()

if option == "average":
    average_radii = residue_radii["AVERAGE"]
    if min_radii < average_radii:
        new_radii = min_radii-cushion
    else:
        new_radii = average_radii
    sum_radii = new_radii * (len(r_flavored)-1) # AVERAGE isn't a bead 
    with open('newradii.dat', 'wb') as f:
        writer = csv.writer(f)
        writer.writerow(("AVERAGE",new_radii))
        writer.writerow(("SUM",sum_radii))
elif option == "flavored":
    r_flavored_max = np.max(r_flavored)
    scale = (min_radii-cushion)/r_flavored_max
    r_scaled = [x*scale for x in r_flavored]
    sum_radii = np.sum(r_scaled)
    avg_radii = np.mean(r_scaled)
    for x in r_keys:
        if r_keys[x] == "AVERAGE":
            r_scaled[x] = avg_radii
    with open('newradii.dat', 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(zip(r_keys,r_scaled))
        writer.writerow(("SUM",sum_radii))

