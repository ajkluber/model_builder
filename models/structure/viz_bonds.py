
def write_bonds_tcl(topology, outfile="bonds.tcl"):
    molid = 0
    bondstring = lambda molid, idx1, idx2: \
'''set sel [atomselect {0} "index {1} {2}"]
lassign [$sel getbonds] bond1 bond2
set id [lsearch -exact $bond1 {2}]
if {{ $id == -1 }} {{
lappend bond1 {2}
}}
set id [lsearch -exact $bond2 {1}]
if {{ $id == -1 }} {{
lappend bond2 {1}
}}
$sel setbonds [list $bond1 $bond2]
$sel delete'''.format(molid, idx1, idx2)

    tclstring = ''
    for atm1, atm2 in topology.bonds:  
        tclstring += bondstring(molid,atm1.index, atm2.index) + "\n"

    with open(outfile, 'w') as fout:
        fout.write(tclstring)

def write_bonds_conect(topology, outfile="conect.pdb"):
    conectstring = ''
    for atm1, atm2 in topology.bonds:  
        conectstring += "CONECT{:5d}{:5d}\n".format(atm1.index + 1, atm2.index + 1)

    with open(outfile, 'w') as fout:
        fout.write(conectstring)
