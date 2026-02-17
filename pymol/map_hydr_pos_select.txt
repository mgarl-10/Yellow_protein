load C_alternans.pdb
hide everything
show cartoon
color gray
select hydro_res, resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO
show sticks, hydro_res
color yellow, hydro_res
select pos_sites, resi 156+189+226+271+276+281+415+700
show spheres, pos_sites
color red, pos_sites

