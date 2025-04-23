import MDAnalysis as mda
import os
import numpy as np
import shutil
import subprocess as sp
import sys
import datetime
from funnel_maker import *

pdb_file = './4erc-4np.eq6.pdb' # MDM2 protein with Nutlin 3A

p0, p1 = make_funnel(pdb_file, ligand_name = '4NP', output_pymol_session=True)
struct_file = './4erc-4np.eq6.rst7' 
top_file = '../4erc-4np.prmtop' # HSP90 protein-ligand system, amber parmtop file

p0, p1 = make_funnel(struct_file, top_file, ligand_name = '4NP', output_pymol_session=True, p0_residue=42)
#p0, p1 = make_funnel(struct_file, top_file, ligand_name = '4NP', output_pymol_session=True)

# get all of the ligand and protein atom IDs
p_ids, l_ids = get_protein_ligand_ids(struct_file, top_file, ligand_name = '4NP')

# write a plumed file in the current working directory
write_plumed_file(p0, p1, p_ids, l_ids)
