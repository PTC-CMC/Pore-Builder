from porebuilder import gph_pore_solv
import mbuild as mb
import numpy as np
import foyer
from foyer import Forcefield
import parmed as pmd

water = '/Users/raymatsumoto/science/il_solvent_local/file_gen/mol2/tip3p.mol2'
ch3cn = '/Users/raymatsumoto/science/il_solvent_local/file_gen/mol2/ch3cn.mol2'
system = gph_pore_solv(x_sheet=4, y_sheet=4, sheets=3, pore_width=1.2,
        x_bulk=3, solvent=[{'SOL': water}, {'ch3cn':ch3cn}], n_solvent=[1000,300])
C_spce = 'C-spce.xml'
C_spce = Forcefield(C_spce)
opls = Forcefield(name='oplsaa')

box = mb.Box(system.box)
totalPM = pmd.Structure()
for child in system.children:
    if child.name in 'Compound':
        totalPM += child.to_parmed(residues='Compound', box=box)
    elif child.name in system.fluid_name[1]:
        totalPM += child.to_parmed(residues='ch3cn', box=box)
    elif child.name in system.fluid_name[0]:
        totalPM += child.to_parmed(residues='SOL', box =box)

ch3cnPM = totalPM['ch3cn',:]
SOLPM = totalPM['SOL',:]
gphPM = totalPM['Compound',:]
import pdb; pdb.set_trace()

SOLPM = C_spce.apply(SOLPM, residues='SOL')
gphPM = C_spce.apply(gphPM, residues='Compound')
ch3cnPM = opls.apply(ch3cnPM, residues='ch3cn')

systemPM = SOLPM + gphPM + ch3cnPM 
systemPM.box = np.empty(6)
systemPM.box[:3] = box.maxs * 10 # convert from nm to angstroms
systemPM.box[3:7] = 90
systemPM.save('init.gro', overwrite=True)
systemPM.save('init.top', overwrite=True)
