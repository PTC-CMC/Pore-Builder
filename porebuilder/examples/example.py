from porebuilder import GraphenePoreSolvent
import mbuild as mb
import numpy as np
import foyer
from foyer import Forcefield
import parmed as pmd

water = 'files/tip3p.mol2'
acn = 'files/acn.mol2'

system = GraphenePoreSolvent(x_sheet=4, y_sheet=4, sheets=3, pore_width=1.2,
        x_bulk=3, solvent=[{'SOL': water}, {'acn':acn}], n_solvent=[1000,300])

C_spce = Forcefield('files/C-spce.xml')
opls = Forcefield(name='oplsaa')

box = mb.Box(system.box)
totalPM = pmd.Structure()
for child in system.children:
    if child.name in 'Compound':
        totalPM += child.to_parmed(residues='Compound', box=box)
    elif child.name in system.fluid_name[1]:
        totalPM += child.to_parmed(residues='acn', box=box)
    elif child.name in system.fluid_name[0]:
        totalPM += child.to_parmed(residues='SOL', box =box)

acnPM = totalPM['acn',:]
SOLPM = totalPM['SOL',:]
gphPM = totalPM['Compound',:]
import pdb; pdb.set_trace()

SOLPM = C_spce.apply(SOLPM, residues='SOL')
gphPM = C_spce.apply(gphPM, residues='Compound')
acnPM = opls.apply(acnPM, residues='acn')

systemPM = SOLPM + gphPM + acnPM
systemPM.box = np.empty(6)
systemPM.box[:3] = box.maxs * 10 # convert from nm to angstroms
systemPM.box[3:7] = 90
systemPM.save('init.gro', overwrite=True)
systemPM.save('init.top', overwrite=True)
