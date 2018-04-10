from porebuilder import GraphenePoreSolvent
import mbuild as mb
import numpy as np
import foyer
from foyer import Forcefield
import parmed as pmd

water = mb.load('files/tip3p.mol2')
acn = mb.load('files/acn.mol2')

C_spce = Forcefield('files/C-spce.xml')
opls = Forcefield(name='oplsaa')

system = GraphenePoreSolvent(x_sheet=4, y_sheet=4, sheets=3, pore_width=1.2,
        x_bulk=3, solvent=[{'SOL': water}, {'acn':acn}], n_solvent=[1000,300])


box = mb.Box(system.box)
system_pmd = pmd.Structure()

for child in system.children:
    if child.name in 'Compound':
        system_pmd += child.to_parmed(residues='Compound', box=box)
    elif child.name in system.fluid_name[1]:
        system_pmd += child.to_parmed(residues='acn', box=box)
    elif child.name in system.fluid_name[0]:
        system_pmd += child.to_parmed(residues='SOL', box =box)

acn_pmd = system_pmd['acn',:]
SOL_pmd = system_pmd['SOL',:]
gph_pmd = system_pmd['Compound',:]

SOL_pmd = C_spce.apply(SOL_pmd, residues='SOL')
gph_pmd = C_spce.apply(gph_pmd, residues='Compound')
acn_pmd = opls.apply(acn_pmd, residues='acn')

system_pmd = SOL_pmd + gph_pmd + acn_pmd
system_pmd.box = np.empty(6)
system_pmd.box[:3] = box.maxs * 10 # convert from nm to angstroms
system_pmd.box[3:7] = 90
system_pmd.save('init.gro', overwrite=True)
system_pmd.save('init.top', overwrite=True)
