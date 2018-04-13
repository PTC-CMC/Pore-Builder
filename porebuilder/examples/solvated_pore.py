from porebuilder import GraphenePoreSolvent
import mbuild as mb
import numpy as np
import foyer
from foyer import Forcefield
import parmed as pmd

water = mb.load('files/tip3p.mol2')
water.name = 'SOL'
acn = mb.load('files/acn.mol2')
acn.name = 'acn'

C = Forcefield('files/carbon.xml')
SPCE = Forcefield('files/spce.xml')
OPLS = Forcefield(name='oplsaa')

system = GraphenePoreSolvent(x_sheet=4, y_sheet=4, sheets=3, pore_width=1.2,
        x_bulk=3, solvent=[water, acn], #[{'SOL': water}, {'acn':acn}],
        n_solvent=[1, 1])

box = mb.Box(system.periodicity)
system_pmd = pmd.Structure()

for child in system.children:
    if child.name in 'GraphenePore':
        system_pmd += child.to_parmed(residues='GraphenePore', box=box)
    elif child.name in 'acn':
        system_pmd += child.to_parmed(residues='acn', box=box)
    elif child.name in 'SOL':
        system_pmd += child.to_parmed(residues='SOL', box =box)

acn_pmd = system_pmd['acn',:]
SOL_pmd = system_pmd['SOL',:]
gph_pmd = system_pmd['GraphenePore',:]

SOL_pmd = SPCE.apply(SOL_pmd)#, residues='SOL')
gph_pmd = C.apply(gph_pmd)#, residues='GraphenePore')
acn_pmd = OPLS.apply(acn_pmd)#, residues='acn')

system_pmd = SOL_pmd + gph_pmd + acn_pmd
system_pmd.box = np.empty(6)
system_pmd.box[:3] = box.maxs * 10 # convert from nm to angstroms
system_pmd.box[3:7] = 90
system_pmd.save('init.gro', overwrite=True)
system_pmd.save('init.top', overwrite=True)
