import mbuild as mb
from foyer import Forcefield
import parmed as pmd

water = mb.load('files/tip3p.mol2')
water.name = 'water'
acn = mb.load('files/acn.mol2')
acn.name = 'acn'

C = Forcefield('files/carbon.xml')
SPCE = Forcefield('files/spce.xml')
OPLS = Forcefield(name='oplsaa')

system = mb.recipes.GraphenePoreSolvent(pore_depth=4, side_dim=4, n_sheets=3,
                             pore_width=1.2, x_bulk=3, solvent=[water, acn],
                             n_solvent=[1000, 1000])
import pdb; pdb.set_trace()

box = mb.Box(system.periodicity)

water_pmd = pmd.Structure()
acn_pmd = pmd.Structure()
pore_pmd = pmd.Structure()
for child in system.children:
    if child.name == 'water':
        water_pmd += SPCE.apply(child, residues='water')
    elif child.name == 'acn':
        acn_pmd += OPLS.apply(child, residues='acn')
    else:
        pore_pmd += C.apply(child)
system = water_pmd + acn_pmd + pore_pmd
system.box = [0, 0, 0, 90, 90, 90]
system.box[:3] = box.maxs * 10.0
system.save('solvated_pore.gro', overwrite=True, combine='all')
system.save('solvated_pore.top', overwrite=True, combine='all')
