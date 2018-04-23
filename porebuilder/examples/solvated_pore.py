from porebuilder import GraphenePoreSolvent
import mbuild as mb
from foyer import Forcefield

water = mb.load('files/tip3p.mol2')
water.name = 'water'
acn = mb.load('files/acn.mol2')
acn.name = 'acn'

C = Forcefield('files/carbon.xml')
SPCE = Forcefield('files/spce.xml')
OPLS = Forcefield(name='oplsaa')

system = GraphenePoreSolvent(x_sheet=4, y_sheet=4, sheets=3, pore_width=1.2,
        x_bulk=3, solvent=[water, acn], #[{'SOL': water}, {'acn':acn}],
        n_solvent=[1, 1])

box = mb.Box(system.periodicity)

for child in system.children:
    if child.name == 'water':
        water_pmd = SPCE.apply(child)
    elif child.name == 'acn':
        acn_pmd = OPLS.apply(child)
    else:
        pore_pmd = C.apply(child)

system = water_pmd + acn_pmd + pore_pmd
system.box[:3] = box.maxs * 10.0
system.save('solvated_pore.gro', overwrite=True)
system.save('solvated_pore.top', overwrite=True, combine='all')
