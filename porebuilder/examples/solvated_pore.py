from porebuilder import GraphenePoreSolvent
import mbuild as mb
from foyer import Forcefield

water = mb.load('files/tip3p.mol2')
water.name = 'water'
na = mb.load('files/na.mol2')
na.name = 'na'
cl = mb.load('files/cl.mol2')
cl.name = 'cl'

C = Forcefield('files/carbon.xml')
OPLS = Forcefield('files/oplsaa.xml')

system = GraphenePoreSolvent(pore_depth=4, side_dim=4, n_sheets=3,
                             pore_width=1.2, x_bulk=3, solvent=[na, cl, water],
                             n_solvent=[100, 100, 4000])

box = mb.Box(system.periodicity)

system = OPLS.apply(system)
system.save('solvated_pore.gro', overwrite=True)
system.save('solvated_pore.top', overwrite=True, combine='all')
