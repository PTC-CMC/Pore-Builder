from porebuilder import GraphenePore
import mbuild as mb
import numpy as np
import foyer
from foyer import Forcefield
import parmed as pmd

C = Forcefield('files/carbon.xml')

system = GraphenePore(x_sheet=3, y_sheet=3, sheets=3, pore_width=1.2)

box = mb.Box(system.periodicity)

system_pmd = system.to_parmed(residues=['Compound'])

system_pmd.box = np.empty(6)
system_pmd.box[:3] = box.maxs * 10 # convert from nm to angstroms
system_pmd.box[3:7] = 90
system_pmd.save('init.gro', overwrite=True)
system_pmd.save('init.top', overwrite=True)
