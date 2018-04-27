from porebuilder import GraphenePore
from foyer import Forcefield

C = Forcefield('files/carbon.xml')

system = GraphenePore(pore_depth=4, side_dim=2, n_sheets=3, pore_width=1.2)

system = C.apply(system)

system.save('dry_pore.gro', overwrite=True)
system.save('dry_pore.top', overwrite=True, combine='all')
