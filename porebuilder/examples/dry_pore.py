from porebuilder import GraphenePore
from foyer import Forcefield

C = Forcefield('files/carbon.xml')

system = GraphenePore(x_sheet=4, y_sheet=2, sheets=3, pore_width=1.2)

system = C.apply(system)

system.save('dry_pore.gro', overwrite=True)
system.save('dry_pore.top', overwrite=True, combine='all')
