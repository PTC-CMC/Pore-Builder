from porebuilder.porebuilder import Pores
import mbuild as mb
import foyer
from foyer import Forcefield

water = '/Users/raymatsumoto/science/il_solvent_local/file_gen/mol2/tip3p.mol2'
system = Pores(x_sheet=4, y_sheet=4, sheets=3, pore_width=1.2,
        x_bulk=3, solvent={'SOL': water}, n_solvent=2500)
ff = 'C-spce.xml'
ff = Forcefield(ff)
system = ff.apply(system.to_parmed(residues=['RES', 'SOL']))

system.save('init.gro', overwrite=True)
print('atomtyping...')
system.save('init.top', overwrite=True, combine='all')
