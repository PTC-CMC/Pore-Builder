import numpy as np
import pytest

import mbuild as mb
from mbuild.utils.geometry import calc_dihedral

class BaseTest:

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def gph_pore_nosolv(self):
        from porebuilder.porebuilder import Pores
        return Pores(x_sheet=3, y_sheet=3, sheets=3, pore_width=1,
                x_bulk=3)

    def gph_pore_solv(self):
        from porebuilder.porebuilder import Pores
        h2o = 'testmolecules/tip3p.mol2'
        return Pores(x_sheet=3, y_sheet=3, sheets=3, pore_width=1,
                x_bulk=3, solvent={'SOL': H2O}, n_solvent=1000)
