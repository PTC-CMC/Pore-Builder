import numpy as np
import pytest
import mbuild as mb
from mbuild.utils.geometry import calc_dihedral
import os
from pkg_resources import resource_filename

TESTFILE_DIR = resource_filename('porebuilder', 'tests/test_molecules')

class BaseTest:

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def gph_pore_nosolv(self):
        from porebuilder.porebuilder import gph_pore
        return gph_pore(x_sheet=3, y_sheet=3, sheets=3, pore_width=1,
                x_bulk=3)

    @pytest.fixture
    def gph_pore_solv(self):
        from porebuilder.porebuilder import gph_pore_solv
        h2o = os.path.join(TESTFILE_DIR, 'tip3p.mol2')
        return gph_pore_solv(x_sheet=3, y_sheet=3, sheets=3, pore_width=1,
                x_bulk=3, solvent={'SOL': h2o}, n_solvent=1000)
