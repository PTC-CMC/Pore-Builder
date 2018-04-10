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
    def GraphenePore_nosolv(self):
        from porebuilder.porebuilder import GraphenePore
        return GraphenePore(x_sheet=3, y_sheet=3, sheets=3, pore_width=1,
                x_bulk=3)

    @pytest.fixture
    def GraphenePoreSolvent(self):
        from porebuilder.porebuilder import GraphenePoreSolvent
        h2o = mb.load(os.path.join(TESTFILE_DIR, 'tip3p.mol2'))
        return GraphenePoreSolvent(x_sheet=3, y_sheet=3, sheets=3, pore_width=1,
                x_bulk=3, solvent={'SOL': h2o}, n_solvent=1000)
