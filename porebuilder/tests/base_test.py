import pytest
import mbuild as mb
import os
from pkg_resources import resource_filename

TESTFILE_DIR = resource_filename('porebuilder', 'tests/test_molecules')


class BaseTest:

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def GraphenePore(self):
        from porebuilder.porebuilder import GraphenePore
        return GraphenePore(pore_depth=3, side_dim=3, n_sheets=3, pore_width=1)

    @pytest.fixture
    def GraphenePoreSolvent(self):
        from porebuilder.porebuilder import GraphenePoreSolvent
        h2o = mb.load(os.path.join(TESTFILE_DIR, 'tip3p.mol2'))
        h2o.name = 'SOL'
        return GraphenePoreSolvent(pore_depth=3, side_dim=3, n_sheets=3,
                                   pore_width=1, x_bulk=3, solvent=[h2o],
                                   n_solvent=10)
