import numpy as np
import pytest
import mbuild as mb
import parmed as pmd
from porebuilder.tests.base_test import BaseTest

class TestPoreBuilder(BaseTest):
    """
    Unit Tests for Pore class functionality.
    """

    def test_save(self, GraphenePore_nosolv, GraphenePoreSolvent):
        GraphenePore_nosolv.save(filename='GraphenePore_nosolv.gro')
        GraphenePoreSolvent.save(filename='GraphenePore.gro')

    def test_porewidth(self, GraphenePoreSolvent):
        bot = np.min(GraphenePoreSolvent.bot_xyz[:,1])
        top = np.max(GraphenePoreSolvent.top_xyz[:,1]) # change naming
        pore_width = bot - top
        np.testing.assert_almost_equal(pore_width, 1, 4)

    def test_sheet_dims(self, GraphenePoreSolvent):
        x_length = np.max(GraphenePoreSolvent.bot_xyz[:,0])
        - np.min(GraphenePoreSolvent.bot_xyz[:,0])
        assert x_length == pytest.approx(3, 0.5)
        y_length = np.max(GraphenePoreSolvent.bot_xyz[:,2])
        - np.min(GraphenePoreSolvent.bot_xyz[:,2])
        assert y_length == pytest.approx(3, 0.5)
    
    def test_n_particles(self, GraphenePoreSolvent, GraphenePore_nosolv):
        if GraphenePoreSolvent:
            assert GraphenePoreSolvent.n_particles == 5016
        elif GraphenePore_nosolv:
            assert GraphenePore_nosolv.n_particles == 2016
    
    def test_particles_in_box(self, GraphenePoreSolvent):
        box = mb.Box(GraphenePoreSolvent.box)
        totalPM = pmd.Structure()
        for child in GraphenePoreSolvent.children:
            if child.name in 'Compound':
                totalPM += child.to_parmed(residues='Compound', box=box)
            elif child.name in GraphenePoreSolvent.fluid_name:
                totalPM += child.to_parmed(residues='SOL', box=box)
        gphPM = totalPM['Compound',:]
        SOLPM = totalPM['SOL', :]
        systemPM = gphPM + SOLPM
        systemPM.box = np.empty(6)
        systemPM.box[:3] = box.maxs * 10 
        systemPM.box[3:7] = 90
        for position in GraphenePoreSolvent.xyz:
            for x in range(3):
                assert position[x] < systemPM.box[x]
                assert position[x] >= 0.0 # May have to change this test

    # TODO: Get correct x-dimension of graphene sheet
    """def test_x_bulk_length(self, GraphenePoreSolvent):
        gph_length = np.max(GraphenePoreSolvent.bot_xyz[0])
        np.min(GraphenePoreSolvent.bot_xyz[0])
        assert (GraphenePoreSolvent.periodicity[0] - gph_length)/2 == # need value"""
