import numpy as np
import pytest
import mbuild as mb
import parmed as pmd
from porebuilder.tests.base_test import BaseTest

class TestPoreBuilder(BaseTest):
    """
    Unit Tests for Pore class functionality.
    """

    def test_save_dry(self, GraphenePore):
        GraphenePore.save(filename='dry_pore.gro')

    def test_save_solvated(self, GraphenePoreSolvent):
        GraphenePoreSolvent.save(filename='solvated_pore.gro')

    def test_hierarchy_dry(self, GraphenePore):
        assert len(GraphenePore.children) == 2

    def test_hierarchy_solvated(self, GraphenePoreSolvent):
        assert len(GraphenePoreSolvent.children) == 11
        assert [len(c.children) for c in GraphenePoreSolvent.children] == [2] + 10 * [3]

    def test_porewidth(self, GraphenePore):
        bot = next(c for c in GraphenePore.children if c.name == 'BOT')
        top = next(c for c in GraphenePore.children if c.name == 'TOP')
        bot_y = np.max(bot.xyz[:, 1])
        top_y = np.min(top.xyz[:, 1])
        assert np.isclose(top_y - bot_y, 1.0, 3)

    def test_sheet_dims(self, GraphenePore):
        bot = next(c for c in GraphenePore.children if c.name == 'BOT')
        top = next(c for c in GraphenePore.children if c.name == 'TOP')
        x_length = np.ptp(bot.xyz[:, 0])
        y_length = np.ptp(bot.xyz[:, 1])
        assert np.isclose(x_length, 3, 1)
        assert np.isclose(y_length, 3, 1)

    def test_n_particles(self, GraphenePoreSolvent, GraphenePore):
        if GraphenePoreSolvent:
            assert GraphenePoreSolvent.n_particles == 5016
        elif GraphenePore:
            assert GraphenePore.n_particles == 2016
    
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
