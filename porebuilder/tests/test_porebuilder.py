import numpy as np
import mbuild as mb
from porebuilder.tests.base_test import BaseTest


class TestPoreBuilder(BaseTest):
    """
    Unit Tests for Pore class functionality.
    """

    def test_save_dry(self, GraphenePore):
        GraphenePore.save(filename='dry_pore.gro')

    def test_save_solvated(self, GraphenePoreSolvent):
        GraphenePoreSolvent.save(filename='solvated_pore.gro')

    def test_save_double(self, DoubleGraphenePore):
        DoubleGraphenePore.save(filename='double_pore.gro')

    def test_hierarchy_dry(self, GraphenePore):
        assert len(GraphenePore.children) == 2

    def test_hierarchy_solvated(self, GraphenePoreSolvent):
        assert len(GraphenePoreSolvent.children) == 11
        lens = [2] + 10 * [3]
        assert [len(c.children) for c in GraphenePoreSolvent.children] == lens

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
            assert GraphenePoreSolvent.n_particles == 2046
        elif GraphenePore:
            assert GraphenePore.n_particles == 2016

    def test_particles_in_box(self, GraphenePoreSolvent):
        box = mb.Box(GraphenePoreSolvent.periodicity)

        for particle in GraphenePoreSolvent.particles():
            assert particle.xyz[0][0] < box.maxs[0]
            assert particle.xyz[0][1] < box.maxs[1]
            assert particle.xyz[0][2] < box.maxs[2]
