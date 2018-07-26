import numpy as np
import mbuild as mb
import pytest

from porebuilder.tests.base_test import BaseTest
import porebuilder as pb


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

    def test_functionalize_surface(self):

        class H(mb.Compound):
            def __init__(self):
                super(H, self).__init__()

                self.add(mb.Particle(name='H'))
                up_port = mb.Port(
                    anchor=self[0], orientation=[
                        0, 1, 0], separation=.075)
                self.add(up_port, "up")

        class O(mb.Compound):
            def __init__(self):
                super(O, self).__init__()

                self.add(mb.Particle(name='O'))
                up_port = mb.Port(
                    anchor=self[0], orientation=[
                        0, 1, 0], separation=.075)
                self.add(up_port, "up")

        with pytest.raises(ValueError):
            pb.GraphenePoreFunctionalized(
                func_groups=[H(), O()], func_percent=[.03, .04, .03])

        for per in range(0, 3):
            per = per * 3 + 1 
            pore = pb.GraphenePoreFunctionalized(
                func_groups=H(), func_percent= per / 10)
            assert(pore.n_particles - 2688 >= (per - .15) / 10 * 864) 
            assert(pore.n_particles - 2688 <= (per + .15) / 10 * 864)

        odd_pore = pb.GraphenePoreFunctionalized(n_sheets=4, pore_width=1.5,pore_depth=2, side_dim=2, func_groups=H(), func_percent=.5)
        
        hydrogens = odd_pore.particles_by_name('H')
        positions = [[],[],[]]
        for h in hydrogens:
            for array in h.xyz:
                for d,i in zip(array,range(0,3)):
                    positions[i].append(d)
        
        assert(np.amin(positions[1]) >= .335 * (3))
        assert(np.amax(positions[1]) <= .336 * (3) + 1.5)
        assert(np.amin(positions[0]) >= -(0.0001))
        assert(np.amax(positions[0]) <= 2)
        assert(np.amin(positions[2]) >= -(0.0001))
        assert(np.amax(positions[2]) <= 2)

        assert(odd_pore.n_bonds == 140)