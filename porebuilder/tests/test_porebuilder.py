import sys

import numpy as np
import mbuild as mb
import pytest
from porebuilder.tests.base_test import BaseTest


class TestPoreBuilder(BaseTest):
    """
    Unit Tests for Pore class functionality.
    """

    def test_imports(self):
        assert "porebuilder" in sys.modules
        assert "GraphenePore" in vars(mb.recipes).keys()
        assert "GraphenePoreSolvent" in vars(mb.recipes).keys()

    def test_save_dry(self, GraphenePore):
        GraphenePore.save(filename="dry_pore.gro", combine="all")

    def test_save_solvated(self, GraphenePoreSolvent):
        GraphenePoreSolvent.save(filename="solvated_pore.gro", combine="all")

    def test_save_surface(self, GrapheneSurface):
        GrapheneSurface.save(filename="graphene_surface.gro", combine="all")

    def test_hierarchy_dry(self, GraphenePore):
        assert len(GraphenePore.children) == 2

    def test_hierarchy_solvated(self, GraphenePoreSolvent):
        assert len(GraphenePoreSolvent.children) == 11
        lens = [2016] + 10 * [3]
        assert [c.n_particles for c in GraphenePoreSolvent.children] == lens

    def test_porewidth(self, GraphenePore):
        bot = next(c for c in GraphenePore.children if c.name == "BOT")
        top = next(c for c in GraphenePore.children if c.name == "TOP")
        bot_y = np.max(bot.xyz[:, 1])
        top_y = np.min(top.xyz[:, 1])
        assert np.isclose(top_y - bot_y, 1.0, 3)

    def test_surfacewidth(self, GrapheneSurface):
        graphene = next(c for c in GrapheneSurface.children)
        x_length = np.ptp(graphene.xyz[:, 0])
        y_length = np.ptp(graphene.xyz[:, 1])

        assert np.isclose(x_length, 3, 1)
        assert np.isclose(y_length, 3, 1)

    def test_surfacevacuum(self, GrapheneSurface):
        graphene = next(c for c in GrapheneSurface.children)
        z_length = np.max(graphene.xyz[:, 2])

        assert GrapheneSurface.periodicity[2] == z_length + 5.0

    def test_sheet_dims(self, GraphenePore):
        bot = next(c for c in GraphenePore.children if c.name == "BOT")
        top = next(c for c in GraphenePore.children if c.name == "TOP")
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

    def test_dimension_error(self):
        from porebuilder.porebuilder import GraphenePore, GrapheneSurface

        with pytest.raises(ValueError):
            GraphenePore(pore_length=0, pore_depth=0, n_sheets=3, pore_width=1)
            GrapheneSurface(x_length=0, y_length=0)

    @pytest.mark.parametrize("dim", (0, 1, 2))
    def test_slitpore_dims(self, dim):
        from porebuilder.porebuilder import GraphenePore, GrapheneSurface

        GraphenePore(slit_pore_dim=dim)

    def test_pore_in_center(self, GraphenePore, GraphenePoreSolvent):
        system = GraphenePore
        assert np.allclose(
            (system.periodicity - np.max(system.xyz, axis=0)),
            (np.min(system.xyz, axis=0) - np.array([0, 0, 0])),
        )
        system = GraphenePoreSolvent
        carbon = system.children[0]
        assert np.isclose(
            (carbon.periodicity[1] - np.max(carbon.xyz[:, 1])),
            (np.min(carbon.xyz[:, 1]) - 0),
        )
