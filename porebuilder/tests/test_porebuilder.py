import numpy as np
import pytest
import mbuild as mb
import base_test import BaseTest

class TestPoreBuilder(BaseTest):
    """
    Unit Tests for Pore class functionality.
    """

    def test_save(self, gph_pore_solv):
        gph_pore_solv.save(filename='gph_pore.gro')

    def test_porewidth(self, gph_pore_solv):
        pore = gph_pore_solv() # may need to correct
        bot = np.min(pore.bottom_sheet.xyz[:,1])
        top = np.max(pore.top_sheet.xyz[:,1]) # change naming
        pore_width = bot - top
        np.testing.assert_almost_equal(pore_width, 1, 4) 
        # may have to change rounding precision
