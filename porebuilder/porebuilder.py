import itertools as it
import mbuild as mb
import numpy as np
from mbuild.compound import Compound
from mbuild.coordinate_transform import force_overlap
from mbuild.utils.validation import assert_port_exists
from mbuild import clone
from copy import deepcopy
import math

__all__ = ['GraphenePoreSolvent', 'GraphenePore']


class GraphenePoreSolvent(mb.Compound):
    """A general slit pore recipe that solvates the system

    Parameters
    ----------
    x_sheet : int
        dimensions of graphene sheet in x-direction [nm]
    y_sheet : int
        dimensions of graphene sheet in y-direction [nm]
    sheets : int
        number of graphene sheets, default=3
    pore_width: int
        width of slit pore [nm]
    x_bulk : int
        length of bulk region in x-direction [nm]
    solvent : dict or an array of 2 dicts
        compound(s) to solvate the system with.
        Note: does not currently support more than 2 solvents
    n_solvent : int
        number of solvents to solvate the system with. Array must
        match size of 'solvent'
    Attributes
    ----------

    """
    def __init__(self, x_sheet, y_sheet, sheets, pore_width, x_bulk,
            solvent, n_solvent):
        super(GraphenePoreSolvent, self).__init__()

        system = GraphenePore(x_sheet=x_sheet, y_sheet=y_sheet,
                              sheets=sheets, pore_width=pore_width)

        system.periodicity[2] += 2 * x_bulk
        box = mb.Box(system.periodicity)

        system = mb.solvate(system, solvent, n_solvent, box=box, overlap=0.2)

        for child in system.children:
            self.add(clone(child))

class GraphenePore(mb.Compound):
    """A general slit pore recipe.  Does not solvate system.  Use
    'GraphenePoreSolvent' instead if you wish to solvate your system.

    Parameters
    ----------
    x_sheet : int
        dimensions of graphene sheet in x-direction [nm]
    y_sheet : int
        dimensions of graphene sheet in y-direction [nm]
    sheets : int
        number of graphene sheets, default=3
    pore_width : int
        width of slit pore [nm]
    x_bulk : int
        length of bulk region in x-direction [nm]
    Attributes
    ----------

    Notes: Match graphene y-dimension with box x-dimension
    """
    def __init__(self, x_sheet, y_sheet, sheets, pore_width):
        super(GraphenePore, self).__init__()

        # Do some math to figure out how much to replicate graphene cell.
        # Multiply replicate[1] by 15/13 to take into account later
        # multiplication
        # TODO: Figure out if rounding is necessary
        factor = np.cos(math.pi/6)
        replicate = [(x_sheet/0.2456),
                (y_sheet/0.2456)*(1/factor)]
        if all(x <= 0 for x in [x_sheet, y_sheet]):
            msg = 'Dimension of graphene sheet must be greater than zero'
            raise ValueError(msg)
        carbon = mb.Compound()
        carbon.name = 'C'
        carbon_locations = [[0, 0, 0], [2/3, 1/3, 0]]
        basis = {carbon.name: carbon_locations}
        graphene_lattice = mb.Lattice(lattice_spacing=[0.2456, 0.2456, 0.335],
                                      angles=[90, 90, 120], lattice_points=basis)
        graphene = graphene_lattice.populate(compound_dict={carbon.name: carbon},
                                             x=replicate[0], y=replicate[1],
                                             z=sheets)
        for particle in graphene.particles():
            if particle.xyz[0][0] < 0:
                particle.xyz[0][0] += graphene.periodicity[0]
        graphene_dims = graphene.periodicity
        graphene_dims[1] *= factor  # cos(30)*.246
        bottom_sheet = mb.clone(graphene)
        bottom_sheet.translate([0, pore_width + (graphene_dims[2] - 0.335), 0])
        bottom_sheet.spin(1.5708, [1, 0, 0])
        bottom_sheet.name = 'BOT'
        top_sheet = mb.clone(graphene)
        top_sheet.spin(1.5708, [1, 0, 0])
        top_sheet.name = 'TOP'
        self.add(top_sheet)
        self.add(bottom_sheet)
