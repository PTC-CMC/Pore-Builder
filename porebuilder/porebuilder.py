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
        self.x_sheet = x_sheet
        self.y_sheet = y_sheet
        self.sheets = sheets
        self.pore_width = pore_width
        self.x_bulk = x_bulk

        # Do some math to figure out how much to replicate graphene cell.
        # Multiply replicate[1] by 15/13 to take into account later
        # multiplication
        # TODO: Figure out if rounding is necessary
        factor = np.cos(math.pi/6)
        replicate = [(self.x_sheet/0.2456),
                (self.y_sheet/0.2456)*(1/factor)]
        if all(x <= 0 for x in [x_sheet, y_sheet]):
            msg = 'Dimension of graphene sheet must be greater than zero'
            raise ValueError(msg)
        self.name = 'C'
        carbon_locations = [[0, 0, 0], [2/3, 1/3, 0]]
        basis = {self.name: carbon_locations}
        graphene_lattice = mb.Lattice(lattice_spacing=[0.2456, 0.2456, 0.335],
                                      angles=[90, 90, 120], lattice_points=basis)
        carbon = mb.Compound(name=self.name)
        graphene = graphene_lattice.populate(compound_dict={self.name: carbon},
                                             x=replicate[0], y=replicate[1],
                                             z=self.sheets)
        for particle in graphene.particles():
            if particle.xyz[0][0] < 0:
                particle.xyz[0][0] += graphene.periodicity[0]
        self.graphene_dims = graphene.periodicity
        self.graphene_dims[1] *= factor  # cos(30)*.246
        bottom_sheet = mb.clone(graphene)
        bottom_sheet.translate([0, self.pore_width + (self.graphene_dims[2] - 0.335), 0])
        bottom_sheet.spin(1.5708, [1, 0, 0])
        top_sheet = mb.clone(graphene)
        top_sheet.spin(1.5708, [1, 0, 0])
        self.bot_xyz = bottom_sheet.xyz
        self.top_xyz = top_sheet.xyz
        system = mb.Compound()
        system.from_parmed(structure=bottom_sheet.to_parmed() + top_sheet.to_parmed())
        if len(solvent) == 1:
            for key, value in solvent.items():
                fluid = value
                fluid.name = key
        elif len(solvent) in [2, 3]:
            for key, value in solvent[0].items():
                fluid_1 = value
                fluid_1.name = key
            for key, value in solvent[1].items():
                fluid_2 = value
                fluid_2.name = key
            fluid = [fluid_1, fluid_2]
            if len(solvent) == 3:
                for key, value in solvent[2].items():
                    fluid_3 = value
                    fluid_3.name = key
                    fluid.append(fluid_3)
        elif len(solvent) > 3:
            msg = '"GraphenePoreSolvent" class currently only supports a maximum of 2 solvents'
            raise ValueError(msg)

        box = [(self.x_bulk*2)+self.graphene_dims[0]+.5,
                self.pore_width+(2 * (self.graphene_dims[2] - 0.335) + 0.5),
                self.graphene_dims[1]]
        system = mb.solvate(system, fluid, n_solvent, box=box, overlap=0.2)
        system.periodicity = box
        self.box = system.periodicity

        if len(solvent) == 1:
            for child in system.children:
                if child.name in fluid.name:
                    self.add(mb.clone(child))
                    self.fluid_name = fluid.name
                elif child.name in 'Compound':
                    self.add(mb.clone(child))
        elif len(solvent) in [2, 3]:
            self.fluid_name = [0 for x in range(2)]
            for child in system.children:
                if child.name in fluid[0].name:
                    self.add(mb.clone(child))
                    self.fluid_name[0] = fluid[0].name
                elif child.name in fluid[1].name:
                    self.add(mb.clone(child))
                    self.fluid_name[1] = fluid[1].name
                elif child.name in 'Compound':
                    self.add(mb.clone(child))
                if len(solvent) == 3:
                    if child.name in fluid[2].name:
                        self.add(mb.clone(child))
                        self.fluid_name.append(fluid[2].name)


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
    def __init__(self, x_sheet, y_sheet, sheets, pore_width, x_bulk):
        super(GraphenePore, self).__init__()
        self.x_sheet = x_sheet
        self.y_sheet = y_sheet
        self.sheets = sheets
        self.pore_width = pore_width
        self.x_bulk = x_bulk

        # Do some math to figure out how much to replicate graphene cell.
        # Multiply replicate[1] by 15/13 to take into account later
        # multiplication
        # TODO: Figure out if rounding is necessary
        factor = np.cos(math.pi/6)
        replicate = [(self.x_sheet/0.2456),
                (self.y_sheet/0.2456)*(1/factor)]
        if all(x <= 0 for x in [x_sheet, y_sheet]):
            msg = 'Dimension of graphene sheet must be greater than zero'
            raise ValueError(msg)
        self.name = 'C'
        carbon_locations = [[0, 0, 0], [2/3, 1/3, 0]]
        basis = {self.name: carbon_locations}
        graphene_lattice = mb.Lattice(lattice_spacing=[0.2456, 0.2456, 0.335], angles=[90, 90, 120], lattice_points=basis)
        carbon = mb.Compound(name=self.name)
        graphene = graphene_lattice.populate(compound_dict={self.name: carbon},
                                             x=replicate[0], y=replicate[1],
                                             z=self.sheets)
        for particle in graphene.particles():
            if particle.xyz[0][0] < 0:
                particle.xyz[0][0] += graphene.periodicity[0]
        self.graphene_dims = graphene.periodicity
        self.graphene_dims[1] *= factor  # cos(30)*.246
        bottom_sheet = mb.clone(graphene)
        bottom_sheet.translate([0, self.pore_width + (self.graphene_dims[2] - 0.335), 0])
        bottom_sheet.spin(1.5708, [1, 0, 0])
        top_sheet = mb.clone(graphene)
        top_sheet.spin(1.5708, [1, 0, 0])
        self.bot_xyz = bottom_sheet.xyz
        self.top_xyz = top_sheet.xyz
        system = mb.Compound()
        system.from_parmed(structure=bottom_sheet.to_parmed() + top_sheet.to_parmed())
        self.add(system)
