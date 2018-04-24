import mbuild as mb
import numpy as np

__all__ = ['GraphenePoreSolvent', 'GraphenePore']


class GraphenePoreSolvent(mb.Compound):
    """A general slit pore recipe that includes baths of fluid.

    Parameters
    ----------
    pore_depth : int, default=4
        dimensions of graphene sheet in x direction in nm
    side_dim : int, default=4
        dimensions of graphene sheet in z direction in nm
    n_sheets : int, default=3
        number of parallel graphene sheets
    pore_width: int, default=1
        width of slit pore in nm
    x_bulk : int, default=3
        length of bulk region in x-direction in nm
    solvent : list of mbuild.Compound
        list of compound(s) to solvate the system with
    n_solvent : list of int
        number of solvents to solvate the system with
    NOTE: length of `solvent` must match length of `n_solvent`

    Attributes
    ----------
    see mbuild.Compound

    """
    def __init__(self, pore_depth, side_dim, n_sheets, pore_width, x_bulk,
            solvent, n_solvent):
        super(GraphenePoreSolvent, self).__init__()

        pore = GraphenePore(pore_depth=pore_depth, side_dim=side_dim,
                              n_sheets=n_sheets, pore_width=pore_width)

        box = mb.Box(pore.periodicity)
        box.maxs[0] += 2 * x_bulk

        system = mb.solvate(pore, solvent, n_solvent, box=box, overlap=0.2)

        for child in system.children:
            self.add(mb.clone(child))

        self.periodicity = box.maxs


class GraphenePore(mb.Compound):
    """A general slit pore recipe.

    Parameters
    ----------
    pore_depth : int, default=4
        dimensions of graphene sheet in x direction in nm
    side_dim : int, default=4
        dimensions of graphene sheet in z direction in nm
    n_sheets : int, default=3
        number of parallel graphene sheets
    pore_width: int, default=1
        width of slit pore in nm

    Attributes
    ----------
    see mbuild.Compound

    """
    def __init__(self, pore_depth, side_dim, n_sheets, pore_width):
        super(GraphenePore, self).__init__()

        # Do some math to figure out how much to replicate graphene cell.
        # Multiply replicate[1] by 15/13 to take into account later
        # multiplication
        # TODO: Figure out if rounding is necessary
        factor = np.cos(np.pi/6)
        replicate = [(pore_depth/0.2456),
                (side_dim/0.2456)*(1/factor)]
        if all(x <= 0 for x in [pore_depth, side_dim]):
            msg = 'Dimension of graphene sheet must be greater than zero'
            raise ValueError(msg)
        carbon = mb.Compound()
        carbon.name = 'C'
        carbon_locations = [[0, 0, 0], [2/3, 1/3, 0]]
        basis = {carbon.name: carbon_locations}
        lattice_spacing = [0.2456, 0.2456, 0.335]
        angles = [90.0, 90.0, 120.0]

        graphene_lattice = mb.Lattice(lattice_spacing=lattice_spacing,
                                      angles=angles, lattice_points=basis)

        graphene = graphene_lattice.populate(compound_dict={carbon.name: carbon},
                                             x=replicate[0], y=replicate[1],
                                             z=n_sheets)

        for particle in graphene.particles():
            if particle.xyz[0][0] < 0:
                particle.xyz[0][0] += graphene.periodicity[0]
        graphene.periodicity[1] *= factor  # cos(30)*.246
        bot_sheet = mb.clone(graphene)
        bot_sheet.spin(1.5708, [1, 0, 0])
        bot_sheet.name = 'BOT'
        top_sheet = mb.clone(graphene)
        top_sheet.spin(1.5708, [1, 0, 0])
        top_sheet.translate([0, pore_width + (graphene.periodicity[2] - 0.335), 0])
        top_sheet.name = 'TOP'
        self.add(top_sheet)
        self.add(bot_sheet)
        self.periodicity[0] = graphene.periodicity[0]
        self.periodicity[1] = 2 * graphene.periodicity[2] - lattice_spacing[2] + pore_width
        self.periodicity[2] = graphene.periodicity[1]
        self.xyz -= np.min(self.xyz, axis=0)
