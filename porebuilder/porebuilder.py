import mbuild as mb
import numpy as np
from copy import deepcopy
from six import string_types

__all__ = ['GraphenePoreSolvent', 'GraphenePoreFunctionalized', 'GraphenePore']


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

    def __init__(self, pore_depth=4, side_dim=3, n_sheets=3, pore_width=1):
        super(GraphenePore, self).__init__()

        factor = np.cos(np.pi / 6)
        # Estimate the number of lattice repeat units
        replicate = [int(pore_depth / 0.2456),
                     (side_dim / 0.2456) * (1 / factor)]
        if all(x <= 0 for x in [pore_depth, side_dim]):
            msg = 'Dimension of graphene sheet must be greater than zero'
            raise ValueError(msg)
        carbon = mb.Compound()
        carbon.name = 'C'
        carbon_locations = [[0, 0, 0], [2 / 3, 1 / 3, 0]]
        basis = {carbon.name: carbon_locations}
        lattice_spacing = [0.2456, 0.2456, 0.335]
        angles = [90.0, 90.0, 120.0]

        graphene_lattice = mb.Lattice(lattice_spacing=lattice_spacing,
                                      angles=angles, lattice_points=basis)

        graphene = graphene_lattice.populate(
            compound_dict={
                carbon.name: carbon},
            x=replicate[0],
            y=replicate[1],
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
        top_sheet.translate(
            [0, pore_width + (graphene.periodicity[2] - 0.335), 0])
        top_sheet.name = 'TOP'
        self.add(top_sheet)
        self.add(bot_sheet)
        self.periodicity[0] = graphene.periodicity[0]
        self.periodicity[1] = 2 * graphene.periodicity[2] - \
            lattice_spacing[2] + pore_width
        self.periodicity[2] = graphene.periodicity[1]
        self.xyz -= np.min(self.xyz, axis=0)


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

    def __init__(self, pore_depth=4, side_dim=3, n_sheets=3, pore_width=1,
                 x_bulk=3, solvent=None, n_solvent=100):

        super(GraphenePoreSolvent, self).__init__()

        pore = GraphenePore(pore_depth=pore_depth, side_dim=side_dim,
                            n_sheets=n_sheets, pore_width=pore_width)

        box = mb.Box(pore.periodicity)
        box.maxs[0] += 2 * x_bulk

        system = mb.solvate(pore, solvent, n_solvent, box=box, overlap=0.2)

        for child in system.children:
            self.add(mb.clone(child))

        self.periodicity = box.maxs


class GraphenePoreFunctionalized(mb.Compound):
    """A general slit pore recipe that functionalizes the inner surfaces

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
    func_groups : mb.Compound or List-like thereof, defualt=None
        compounds to attach to the surface of a slit pore
    func_ports : string or List-like thereof, defualt='up'
        ports directly accesible by the func_groups that are the point of
         attachment to the pore
    func_percent : float or List-like thereof, defualt=.03
        percentage of surface to be functionalized by groups in func_groups
    NOTE: Both func_ports and func_percent must be in the same order as their
        corresponding functional group in func_groups. If only one value is given
        for either func_ports or func_percent that value will be used for all
        compounds in func_groups.
    Attributes
    ----------
    see mbuild.Compound

    """

    def __init__(self, pore_depth=4, side_dim=3, n_sheets=3, pore_width=1,
                 func_groups=None, func_percent=.03, func_ports='up'):

        super(GraphenePoreFunctionalized, self).__init__()

        pore = GraphenePore(pore_depth=pore_depth, side_dim=side_dim,
                            n_sheets=n_sheets, pore_width=pore_width)

        if isinstance(func_groups, mb.Compound):
            func_groups = [func_groups]
        if isinstance(func_percent, float):
            func_percent = [func_percent] * len(func_groups)
        if isinstance(func_ports, string_types):
            func_ports = [func_ports] * len(func_groups)

        if len(func_groups) != len(func_percent) or len(
                func_groups) != len(func_ports):
            raise ValueError(
                "If more than one port name or percent is given then "
                "it must be specifeid for all functional groups")

        cdf = deepcopy(func_percent)
        for i in range(0, len(func_percent)):
            if i != 0:
                cdf[i] = sum(func_percent[0:i + 1])

        Top = pore.children[0]
        Bot = pore.children[1]

        t_surface = []
        b_surface = []

        for C in Top:
            if C.pos[0] >= 0 and C.pos[1] <= .336 * \
                    (n_sheets - 1) + pore_width:
                t_surface.append(C)

        for C in t_surface:
            roll = np.random.rand()
            for chance, group, port in zip(cdf, func_groups, func_ports):
                if roll <= chance:
                    down_port = mb.Port(
                        anchor=C, orientation=[
                            0, -1, 0], separation=0.075)
                    C.add(down_port, 'down', containment=False)
                    new_group = deepcopy(group)
                    Top.add(new_group)
                    mb.force_overlap(
                        new_group,
                        new_group.labels[port],
                        C.labels['down'])
                    break

        for C in Bot:
            if C.pos[0] >= 0 and C.pos[1] >= .335 * (n_sheets - 1):
                b_surface.append(C)

        for C in b_surface:
            roll = np.random.rand()
            for chance, group, port in zip(cdf, func_groups, func_ports):
                if roll <= chance:
                    up_port = mb.Port(
                        anchor=C, orientation=[
                            0, 1, 0], separation=0.075)
                    C.add(up_port, 'up', containment=False)
                    new_group = deepcopy(group)
                    Bot.add(new_group)
                    mb.force_overlap(
                        new_group, new_group.labels[port], C.labels['up'])
                    break

        for child in pore.children:
            self.add(mb.clone(child))
