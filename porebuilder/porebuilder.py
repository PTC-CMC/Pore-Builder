import mbuild as mb
import numpy as np
from random import shuffle
from six import string_types

__all__ = ['GraphenePoreSolvent', 'GraphenePore']


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

        factor = np.cos(np.pi/6)
        # Estimate the number of lattice repeat units
        replicate = [int(pore_depth/0.2456), (side_dim/0.2456)*(1/factor)]
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
    func_groups : mb.Compound or List-like thereof, default=None
        compounds to attach to the surface of a slit pore
    func_ports : string or List-like thereof, default='up'
        ports directly accesible by the func_groups that are the point of
         attachment to the pore
    func_percent : float or List-like thereof, default=.03
        percentage of surface to be functionalized by groups in func_groups
    NOTE: Both func_ports and func_percent must be in the same order as their
        corresponding functional group in func_groups. If only one value is
        given for either func_ports or func_percent that value will be used for all
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

        if not len(func_groups) == len(func_percent) == len(func_ports):
            raise ValueError(
                "If more than one port name or percent is given then "
                "it must be specifeid for all functional groups")

        top = pore.children[0]
        bot = pore.children[1]

        t_surface = []
        b_surface = []

        size = len(top.xyz.T[1])

        #These two for loops selected the inner surfaces of the graphene pore
        #by searching through the array of all y positions. 

        for pos, i in zip(top.xyz.T[1], range(0,size)):
            if pos <= .336 * (n_sheets - 1) + pore_width:
                t_surface.append(top.children[i])
        
        for pos, i in zip(bot.xyz.T[1], range(0,size)):
            if pos >= .335 * (n_sheets - 1):
                b_surface.append(bot.children[i])

        #The pore is then functionalized randomly, the rng being taken care of 
        #by random.shuffle. THis method should get as close to the desired 
        #percent functionalization as possible.

        for pore_wall, surface, orientation_factor in zip((top,bot),(t_surface,
        b_surface),(-1,1)):
        
            shuffle(surface)

            #The queue is a list containing the amount of each functional group 
            #to be added to the pore's surface. The start list will contain 0, 
            #the first index to work with and the numbers between function 
            #groups. i.e. with a queue of [5,17,14] the start list should be 
            #[0,5,22]

            queue = np.multiply(np.array(func_percent),(len(surface)))
            queue = queue.astype(int)
            start = [0]
            for i in range(len(queue)-1):
                start.append(sum(start) + queue[i])

            for prev, n, group, port in zip(start, queue, func_groups, 
            func_ports):                   
                for i in range (prev,n+prev):
                    down_port = mb.Port(anchor=surface[i],orientation=[0, 
                    orientation_factor, 0], separation=0.075)
                    surface[i].add(down_port, 'down', containment=False)
                    new_group = mb.clone(group)
                    pore_wall.add(new_group)
                    mb.force_overlap(new_group, new_group.labels[port], 
                    surface[i].labels['down'])

        for child in pore.children:
            self.add(mb.clone(child))