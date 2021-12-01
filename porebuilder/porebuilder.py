import mbuild as mb
import numpy as np

__all__ = ['GraphenePoreSolvent', 'GraphenePore', 'GrapheneSurface']


class GraphenePore(mb.Compound):
    """A general slit pore recipe.

    Parameters
    ----------
    pore_length : int, default=4
        dimensions of graphene sheet length in nm
    pore_depth : int, default=4
        dimensions of graphene sheet depth in nm
    n_sheets : int, default=3
        number of parallel graphene sheets
    pore_width: int, default=1
        width of slit pore in nm
    slit_pore_dim : int, default=1
        dimension slit pore, default is in the y-axis

    Attributes
    ----------
    see mbuild.Compound

    """
    def __init__(self, pore_length=4, pore_depth=3, n_sheets=3, pore_width=1, slit_pore_dim=1):
        super(GraphenePore, self).__init__()

        factor = np.cos(np.pi/6)
        # Estimate the number of lattice repeat units
        replicate = [int(pore_length/0.2456), (pore_depth/0.2456)*(1/factor)]
        if all(x <= 0 for x in [pore_length, pore_depth]):
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
                particle.xyz[0][0] += graphene.box.Lx

        bot_sheet = mb.clone(graphene)
        bot_sheet.name = 'BOT'
        top_sheet = mb.clone(graphene)
        if slit_pore_dim == 0:
            bot_sheet.spin(1.5708, [0, 1, 0])
            top_sheet.spin(1.5708, [0, 1, 0])
            top_sheet.translate([pore_width + (graphene.box.Lz - 0.335), 0, 0])
        elif slit_pore_dim == 1:
            bot_sheet.spin(1.5708, [1, 0, 0])
            top_sheet.spin(1.5708, [1, 0, 0])
            top_sheet.translate([0, pore_width + (graphene.box.Lz - 0.335), 0])
        elif slit_pore_dim == 2:
            top_sheet.translate([0, 0, pore_width + (graphene.box.Lz - 0.335)])
        top_sheet.name = 'TOP'
        self.add(top_sheet)
        self.add(bot_sheet)
        if slit_pore_dim == 0:
            new_Lx = 2 * graphene.box.Lz - lattice_spacing[2] + pore_width
            new_Ly =  graphene.box.Lx
            new_Lz = factor * graphene.box.Ly

            new_box = mb.Box((new_Lx,
                              new_Ly,
                              new_Lz),
                             (90, 90, 90)
                             )
            graphene.box = new_box

        elif slit_pore_dim == 1:
            new_Lx = graphene.box.Lx
            new_Ly =  2 * graphene.box.Lz - lattice_spacing[2] + pore_width
            new_Lz = factor * graphene.box.Ly

            new_box = mb.Box((new_Lx,
                              new_Ly,
                              new_Lz),
                             (90, 90, 90)
                             )
            graphene.box = new_box

        elif slit_pore_dim == 2:
            new_Lx = graphene.box.Lx
            new_Ly = factor * graphene.box.Ly
            new_Lz = 2 * graphene.box.Lz - lattice_spacing[2] + pore_width

            new_box = mb.Box((new_Lx,
                              new_Ly,
                              new_Lz),
                             (90, 90, 90)
                             )
            graphene.box = new_box

        self.xyz -= np.min(self.xyz, axis=0)



class GraphenePoreSolvent(mb.Compound):
    """A general slit pore recipe that includes baths of fluid.

    Parameters
    ----------
    pore_length : int, default=4
        dimensions of graphene sheet length in nm
    pore_depth : int, default=4
        dimensions of graphene sheet depth in nm
    n_sheets : int, default=3
        number of parallel graphene sheets
    pore_width: int, default=1
        width of slit pore in nm
    slit_pore_dim : int, default=1
        dimension slit pore, default is in the y-axis
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
    def __init__(self, pore_length=4, pore_depth=3, n_sheets=3, pore_width=1,
                 x_bulk=3, slit_pore_dim=1, solvent=None, n_solvent=100):

        super(GraphenePoreSolvent, self).__init__()

        pore = GraphenePore(pore_length=pore_length, pore_depth=pore_depth,
                            n_sheets=n_sheets, pore_width=pore_width,
                            slit_pore_dim=slit_pore_dim)

        box = mb.Box(lengths=[pore.box.Lx, pore.box.Ly, pore.box.Lz])
        if x_bulk != 0:
            box_max_0_direction = box.from_lo_hi_tilt_factors
            new_box = mb.Box((box_max_0_direction.Lx + 2 * x_bulk,
                              pore.box.Ly,
                              pore.box.Lz),
                             (90, 90, 90)
                             )
            box = new_box

        system = mb.solvate(pore, solvent, n_solvent, box=[pore.box.Lx, pore.box.Ly, pore.box.Lz], overlap=0.2)

        for child in system.children:
            self.add(mb.clone(child))

        # reset box dimenstions to box maxes
        box_max_all_direction = box.from_lo_hi_tilt_factors

        new_box = mb.Box((box_max_all_direction.box.Lx,
                          box_max_all_direction.box.Ly,
                          box_max_all_direction.box.Lz),
                         (90, 90, 90)
                         )
        system.box = new_box


class GrapheneSurface(mb.Compound):
    """A general graphene surface recipe exposed to vacuum.

    Parameters
    ----------
    x_length : int, default=4
        dimensions of graphene sheet length in nm
    y_length : int, default=4
        dimensions of graphene sheet depth in nm
    n_sheets : int, default=3
        number of parallel graphene sheets
    vacuum : float, default=10.0
        Dimension of vacuum space in z-direction in nm

    Attributes
    ----------
    see mbuild.Compound

    """
    def __init__(self, x_length=4, y_length=4, n_sheets=3, vacuum=10.0):
        super(GrapheneSurface, self).__init__()

        factor = np.cos(np.pi/6)
        # Estimate the number of lattice repeat units
        replicate = [int(x_length/0.2456), (y_length/0.2456)*(1/factor)]
        if all(x <= 0 for x in [x_length, y_length]):
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
                particle.xyz[0][0] += graphene.box.Lx

        new_Ly = graphene.box.Ly * factor
        new_Lz = graphene.box.Lz - lattice_spacing[2] + vacuum

        new_box = mb.Box((graphene.box.Lx,
                          new_Ly,
                          new_Lz),
                         (90, 90, 90)
                         )
        graphene.box = new_box

        self.add(graphene)
        self.xyz -= np.min(self.xyz, axis=0)
