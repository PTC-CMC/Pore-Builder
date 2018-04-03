import itertools as it
import mbuild as mb
import numpy as np
from mbuild.compound import Compound
from mbuild.coordinate_transform import force_overlap
from mbuild.utils.validation import assert_port_exists
from mbuild import clone
from copy import deepcopy
import math

__all__ = ['Pores']

class Pores(mb.Compound):
    """A general slit pore recipe
    
    Parameters
    ----------
    x_sheet: int
        dimensions of graphene sheet in x-direction [nm]
    y_sheet: int
        dimensions of graphene sheet in y-direction [nm]
    sheets: int
        number of graphene sheets, default=3
    pore_width: int
        width of slit pore [nm]
    x_bulk: int
        length of bulk region in x-direction [nm]
    solvent: dict {'name': 'solvent file'}, optional, Default=None
        compound to solvate the system with.  If not provided, system will not be solvated
    n_solvent: int, optional, Default=None
        number of solvents to solvate the system with.  If not
        provided, system will not be solvated
    Attributes
    ----------
    
    Notes: Match graphene y-dimension with box x-dimension
    """
    def __init__(self,x_sheet, y_sheet, sheets, pore_width, x_bulk,
            solvent=None, n_solvent=None):
        super(Pores,self).__init__()
        self.x_sheet = x_sheet
        self.y_sheet = y_sheet
        self.sheets = sheets
        self.pore_width = pore_width
        self.x_bulk = x_bulk
        
        # Do some math to figure out how much to replicate graphene cell. TODO: Figure out if rounding is necessary
        # Multiply replicate[1] by 15/13 to take into account later multiplication
        factor = np.cos(math.pi/6)
        print(factor)
        replicate = [(self.x_sheet/0.2456),
                (self.y_sheet/0.2456)*(1/factor)]
        if all(x <= 0 for x in [x_sheet, y_sheet]):
            raise ValueError('Dimension of graphene sheet must be greater than zero')
        self.name = 'C'
        carbon_locations = [[0,0,0], [2/3,1/3,0]]
        basis = {self.name: carbon_locations}
        graphene_lattice = mb.Lattice(lattice_spacing=[.2456,.2456,.335], angles=[90,90,120], lattice_points=basis)
        carbon = mb.Compound(name=self.name)
        graphene = graphene_lattice.populate(compound_dict={self.name: carbon},
                                         x=replicate[0],y=replicate[1],z=self.sheets)
        for particle in graphene.particles():
            if particle.xyz[0][0] < 0:
                particle.xyz[0][0] += graphene.periodicity[0]
        self.graphene_dims = graphene.periodicity
        self.graphene_dims[1] *= factor # cos(30)*.246
        bottom_sheet = mb.clone(graphene)
        bottom_sheet.translate([0, self.pore_width+(self.graphene_dims[2]-.335), 0]) 
        bottom_sheet.spin(1.5708,[1,0,0])
        top_sheet = mb.clone(graphene)
        top_sheet.spin(1.5708,[1,0,0])
        self.bot_xyz = bottom_sheet.xyz
        self.top_xyz = top_sheet.xyz
        system = mb.Compound()
        system.from_parmed(structure=bottom_sheet.to_parmed() + top_sheet.to_parmed())
        if solvent:
            self._solvate(solvent=solvent, n_solvent=n_solvent, system=system)
        else:
            self.add(system)
    def _solvate(self, system, solvent, n_solvent):
        """Solvate slit pore box
        Parameters
        ----------
        solvent: compatible file 
        n_solvent: int
            Number of compounds to solvate with
        """
        if len(solvent) == 1:
            fluid = mb.load(solvent)
            fluid.name = '{}'.format(solvent)
        elif len(solvent) == 2:
            fluid_1 = mb.load(solvent[0])
            fluid_1.name = '{}'.format(solvent[0])
            fluid_2 = mb.load(solvent[1])
            fluid_2.name = '{}'.format(solvent[1])
            fluid = [fluid_1, fluid_2]
        elif len(solvent) > 2:
            raise ValueError('"gph_pore_solv" class currently only supports a maximum of 2 solvents')
        box = [(self.x_bulk*2)+self.graphene_dims[0]+.5,
                self.pore_width+(2*(self.graphene_dims[2]-0.335)+.5),
                self.graphene_dims[1]]
        print(box)
        system = mb.solvate(system, fluid, n_solvent, box=box, overlap=0.2)
        system.periodicity = box
        self.add(system)
