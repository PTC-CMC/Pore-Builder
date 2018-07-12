from porebuilder import GraphenePoreFunctionalized
import mbuild as mb


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


class OH(mb.Compound):
    def __init__(self):
        super(OH, self).__init__()

        self.add(O())
        self.add(H())
        mb.force_overlap(
            self.children[1],
            self.children[1].labels['up'],
            self.children[0].labels['up'])
        up_port = mb.Port(
            anchor=self[0], orientation=[
                0, -1, 0], separation=.075)
        self.add(up_port, "up")


system = GraphenePoreFunctionalized(
    func_groups=[H(), O(), OH()], func_percent=[.1, .05, .03], func_ports='up')
    
system.save("func_dry_pore.mol2", overwrite=True)
