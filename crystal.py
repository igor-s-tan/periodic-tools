import numpy as np
import networkx as nx

COVALENT_RADII = {1: 0.37, 'H': 0.37,
                  6: 0.77, 'C': 0.77,
                  7: 0.75, 'N': 0.75,
                  8: 0.73, 'O': 0.73}


class Crystal:
    
    def __init__(self, atoms: list, cell: list, types: list, graph: nx.classes.graph.Graph):
        self.atoms = [np.array(atom) for atom in atoms]
        self.cell = np.array(cell)
        self.types = types
        self.natoms = len(atoms)
        self.graph = graph
        
        self.natoms_in_molecule = len(list(self.graph.nodes))
        self.nmolecules = int(self.natoms / self.natoms_in_molecule)
        
        self.crystal_atoms_types = dict()
        self.crystal_molecules = []


    def parametrize(self) -> None:
        # expanding cell to eliminate fragments
        self.__make_super_cell()
        print("Number of molecules: " + str(self.nmolecules))
        big_graph = nx.Graph()
        print("Searching expanded graph...")
        big_graph.add_nodes_from([(i, {"type": self.types[i]}) for i in range(len(self.atoms))])
        for i in range(len(self.atoms)):
            for j in range(i+1, len(self.atoms)):
                if np.linalg.norm(self.atoms[i] -
                                  self.atoms[j]) < COVALENT_RADII[self.types[i]] + COVALENT_RADII[self.types[j]]:
                    big_graph.add_edge(i, j)
        fragments = []
        count = 0
        viewed = set()
        for i in range(self.natoms):
            if i in viewed:
                continue
            molecule = nx.node_connected_component(big_graph, i)
            
            viewed.update(molecule)
            isodict = nx.vf2pp_isomorphism(big_graph.subgraph(molecule), self.graph)
            for k in isodict:
                self.crystal_atoms_types[k - (k // self.natoms * self.natoms)] = isodict[k]
            if len(self.crystal_atoms_types) == self.natoms:
                break
            
        # for frag in fragments:
            # for atom in frag:
                # print(np.matmul(self.cell, self.atoms[atom]))
        

    def __make_super_cell(self) -> None:
        # one vector
        for k in [-1, 1]:
            for j in range(3):
                for i in range(self.natoms):
                    self.atoms.append(np.add(self.atoms[i], k * self.cell[j]))

        # two vectors
        for k in [-1, 1]:
            for h in [-1, 1]:
                for j in range(3):
                    for i in range(self.natoms):
                        self.atoms.append(np.add(self.atoms[i], k * np.add(self.cell[j], h * self.cell[(j + 1) % 3])))

        # three vectors
        for k in [-1, 1]:
            for h in [-1, 1]:
                for j in [-1, 1]:
                    for i in range(self.natoms):
                        self.atoms.append(
                            np.add(self.atoms[i], k * np.add(self.cell[0], h * np.add(self.cell[1], j * self.cell[2]))))

        self.types *= 27



    def get_crystal_atoms_types(self) -> dict:
        return self.crystal_atoms_types
        
    def get_crystal_molecules(self) -> list:
        return self.crystal_molecules
        