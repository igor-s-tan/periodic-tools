import numpy as np
import networkx as nx
import sys

import PTTools


class Crystal:
    
    def __init__(self, atoms: list, cell: list, types: list, graphs: list):
        self.atoms = atoms.copy()
        self.cell = cell
        self.types = types
        self.natoms = len(atoms)
        self.graphs = graphs
        
        # atom_id-atom_type dictionary
        self.crystal_atoms_types = dict()
        
        # molecule_number-molecule_amount dictionary
        self.molamounts = dict()


    def parametrize(self):
        
        # expanding the cell to eliminate fragments
        self.__make_super_cell()
        
        # scanning the cell for molecules
        big_graph = nx.Graph()
        big_graph.add_nodes_from([(i, {"type": self.types[i]}) for i in range(len(self.atoms))])
        for i in range(len(self.atoms)):
            for j in range(i+1, len(self.atoms)):
                if np.linalg.norm(self.atoms[i] -
                                  self.atoms[j]) < PTTools.COVALENT_RADII[self.types[i]] + PTTools.COVALENT_RADII[self.types[j]]:
                    big_graph.add_edge(i, j)
        
        # making graphs
        viewed = set()
        for i in range(self.natoms):
            # skip if current node has already been viewed
            if i in viewed:
                continue
                
            molecule = nx.node_connected_component(big_graph, i)
            
            # update viewed nodes set
            viewed.update(molecule)
            
            # making table of shifts to assign types correctly
            table = []
            shift = 0
            for graph in self.graphs:
                table.append(shift)
                shift += len(graph.nodes)
            
            # solving isomorphism problem
            for k, graph in enumerate(self.graphs):
                isodict = nx.vf2pp_isomorphism(big_graph.subgraph(molecule), graph)
                if isodict is not None:
                    if k not in self.molamounts:
                        self.molamounts[k] = 0
                    for key in isodict:
                        if key < self.natoms:
                            self.molamounts[k] += 1 
                        isodict[key] += table[k]
                    break
            
            # in case we caught a fragment on the egde of the cell
            if isodict is None:
                continue
              
            for k in isodict:
                self.crystal_atoms_types[k - (k // self.natoms * self.natoms)] = isodict[k]            
            
            # stop if full atoms pack has already been created
            if len(self.crystal_atoms_types) >= self.natoms * 2:
                break
        
    def __make_super_cell(self):
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
        
    def get_crystal_molamounts(self) -> dict:
        return {key: self.molamounts[key] // len(self.graphs[key].nodes) for key in self.molamounts}
        