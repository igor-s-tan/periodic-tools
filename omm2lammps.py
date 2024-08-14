import argparse
import numpy as np

from writer import FFWriter
from reader import OMMReader
from reader import POSCARReader
from crystal import Crystal
import networkx as nx


def omm2lammps(pairs, bonds, angles, dihedrals, impropers):
    for i in range(len(pairs)):
        pairs[i]['epsilon'] /= 4.18
        pairs[i]['sigma'] *= 10
    for i in range(len(bonds)):
        bonds[i]['length'] *= 10
        bonds[i]['k'] /= 100 * 4.184 * 2
    for i in range(len(angles)):
        angles[i]['angle'] *= 180 / np.pi
        angles[i]['k'] /= 4.184 * 2
    for i in range(len(dihedrals)):
        dihedrals[i]['k1'] /= 4.184 / 2
        dihedrals[i]['k2'] /= 4.184 / 2
        dihedrals[i]['k3'] /= 4.184 / 2
        dihedrals[i]['k4'] /= 4.184 / 2
    for i in range(len(impropers)):
        impropers[i]['k2'] /= 4.184

    return pairs, bonds, angles, dihedrals, impropers


omm = OMMReader("melamine.xml")

pairs, bonds, angles, dihedrals, impropers = omm2lammps(omm.get_pairs(), 
                                                        omm.get_bonds(),
                                                        omm.get_angles(),
                                                        omm.get_dihedrals(), 
                                                        omm.get_impropers())

types = omm.get_types()

molecule = nx.Graph()

molecule.add_nodes_from([(i, {"type": types[i]}) for i in range(len(pairs))])
molecule.add_edges_from([(bond["atom1"], bond["atom2"]) for bond in bonds])

poscar = POSCARReader("TEST.vasp")
atoms = poscar.get_atoms()
cell = poscar.get_lattice()

crystal = Crystal(atoms, cell, poscar.get_types(), molecule)
crystal.parametrize()

ctypes = [(k, v) for k, v in crystal.get_crystal_atoms_types().items()]

atoms_rearr = []

for mol in range(len(atoms) // len(molecule.nodes)):
    tempsorted = sorted(ctypes[mol*len(molecule.nodes):(mol+1)*len(molecule.nodes)], key=lambda x: x[1])
    atoms_rearr += [atoms[tempsorted[i][0]] for i in range(len(molecule.nodes))]


ff = FFWriter("melamine.data", omm.get_charges(), omm.get_masses(), pairs, bonds, angles, dihedrals, impropers, atoms_rearr, cell)

ff.writeLammps()

