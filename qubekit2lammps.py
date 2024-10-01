import argparse
import numpy as np
import networkx as nx
import sys

from PTWriter import FFWriter
from PTReader import OpenMMReader
from PTReader import POSCARReader
from PTCrystal import Crystal
from PTTools import ForceField
from PTTools import convert_qubekit2lammps


parser = argparse.ArgumentParser(description="QubeKit to LAMMPS force field converter")

parser.add_argument("-m","--molecules", help="Configurated molecules from QubeKit in .xml format", required=True, nargs='*')
parser.add_argument("-c","--crystal", help="POSCAR file", required=False)
args = parser.parse_args()

omms = [OpenMMReader(name) for name in args.molecules]

forcefields = [convert_qubekit2lammps(omm.get_forcefield()) for omm in omms]

if args.crystal is None:
    if len(args.molecules) > 1:
        sys.exit("Without a crystal, only one molecule can be configured at once")
    
    writer = FFWriter(args.molecules[0] + ".data", [omms[0].get_charges()], [omms[0].get_masses()], [[i for i in range(len(omms[0].get_charges()))]], forcefields, [np.zeros((len(omms[0].get_charges()), 3))], [1], [0])

    writer.writeLammps()
else:
    graphs = [nx.Graph() for i in args.molecules]

    [graph.add_nodes_from([(i, {"type": type}) for i, type in enumerate(forcefields[j].pairs)]) for j, graph in enumerate(graphs)]
    [graph.add_edges_from([(bond["atom1"], bond["atom2"]) for bond in forcefields[j].bonds]) for j, graph in enumerate(graphs)]

    poscar = POSCARReader(args.crystal)
    atoms = poscar.get_atoms()
    cell = poscar.get_lattice()
    
    crystal = Crystal(atoms, cell, poscar.get_types(), graphs)
    crystal.parametrize()

    ctypes = [(k, v) for k, v in crystal.get_crystal_atoms_types().items()]

    atoms_rearr = []
    types_rearr = []
    molecules_rearr = []
    charges = []
    masses = []

    #cmolecules = dict(sorted(crystal.get_crystal_molamounts().items()))
    cmolecules = crystal.get_crystal_molamounts()
    left = 0
    right = 0

    for graph, amount in cmolecules.items():

        for k in range(amount):
            right += len(graphs[graph].nodes)

            # this thing here takes molecules one by one from an array and sorts its atoms
            tempsorted = sorted(ctypes[left:right], key=lambda x: x[1])
            
            atoms_rearr.append([atoms[tempsorted[i][0]] for i in range(len(graphs[graph].nodes))])
            types_rearr.append([tempsorted[i][1] for i in range(len(graphs[graph].nodes))])
                
            molecules_rearr.append(graph)

            
            charges.append(omms[graph].get_charges())
            
            left += len(graphs[graph].nodes)

        masses.append(omms[graph].get_masses())


    writer = FFWriter(args.crystal + ".data", charges, masses, types_rearr, forcefields, atoms_rearr, cmolecules, molecules_rearr, cell)

    writer.writeLammps()
    