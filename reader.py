from bs4 import BeautifulSoup


class OMMReader:

    def __init__(self, filename):
        with open(filename, 'r', encoding='utf-8') as file:
            data = file.read()

        self.bs = BeautifulSoup(data, 'xml')

        self.atom_types = self.get_types()

    # [type]
    def get_types(self) -> list:
        types_data = self.bs.select('AtomTypes > Type')
        types_list = [atom_type.get('element') for atom_type in types_data]
        return types_list

    # [mass]
    def get_masses(self) -> list:
        types_data = self.bs.select('AtomTypes > Type')
        masses_list = [float(atom_type.get('mass')) for atom_type in types_data]
        return masses_list

    # [charge]
    def get_charges(self) -> list:
        atoms_data = self.bs.select('NonbondedForce > Atom')
        charges_list = [float(atom_type.get('charge')) for atom_type in atoms_data]
        return charges_list

    # [sigma, epsilon]
    def get_pairs(self) -> list:
        atoms_data = self.bs.select('NonbondedForce > Atom')
        pairs_list = [{'epsilon': float(atom.get('epsilon')),
                       'sigma': float(atom.get('sigma'))} for atom in atoms_data]
        return pairs_list

    # [atom1, atom2, length, k]
    def get_bonds(self) -> list:
        bonds_data = self.bs.select('HarmonicBondForce > Bond')
        bonds_list = [{'atom1': int(bond.get('class1')),
                       'atom2': int(bond.get('class2')),
                       'length': float(bond.get('length')),
                       'k': float(bond.get('k'))} for bond in bonds_data]
        return bonds_list

    # [atom1, atom2, atom3, angle, k]
    def get_angles(self) -> list:
        angles_data = self.bs.select('HarmonicAngleForce > Angle')
        angles_list = [{'atom1': int(angle.get('class1')),
                        'atom2': int(angle.get('class2')),
                        'atom3': int(angle.get('class3')),
                        'angle': float(angle.get('angle')),
                        'k': float(angle.get('k'))} for angle in angles_data]
        return angles_list

    # [atom1, atom2, atom3, atom4, k1, k2, k3, k4, phase1, phase2, phase3, phase4]
    def get_dihedrals(self) -> list:
        dihedrals_data = self.bs.select('PeriodicTorsionForce > Proper')
        dihedrals_list = [{'atom1': int(dihedral.get('class1')),
                           'atom2': int(dihedral.get('class2')),
                           'atom3': int(dihedral.get('class3')),
                           'atom4': int(dihedral.get('class4')),
                           'k1': float(dihedral.get('k1')),
                           'k2': float(dihedral.get('k2')),
                           'k3': float(dihedral.get('k3')),
                           'k4': float(dihedral.get('k4')),
                           'phase1': float(dihedral.get('phase1')),
                           'phase2': float(dihedral.get('phase2')),
                           'phase3': float(dihedral.get('phase3')),
                           'phase4': float(dihedral.get('phase4'))} for dihedral in dihedrals_data]
        return dihedrals_list

    # [atom1, atom2, atom3, atom4, k1, k2, k3, k4, phase1, phase2, phase3, phase4]
    def get_impropers(self) -> list:
        impropers_data = self.bs.select('PeriodicTorsionForce > Improper')
        impropers_list = [{'atom1': int(improper.get('class1')),
                           'atom2': int(improper.get('class2')),
                           'atom3': int(improper.get('class3')),
                           'atom4': int(improper.get('class4')),
                           'k1': float(improper.get('k1')),
                           'k2': float(improper.get('k2')),
                           'k3': float(improper.get('k3')),
                           'k4': float(improper.get('k4')),
                           'phase1': float(improper.get('phase1')),
                           'phase2': float(improper.get('phase2')),
                           'phase3': float(improper.get('phase3')),
                           'phase4': float(improper.get('phase4'))} for improper in impropers_data]
        return impropers_list



class POSCARReader:
    
    def __init__(self, filename):
        self.lattice = []
        self.atoms = []
        self.types = []
        self.__intypes = []
        self.__inamounts = []
        self.is_fractional = False
        with open(filename, 'r', encoding='utf-8') as file:
            count = 0
            typecount = 0
            typeswitch = 0
            for line in file:
                if count <= 1:
                    count += 1
                elif count >= 2 and count <= 4:
                    self.lattice.append([float(i) for i in line.split()])
                    count += 1
                elif count == 5:
                    self.__intypes = line.split()
                    count += 1
                elif count == 6:
                    self.__inamounts = [int(i) for i in line.split()]
                    count += 1
                elif count == 7:
                    if "Fractional" in line:
                        self.is_fractional = True
                    count += 1
                elif count >= 7:
                    self.atoms.append([float(i) for i in line.split()])
                    self.types.append(self.__intypes[typeswitch])
                    typecount += 1
                    if self.__inamounts[typeswitch] == typecount:
                        typeswitch += 1
                        typecount = 0
                    count += 1
    
    
    def get_types(self) -> list:
        return self.types
    
    
    def get_atoms(self) -> list:
        return self.atoms
        
    
    def get_lattice(self) -> list:
        return self.lattice
        