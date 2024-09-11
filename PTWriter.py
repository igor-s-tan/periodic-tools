import numpy as np
import PTTools

class FFWriter:
    def __init__(self, 
                 filename: str,
                 charges: list,
                 masses: list,
                 types: list,
                 forcefields: list,
                 atoms: list,
                 molamounts: dict,
                 molecules: list,
                 cell=None,
                 virtual_sites=None):
                     
        self.filename = filename
        
        self.charges = charges
        self.masses = masses
        self.types = types
        
        self.forcefields = forcefields
        
        self.atoms = atoms
        self.molamounts = molamounts
        self.molecules = molecules
        self.cell = cell
        
        # NOT IMPLEMENTED YET
        self.virtual_sites = virtual_sites

    def writeLammps(self):
        
        # Calculation of all possible amounts
        atoms_numbers, bonds_numbers, angles_numbers, dihedrals_numbers, impropers_numbers = 0, 0, 0, 0, 0
        pairs_types, bonds_types, angles_types, dihedrals_types, impropers_types = 0, 0, 0, 0, 0
        
        for cluster in self.atoms:
            atoms_numbers += len(cluster)

        for i, ff in enumerate(self.forcefields):
            bonds_numbers += len(ff.bonds) * self.molamounts[i]
            angles_numbers += len(ff.angles) * self.molamounts[i]
            dihedrals_numbers += len(ff.dihedrals) * self.molamounts[i]
            impropers_numbers += len(ff.impropers) * self.molamounts[i]
            
            pairs_types += len(ff.pairs)
            bonds_types += len(ff.bonds)
            angles_types += len(ff.angles)
            dihedrals_types += len(ff.dihedrals)
            impropers_types += len(ff.impropers)

        # Use default vectors when cell is not passed
        if self.cell is None:
            self.cell = PTTools.NullLMPCell
        else:
            self.cell = PTTools.vec2lammps(self.cell)
        
        with open(self.filename[:self.filename.find('.')] + '.data', 'w', encoding='utf-8') as file:
            result = f'Periodic-Tools output\n\n'

            result += f'{atoms_numbers} atoms\n' \
                      f'{bonds_numbers} bonds\n' \
                      f'{angles_numbers} angles\n' \
                      f'{dihedrals_numbers} dihedrals \n' \
                      f'{impropers_numbers} impropers\n\n' \
                      f'{pairs_types} atom types\n' \
                      f'{bonds_types} bond types\n' \
                      f'{angles_types} angle types\n' \
                      f'{dihedrals_types} dihedral types\n' \
                      f'{impropers_types} improper types\n\n'

            result += f'{self.cell[0][0]} {self.cell[0][1]} xlo xhi\n' \
                      f'{self.cell[1][0]} {self.cell[1][1]} ylo yhi\n' \
                      f'{self.cell[2][0]} {self.cell[2][1]} zlo zhi\n' \
                      f'{self.cell[3][0]} {self.cell[3][1]} {self.cell[3][2]} xy xz yz\n\n'

            result += f'Masses\n\n'
            count = 1
            for cluster in self.masses:
                for mass in cluster:
                    result += f'{count} {mass}\n'
                    count += 1

            result += f'\nPair Coeffs\n\n'
            count = 1
            for ff in self.forcefields:
                for atom_type in ff.pairs:
                    result += f'{count} {atom_type["epsilon"]} {atom_type["sigma"]}\n'
                    count += 1
                    
            result += f'\nBond Coeffs\n\n'
            count = 1
            for ff in self.forcefields:
                for bond_type in ff.bonds:
                    result += f'{count} {bond_type["k"]} {bond_type["length"]}\n'
                    count += 1

            result += f'\nAngle Coeffs\n\n'
            count = 1
            for ff in self.forcefields:
                for angle_type in ff.angles:
                    result += f'{count} {angle_type["k"]} {angle_type["angle"]}\n'
                    count += 1

            result += f'\nDihedral Coeffs\n\n'
            count = 1
            for ff in self.forcefields:
                for dihedral_type in ff.dihedrals:
                    result += f'{count} {dihedral_type["k1"]} {dihedral_type["k2"]} ' \
                              f'{dihedral_type["k3"]} {dihedral_type["k4"]}\n'
                    count += 1
            
            result += f'\nImproper Coeffs\n\n'
            count = 1
            for ff in self.forcefields:
                for improper_type in ff.impropers:
                    result += f'{count} {improper_type["k2"]} -1 2\n'
                    count += 1

            result += f'\nAtoms\n\n'
            count = 1
            for i, an in enumerate(self.atoms):
                for j, atom in enumerate(an):
                    result += f'{count} ' \
                              f'{i + 1} ' \
                              f'{self.types[i][j] + 1} ' \
                              f'{self.charges[i][j]} ' \
                              f'{atom[0]} ' \
                              f'{atom[1]} ' \
                              f'{atom[2]}\n'
                    count += 1

            result += f'\nBonds\n\n'
            count = 1
            atomshift = 0
            for mn in self.molecules:
                typeshift = sum([len(self.forcefields[i].bonds) for i in range(mn)])
                for j, bond in enumerate(self.forcefields[mn].bonds):
                    result += f'{count} ' \
                              f'{j + 1 + typeshift} ' \
                              f'{bond["atom1"] + 1 + atomshift} ' \
                              f'{bond["atom2"] + 1 + atomshift}\n'
                    count += 1
                atomshift += len(self.forcefields[mn].pairs)

            result += f'\nAngles\n\n'
            count = 1
            atomshift = 0
            for mn in self.molecules:
                typeshift = sum([len(self.forcefields[i].angles) for i in range(mn)])
                for j, angle in enumerate(self.forcefields[mn].angles):
                    result += f'{count} ' \
                              f'{j + 1 + typeshift} ' \
                              f'{angle["atom1"] + 1 + atomshift} ' \
                              f'{angle["atom2"] + 1 + atomshift} ' \
                              f'{angle["atom3"] + 1 + atomshift}\n'
                    count += 1
                atomshift += len(self.forcefields[mn].pairs)

            result += f'\nDihedrals\n\n'
            count = 1
            atomshift = 0
            for mn in self.molecules:
                typeshift = sum([len(self.forcefields[i].dihedrals) for i in range(mn)])
                for j, dihedral in enumerate(self.forcefields[mn].dihedrals):
                    result += f'{count} ' \
                              f'{j + 1 + typeshift} ' \
                              f'{dihedral["atom1"] + 1 + atomshift} ' \
                              f'{dihedral["atom2"] + 1 + atomshift} ' \
                              f'{dihedral["atom3"] + 1 + atomshift} ' \
                              f'{dihedral["atom4"] + 1 + atomshift}\n'
                    count += 1
                atomshift += len(self.forcefields[mn].pairs)

            result += f'\nImpropers\n\n'
            count = 1
            atomshift = 0
            for mn in self.molecules:
                typeshift = sum([len(self.forcefields[i].impropers) for i in range(mn)])
                for j, improper in enumerate(self.forcefields[mn].impropers):
                    result += f'{count} ' \
                              f'{j + 1 + typeshift} ' \
                              f'{improper["atom1"] + 1 + atomshift} ' \
                              f'{improper["atom2"] + 1 + atomshift} ' \
                              f'{improper["atom3"] + 1 + atomshift} ' \
                              f'{improper["atom4"] + 1 + atomshift}\n'
                    count += 1
                atomshift += len(self.forcefields[mn].pairs)

            file.write(result)
