import numpy as np


""" Wrapper structure for force field parameters """
class ForceField:
    
    def __init__(self, pairs=None, bonds=None, angles=None, dihedrals=None, impropers=None):
        self._pairs = self.list_init(pairs)                 
        self._bonds = self.list_init(bonds)                 
        self._angles = self.list_init(angles)               
        self._dihedrals = self.list_init(dihedrals)         
        self._impropers = self.list_init(impropers)         
    
    def __add__(self, ff):
        return ForceField(self._pairs + ff._pairs, 
                          self._bonds + ff._bonds, 
                          self._angles + ff._angles, 
                          self._dihedrals + ff._dihedrals, 
                          self._impropers + ff._impropers)

    def __iadd__(self, ff):
        self = self + ff
        return self
    
    def list_init(self, lst):
        if lst is None:
            return []
        else:
            return lst
    
    @property    
    def pairs(self):
        return self._pairs
    
    @property    
    def bonds(self):
        return self._bonds
        
    @property    
    def angles(self):
        return self._angles
        
    @property    
    def dihedrals(self):
        return self._dihedrals
        
    @property    
    def impropers(self):
        return self._impropers


""" Default cell vectors for LAMMPS """
NullLMPCell = np.array([[.0, 20.0],
                      [.0, 20.0],
                      [.0, 20.0],
                      [.0, .0, .0]], dtype=object)



""" Covalent radii dictionary """
COVALENT_RADII = {1: 0.38, 'H': 0.38,
                  6: 0.77, 'C': 0.77,
                  7: 0.75, 'N': 0.75,
                  8: 0.73, 'O': 0.73}

""" Energy conversion constant """
CAL2J = 4.184

""" Cell vectors to parameters converter """
def vec2param(cell: np.array) -> np.array:
    
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    
    alpha = np.arccos(np.dot(cell[0], cell[1]) / a / b)
    beta = np.arccos(np.dot(cell[0], cell[2]) / a / c)
    gamma = np.arccos(np.dot(cell[1], cell[2]) / c / b)
    
    return np.array([a, b, c, alpha, beta, gamma])


""" Cell parameters to vectors converter """
def param2vec(param: np.array) -> np.array:
    
    gammac = np.cos(param[5])
    gammas = np.sin(param[5])
    
    cx = param[2] * np.cos(param[4])
    cy = (param[2] * np.cos(param[3]) - cx * gammac) / gammas
        
    cell = np.array([[param[0], 0.0, 0.0],
                     [param[1] * gammac, param[1] * gammas, 0.0],
                     [cx, cy, np.sqrt(param[2] ** 2 - cx ** 2 - cy ** 2)]])
                     
    return cell
    
    
""" Cell parameters to LAMMPS """ 
def param2lammps(param: np.array) -> np.array:
    
    lx = param[0]
    xy = param[1] * np.cos(param[5])
    xz = param[2] * np.cos(param[4])
    ly = np.sqrt(param[1] ** 2 - xy ** 2)
    yz = (param[1] * param[2] * np.cos(param[3]) - xy * xz) / ly
    lz = np.sqrt(param[2] ** 2 - xz ** 2 - yz ** 2)
    
    return np.array([[0, lx],
                      [0, ly],
                      [0, lz],
                      [xy, xz, yz]], dtype=object)


  
""" Cell vectors to LAMMPS """ 
def vec2lammps(cell: np.array) -> np.array:
    
    return param2lammps(vec2param(cell))
    

""" Convert FF units from QubeKit(OpenMM) to LAMMPS """
def convert_qubekit2lammps(ff: ForceField):
    for i in ff.pairs:
        i['epsilon'] /= CAL2J
        i['sigma'] *= 10
    for i in ff.bonds:
        i['length'] *= 10
        i['k'] /= 100 * CAL2J * 2
    for i in ff.angles:
        i['angle'] *= 180 / np.pi
        i['k'] /= CAL2J * 2
    for i in ff.dihedrals:
        i['k1'] /= CAL2J / 2
        i['k2'] /= CAL2J / 2
        i['k3'] /= CAL2J / 2
        i['k4'] /= CAL2J / 2
    for i in ff.impropers:
        i['k2'] /= CAL2J
    return ff
