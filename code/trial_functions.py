# Trial functions based on Legendre polynomials

import numpy as np
import matplotlib.pyplot as plt

from plate import Plate

class Legendre():
    """
    A class for computing Legendre Polynomials.
    """
    def __init__(self, x, max_order: int):
        """
        Initialise Legendre class.
        
        Args:
            x (numpy.ndarray): Generic coordinate.
            max_order (int): Highest order of Legendre Polynomials to compute.
        """
        self.x = x
        self.max_order = max_order
        self.P = np.ones((max_order + 1, len(x))) # Initialise array storing Legendre polynomial
        self.dP = np.zeros((max_order + 1, len(x))) # Initialise array storing its first derivatives
        self.d2P = np.zeros((max_order + 1, len(x))) # Initialise array storing its second derivatives

    def get_legendre(self, n: int):
        """
        Compute Legendre polynomials and their first and second derivatives
        using Bonnet's recursion formula.
        
        Args:
            n (int): Order of Legendre polynomial.
        """
        if n == 0:
            self.P[n + 1] = self.x
        else:
            self.P[n + 1] = ((2 * n + 1) * self.x * self.P[n] - n * self.P[n - 1]) / (n + 1)
        
        self.dP[n + 1] = (n + 1) * self.P[n] + self.x * self.dP[n]
        self.d2P[n + 1] = (n + 2) * self.dP[n] + self.x * self.d2P[n]

        return self

    def compute_legendre(self):
        """
        Update Legendre polynomials' order and their derivatives.
        """
        for order in range(self.max_order):
            self.get_legendre(order)
        
        return self
    
class TrialFunctions(Legendre):
    """
    A class for computing trial functions based on Legendre Polynomials.
    """
    def __init__(self, x, max_order):
        super().__init__(x, max_order)

    TRIAL_FUNCTION_BASE_DOCSTRING = """
    Calculate the trial function for {function_type}.
    
    Args:
        i_index (int): BC index (0 for 'F', 1 for 'S' and 'C').
        j_index (int): BC index (0 for 'F', 1 for 'S' and 'C').
        k_index (int, optional): BC index (0 for 'F' and 'S', 1 for 'C'); only used in trial_phi() and trial_dphi().
        l_index (int, optional): BC index (0 for 'F' and 'S', 1 for 'C'); only used in trial_phi() and trial_dphi().
        
    Returns:
        float: Value of the trial function, its first or its second derivative.
    """
        
    # Displacements (X, Y)
    def get_x(self, i_index, j_index):
        self.TRIAL_FUNCTION_BASE_DOCSTRING.format(function_type="displacements")
        return (1 - self.x) ** i_index * (1 + self.x) ** j_index * self.P[self.max_order]
    
    def get_dx(self, i_index, j_index):
        self.TRIAL_FUNCTION_BASE_DOCSTRING.format(function_type="displacements")
        first_term  = - (1 + self.x) ** j_index * self.P[self.max_order] if i_index else 0
        second_term =   (1 - self.x) ** i_index * self.P[self.max_order] if j_index else 0
        third_term  =   (1 - self.x) ** i_index * (1 + self.x) ** j_index * self.dP[self.max_order]

        return first_term + second_term + third_term

    def get_d2x(self, i_index, j_index):
        self.TRIAL_FUNCTION_BASE_DOCSTRING.format(function_type="displacements")
        if i_index and j_index:
            return - 2 * self.P[self.max_order] - 4 * self.x * self.dP[self.max_order] + \
                    (1 - self.x) * (1 + self.x) * self.d2P[self.max_order]
        elif i_index:
            return - 2 * self.dP[self.max_order] + (1 - self.x) * self.d2P[self.max_order]
        elif j_index:
            return 2 * self.dP[self.max_order] + (1 + self.x) * self.d2P[self.max_order]
        else:
            return self.d2P[self.max_order]
    
    # Section rotations (PHI, PSI)
    def get_phi(self, i_index, j_index, k_index=None, l_index=None):
        self.TRIAL_FUNCTION_BASE_DOCSTRING.format(function_type="section rotations")
        if k_index is None and l_index is None:
            return self.get_x(i_index, j_index)
        else:
            return (1 - self.x) ** (i_index * k_index) * (1 + self.x) ** (j_index * l_index) * self.P[self.max_order]
        
    def get_dphi(self, i_index, j_index, k_index=None, l_index=None):
        self.TRIAL_FUNCTION_BASE_DOCSTRING.format(function_type="section rotations")
        if (k_index is None and l_index is None) or all((i_index, j_index, k_index, l_index)):
            return self.get_dx(i_index, j_index)
        elif i_index and k_index:
            return - self.P[self.max_order] + (1 - self.x) * self.dP[self.max_order]
        elif j_index and l_index:
            return self.P[self.max_order] + (1 + self.x) * self.dP[self.max_order]
        else:
            return self.dP[self.max_order]
    
    # Airy function (F, G, nu, psi)
    def get_airy(self):
        return (1 - self.x ** 2) ** 2 * self.P[self.max_order]
    
    def get_dairy(self):
        return - 4 * self.x * (1 - self.x ** 2) * self.P[self.max_order] + \
                (1 - self.x ** 2) ** 2 * self.dP[self.max_order]
        
    def get_d2airy(self):
        return   4 * (3 * self.x ** 2 - 1) * self.P[self.max_order] \
               - 8 * self.x * (1 - self.x ** 2) * self.dP[self.max_order] + \
                (1 - self.x ** 2) ** 2 * self.d2P[self.max_order]
    
    def get_nu(self, bc_in_plane):
        if 'immovable' in bc_in_plane:
            return self.P[self.max_order]
        elif 'movable' in bc_in_plane:
            return [self.P[i] for i in range(1, self.max_order)]
        else:
            return None

# TESTING
if __name__ == '__main__':

    # Handle boundary conditions
    """
    Boundary conditions, order of edges below:

    _____4_____
    |         |
    |   y^    |
    1|    |-> x|3
    |         |
    |_________|
        2

    Out-of-plane boundary conditions:
    - Free: 'F'.
    - Simply supported: 'S'.
    - Clamped: 'C'.

    In-plane boundary conditions:
    - Completely free.
    - Movable.
    - Immovable.
    """

    bc = Plate.BoundaryConditions('CCCC', 'movable')
    bc_op = bc.out_of_plane
    bc_ip = bc.in_plane

    bc_op_cases = {'F': (0, 0), 'S': (1, 0), 'C': (1, 1)}
    i_bc = []
    k_bc = []
    for bc in bc_op:
        i, k = bc_op_cases.get(bc)
        i_bc.append(i)
        k_bc.append(k)

    i1, i2, j1, j2 = i_bc
    k1, k2, l1, l2 = k_bc

    # Create trial function instances and display
    L = 1000
    N = 5
    x_values = np.linspace(-1.0, 1.0, num = L)
    plt.figure(figsize=(20, 12))
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.size"] = 25.0
    for i in range(N + 1):
        trial = TrialFunctions(x_values, i)
        legendre_0 = trial.compute_legendre().P[i]
        trial_x = trial.get_x(i1, j1)
        trial_nu = trial.get_nu(bc_ip)
        # print(legendre_0)
        # plt.plot(x_values, legendre_0, label=fr'$P_{i}$', linewidth=3.0)
        print(trial_nu)
        plt.plot(x_values, trial_x, label=fr'$X_{i}$', linewidth=3.0)

    plt.xlabel(r'$x$', fontsize=25.0)
    plt.ylabel(r'$P_i(x)$', fontsize=25.0)
    plt.legend(loc='lower right', fontsize=25.0)
    plt.grid()
    plt.savefig('trial_functions.png')
