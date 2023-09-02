# Trial functions based on Legendre polynomials

import numpy as np

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

    def _get_legendre(self, n: int):
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
            self._get_legendre(order)
        
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

    AIRY_FUNCTION_BASE_DOCSTRING = """
    Calculate the Airy stress trial function.
    
    Args:
        bc_in_plane (str): In-plane BC; only used in trial_nu().

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
        self.AIRY_FUNCTION_BASE_DOCSTRING
        return (1 - self.x ** 2) ** 2 * self.P[self.max_order]
    
    def get_dairy(self):
        self.AIRY_FUNCTION_BASE_DOCSTRING
        return - 4 * self.x * (1 - self.x ** 2) * self.P[self.max_order] + \
                (1 - self.x ** 2) ** 2 * self.dP[self.max_order]
        
    def get_d2airy(self):
        self.AIRY_FUNCTION_BASE_DOCSTRING
        return   4 * (3 * self.x ** 2 - 1) * self.P[self.max_order] \
               - 8 * self.x * (1 - self.x ** 2) * self.dP[self.max_order] + \
                (1 - self.x ** 2) ** 2 * self.d2P[self.max_order]
    
    def get_nu(self, bc_in_plane):
        self.AIRY_FUNCTION_BASE_DOCSTRING
        if 'immovable' in bc_in_plane:
            return self.P[self.max_order]
        elif 'movable' in bc_in_plane:
            return self.P[self.max_order] if self.max_order > 0 else None
        else:
            None
        