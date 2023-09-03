import sys
import matplotlib.pyplot as plt

sys.path.append("/hpc/ac48390/git/nonlinear_vibrations/code")
from plate import *
from trial_functions import *

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
    print(trial_x)
    plt.plot(x_values, trial_x, label=fr'$X_{i}$', linewidth=3.0)
    # print(trial_nu)
    # plt.plot(x_values, trial_nu, label=fr'$X_{i}$', linewidth=3.0) if trial_nu is not None else None

plt.xlabel(r'$x$', fontsize=25.0)
plt.ylabel(r'$P_i(x)$', fontsize=25.0)
plt.legend(loc='lower right', fontsize=25.0)
plt.grid()
plt.savefig('../results/trial_functions.png')
