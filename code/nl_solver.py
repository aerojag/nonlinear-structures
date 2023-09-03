# Nonlinear solver

import os
import json

from plate import *
from trial_functions import *

module_dir = os.path.dirname(__file__)
print(module_dir)
res_path = os.path.join(module_dir, '../resources/')
with open(res_path + 'data.json', 'r') as data_file:
    data = json.load(data_file)

a = 1
b = 1
h = 0.001
plate1 = Plate()
plate1.Geometry(a, b, h)
plate1.BoundaryConditions('CCCC', 'immovable')

bc_op = plate1.out_of_plane
bc_ip = plate1.in_plane

bc_op_cases = {'F': (0, 0), 'S': (1, 0), 'C': (1, 1)}
i_bc = []
k_bc = []
for bc in bc_op:
    i, k = bc_op_cases.get(bc)
    i_bc.append(i)
    k_bc.append(k)

i1, i2, j1, j2 = i_bc
k1, k2, l1, l2 = k_bc

class NonlinearSolver():

    def __init__(self):
        pass