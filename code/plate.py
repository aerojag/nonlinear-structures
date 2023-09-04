# Plate definition

import logging

class Plate:

    class Geometry:
        def __init__(self, a: float, b: float, h: float):
            self.a = a
            self.b = b
            self.h = h

        def _set_geometry(self):
            try:
                self.r = self.a / self.b
            except ZeroDivisionError:
                print("Error: division by zero is not allowed, 'b' set equal to 1.0")
                self.r = self.a

    class BoundaryConditions:
        def __init__(self, out_of_plane: str, in_plane: str):
            self.out_of_plane = out_of_plane
            self.in_plane = in_plane

        def _set_boundary_conditions(self, out_of_plane_values, in_plane_values):
            if (all(side_bc in out_of_plane_values for side_bc in self.out_of_plane) and len(self.out_of_plane) == 4) and \
               (self.in_plane.lower() in in_plane_values):
                return self.out_of_plane, self.in_plane.lower()
            else:
                self.out_of_plane = "CCCC"
                self.in_plane = "immovable"
                logging.error("Incorrect boundary conditions, default ones are set: 'CCCC' and 'immovable'.")
                return self.out_of_plane, self.in_plane
    
    def __init__(self):
        self.Geometry = None
        self.BoundaryConditions = None
