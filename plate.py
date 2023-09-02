# Plate definition

class Plate:

    class Geometry:
        def __init__(self, a: float, b: float, h: float):
            self.a = a
            self.b = b
            self.h = h
            try:
                self.r = a / b
            except ZeroDivisionError:
                print("Error: division by zero is not allowed, b set equal to 1.0")
                self.r = a

    class BoundaryConditions:
        def __init__(self, out_of_plane: str, in_plane: str):
            try:
                if isinstance(out_of_plane, str):
                    self.out_of_plane = out_of_plane
                else:
                    raise TypeError("Attribute 'out_of_plane' must be a string.")

                if isinstance(in_plane, str):
                    self.in_plane = in_plane
                else:
                    raise TypeError("Attribute 'in_plane' must be a string.")
            except TypeError:
                print("Error: boundary conditions must be defined as strings, \
                set default boundary conditions 'SSSS' & 'immovable'.")
                self.out_of_plane = 'SSSS'
                self.in_plane = 'immovable'
    
    def __init__(self, a, b, h, out_of_plane, in_plane):
        self.Geometry = self.Geometry(a, b, h)
        self.BoundaryConditions = self.BoundaryConditions(out_of_plane, in_plane)
