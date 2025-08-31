import numpy as np
from . import electric
from . import magnetic







class ElectroMagneticField:
    """A container to hold all electric and magnetic objects and calculate total fields."""
    def __init__(self):
        self.electric_objects = []
        self.magnetic_objects = []

    def add_electric_object(self, obj):
        if isinstance(obj, electric.ElectricObject):
            self.electric_objects.append(obj)

    def add_magnetic_object(self, obj):
        if isinstance(obj, magnetic.MagneticObject):
            self.magnetic_objects.append(obj)
            
    def get_electric_field(self, point):
        """Calculates the total E-field at a point by superposition."""
        total_e_field = np.array([0.0, 0.0, 0.0])
        for obj in self.electric_objects:
            total_e_field += obj.get_field_at(point)
        return total_e_field
    

    def get_magnetic_field(self, point):
        """Calculates the total B-field at a point by superposition."""
        total_b_field = np.array([0.0, 0.0, 0.0])
        for obj in self.magnetic_objects:
            total_b_field += obj.get_field_at(point)
        return total_b_field