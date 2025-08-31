import numpy as np

class ElectricObject:
    """Base class for any object that generates an electric field."""
    def get_field_at(self, point):
        """Returns the electric field vector (E) at a given 3D point."""
        # This method must be implemented by all subclasses.
        raise NotImplementedError("Subclasses must implement this method.")
    
class ElectricSurface(ElectricObject):
    def __init__(self, position, voltage):
        self.position = np.array(position, dtype=float)
        self.voltage = float(voltage)
        self.points = []

    def get_field_at(self, point):
        total_field = np.zeros(3, dtype=float)
        for electric_point in self.points:
            total_field += electric_point.get_field_at(point)
        return total_field

class ElectricPoint(ElectricObject):
    """Represents a point charge in 3D space."""
    EPSILON_0 = 8.854187817e-12  # Electric constant (F/m)
    COULOMB_CONSTANT = 1 / (4 * np.pi * EPSILON_0)  # Coulomb's constant

    def __init__(self, position, charge):
        self.position = np.array(position, dtype=float)
        self.charge = float(charge)

    def get_field_at(self, point):
        """Returns the electric field vector at a given point due to this point charge."""
        r_vector = np.array(point, dtype=float) - self.position
        r_magnitude = np.linalg.norm(r_vector)
        if r_magnitude < 5e-5:
            return np.zeros(3)  # Field is undefined at the charge's position
        return self.COULOMB_CONSTANT * self.charge * (r_vector / r_magnitude**3)

class ElectricPlate(ElectricSurface):
    """Represents a charged plate in 3D space."""
    def __init__(self, position, voltage, Dx, Dy, Nx=10, Ny=10, normal=[0,0,0], rotation=0):
        super().__init__(position, voltage)
        self.Dx = Dx  # Store dimensions
        self.Dy = Dy
        self.points = []  # Initialize as empty list first
        self.normal = np.array(normal, dtype=float)
        if np.any(self.normal):
            self.normal = self.normal / np.linalg.norm(self.normal)
        self.rotation = rotation  # Store rotation angle
        dA = (Dx * Dy) / (Nx * Ny)  # Area per element
        
        # Constants
        EPSILON_0 = 8.854187817e-12  # Electric constant (F/m)
        charge_density = EPSILON_0 * self.voltage  # Simplified parallel plate approximation
        element_charge = charge_density * dA  # Total charge for this element

        # Create grid of points
        for i in range(Nx):
            for j in range(Ny):
                x = (i + 0.5) * Dx/Nx - Dx/2
                y = (j + 0.5) * Dy/Ny - Dy/2
                local_position = np.array([x, y, 0.0], dtype=float)
                
                # Apply rotation if specified
                if rotation != 0:
                    rot_z = np.array([
                        [np.cos(rotation), -np.sin(rotation), 0],
                        [np.sin(rotation), np.cos(rotation), 0],
                        [0, 0, 1]
                    ])
                    local_position = rot_z @ local_position

                # Align with normal vector
                if np.any(normal):
                    z_axis = np.array([0, 0, 1])
                    normal = normal / np.linalg.norm(normal)
                    v = np.cross(z_axis, normal)
                    s = np.linalg.norm(v)
                    c = np.dot(z_axis, normal)
                    if s != 0:
                        v_cross = np.array([
                            [0, -v[2], v[1]],
                            [v[2], 0, -v[0]],
                            [-v[1], v[0], 0]
                        ])
                        rot_normal = np.eye(3) + v_cross + v_cross @ v_cross * (1 - c) / (s * s)
                        local_position = rot_normal @ local_position

                # Add point charge at the calculated position
                point_position = self.position + local_position
                self.points.append(ElectricPoint(point_position, element_charge))
    def get_field_at(self, point):
        """Returns the electric field vector at a given point due to this plate."""
        return super().get_field_at(point)
    

class ConstantElectricField(ElectricObject):
    """Represents a constant electric field in 3D space."""
    def __init__(self, field_vector):
        self.field_vector = np.array(field_vector, dtype=float)

    def get_field_at(self, point):
        """Returns the constant electric field vector at any point."""
        return self.field_vector


