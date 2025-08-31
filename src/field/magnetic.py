import numpy as np

class MagneticObject:
    """Base class for any object that generates a magnetic field."""
    def get_field_at(self, point):
        """Returns the magnetic field vector (B) at a given 3D point."""
        # This method must be implemented by all subclasses.
        raise NotImplementedError("Subclasses must implement this method.")
    
class CurrentWire(MagneticObject):
    """
    Represents a current-carrying wire defined by a series of points.
    It can be an open wire or a closed loop.
    """
    # Define MU_0 as a class constant for clarity and efficiency.
    MU_0 = 4e-7 * np.pi  # Permeability of free space

    # NEW: Added 'closed' parameter to the constructor.
    def __init__(self, points, current, closed=False):
        """
        Initializes the wire.

        Args:
            points (list): A list of 3D points (as lists or numpy arrays) defining the wire's path.
            current (float): The current flowing through the wire in Amperes.
            closed (bool): If True, the wire forms a closed loop from the last point to the first.
        """
        self.points = [np.array(p, dtype=float) for p in points]
        self.current = float(current)
        self.closed = closed


class CurrentWire(MagneticObject):
    """
    Represents a current-carrying wire defined by a series of points.
    It can be an open wire or a closed loop.
    """
    # Define MU_0 as a class constant for clarity and efficiency.
    MU_0 = 4e-7 * np.pi  # Permeability of free space
    BS_CONSTANT = MU_0 / (4 * np.pi)  # Biot-Savart constant

    # NEW: Added 'closed' parameter to the constructor.
    def __init__(self, points, current, closed=False):
        """
        Initializes the wire.

        Args:
            points (list): A list of 3D points (as lists or numpy arrays) defining the wire's path.
            current (float): The current flowing through the wire in Amperes.
            closed (bool): If True, the wire forms a closed loop from the last point to the first.
        """
        self.points = [np.array(p, dtype=float) for p in points]
        self.current = float(current)
        self.closed = closed
    
    # Helper function for a single segment's contribution.
    def _biot_savart_segment(self, p1, p2, point):
        """
        Calculates the magnetic field contribution from a single straight wire segment.
        
        This uses a midpoint approximation for the vector r, which is a common
        simplification for numerical simulations.
        """
        dl = p2 - p1  # Vector representing the wire segment
        
        # Avoid calculation for zero-length segments
        #if np.linalg.norm(dl) < 1e-10:
            #return np.zeros(3)

        segment_midpoint = (p1 + p2) / 2
        r = point - segment_midpoint  # Vector from segment to the point
        r_mag = np.linalg.norm(r)

        # Avoid division by zero if the point is on the segment's midpoint
        if r_mag < 1e-10:
            return np.zeros(3)
            
        # Biot-Savart Law: dB = (μ₀ * I) / (4π) * (dL x r) / |r|³
        dB = self.BS_CONSTANT * self.current * np.cross(dl, r) / (r_mag ** 3)
        return dB


    # REFACTORED: The main method is now much cleaner.
    def get_field_at(self, point):
        """Calculates the total magnetic field at a point using the Biot-Savart law."""
        point = np.array(point, dtype=float)
        total_field = np.zeros(3, dtype=float)

        # Iterate through the main wire segments
        for i in range(len(self.points) - 1):
            p1 = self.points[i]
            p2 = self.points[i+1]
            total_field += self._biot_savart_segment(p1, p2, point)

        # If the wire is a closed loop, add the contribution from the last segment to the first
        if self.closed and len(self.points) > 1:
            p1 = self.points[-1] # Last point
            p2 = self.points[0]  # First point
            total_field += self._biot_savart_segment(p1, p2, point)
            
        return total_field
    

class CurrentRing(CurrentWire):
    """Represents a simple, circular ring of current oriented in 3D space."""
    def __init__(self, center, radius, direction, current, num_points=100):
        self.center = np.array(center, dtype=float)
        self.radius = float(radius)
        self.current = float(current)
        self.direction = np.array(direction, dtype=float) / np.linalg.norm(direction)

        # Create a coordinate system for the ring
        # First basis vector (can be any vector perpendicular to direction)
        if not np.allclose(self.direction, [0, 0, 1]):
            v1 = np.cross([0, 0, 1], self.direction)
        else:
            v1 = np.cross([1, 0, 0], self.direction)
        v1 = v1 / np.linalg.norm(v1)
        
        # Second basis vector
        v2 = np.cross(self.direction, v1)

        # Generate points around the circle
        angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
        self.points = [
            self.center + 
            radius * (v1 * np.cos(angle) + v2 * np.sin(angle))
            for angle in angles
        ]
        
        # Call the parent constructor with closed=True
        super().__init__(self.points, current, closed=True)

class ConstantMagneticField(MagneticObject):
    """Represents a constant magnetic field in 3D space."""
    def __init__(self, field_vector):
        self.field_vector = np.array(field_vector, dtype=float)

    def get_field_at(self, point):
        """Returns the constant magnetic field vector at any point."""
        return self.field_vector

class MagneticMoment(MagneticObject):
    """Represents a magnetic moment in 3D space."""
    def __init__(self, position, moment):
        """
        Initializes the MagneticMoment.

        Args:
            position (list or np.ndarray): The 3D position of the magnetic moment.
            moment (list or np.ndarray): The 3D vector of the magnetic moment.
        """
        self.position = np.array(position, dtype=float)
        self.moment = np.array(moment, dtype=float)

    def get_field_at(self, point):
        """
        Returns the magnetic field vector at a given point due to this magnetic moment.

        Args:
            point (list or np.ndarray): The 3D point at which to calculate the field.

        Returns:
            np.ndarray: The magnetic field (B) vector at the given point.
        """
        r_vector = np.array(point, dtype=float) - self.position
        r_magnitude = np.linalg.norm(r_vector)

        if r_magnitude < 5e-5:
            return np.zeros(3)  # Field is singular at the moment's position

        # The magnetic field from a dipole is given by:
        # B = (μ₀ / (4π)) * (3(m·r̂)r̂ - m) / |r|³
        # where r̂ is the unit vector in the direction of r.

        mu_0 = 4 * np.pi * 1e-7  # Permeability of free space (T·m/A)
        r_hat = r_vector / r_magnitude
        m_dot_r_hat = np.dot(self.moment, r_hat)

        # Correct implementation of the dipole field formula
        B = (mu_0 / (4 * np.pi)) * (
            (3 * m_dot_r_hat * r_hat - self.moment) / (r_magnitude**3)
        )

        return B

class MovingPointCharge(MagneticObject):
    """
    Represents a single point charge moving with a constant velocity.
    
    This class calculates the magnetic field generated by the moving charge
    using the Biot-Savart law for a point charge. This is the physically
    correct way to model the magnetic field from a single moving charged particle,
    as opposed to treating it as a magnetic dipole.
    """
    def __init__(self, position, velocity, charge):
        """
        Initializes the MovingPointCharge object.

        Args:
            position (list, tuple, or np.ndarray): The initial 3D position vector 
                                                    [x, y, z] of the charge.
            velocity (list, tuple, or np.ndarray): The 3D velocity vector 
                                                    [vx, vy, vz] of the charge.
            charge (float): The electric charge of the particle (in Coulombs).
        """
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.charge = float(charge)

    def get_field_at(self, point):
        """
        Calculates the magnetic field at a specific point due to this moving charge.

        The calculation is based on the Biot-Savart law for a point charge:
        B(r') = (μ₀ / 4π) * (q * v × R) / |R|³
        where R is the vector from the charge's position to the observation point.

        Args:
            point (list, tuple, or np.ndarray): The 3D coordinates [x, y, z] of the 
                                                observation point.

        Returns:
            np.ndarray: A 3D numpy array for the magnetic field vector B (in Teslas)
                        at the specified point. Returns a zero vector if the point
                        is at the charge's exact position.
        """
        # Define the physical constant, permeability of free space (in T·m/A)
        mu_0 = 4 * np.pi * 1e-7

        # Calculate the vector R from the source (charge's position) to the target (point)
        R_vector = np.array(point, dtype=float) - self.position
        
        # Calculate the magnitude of the R vector
        R_magnitude = np.linalg.norm(R_vector)

        # To avoid division by zero, if the point is at the charge's location,
        # the field is singular. We return a zero vector for this case.
        if R_magnitude < 1e-7:
            return np.zeros(3)

        # Calculate the cross product of the particle's velocity and the R vector
        v_cross_R = np.cross(self.velocity, R_vector)

        # Apply the Biot-Savart formula
        B_field = (mu_0 / (4 * np.pi)) * self.charge * (v_cross_R / (R_magnitude**3))

        return B_field