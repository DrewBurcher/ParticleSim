import numpy as np
from numpy import linalg

class Particle:
    def __init__(self, field, particle_field = None, position=None, velocity=None, bound=True, mass=1.0, 
                 charge=-1.0, classical=True, traced=True, color=None, interaction=True):
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.field = field  # Reference to the electromagnetic field
        self.bound = bound
        self.mass = mass
        self.charge = charge
        self.classical = classical
        self.traced = traced
        self.interaction = interaction
        self.q_m = self.charge / self.mass
        self.particle_field = particle_field
        
        # Set color based on charge if not specified
        if color is None:
            if charge > 0:
                self.color = (1.0, 0.0, 0.0)  # Red for positive
            elif charge < 0:
                self.color = (0.0, 0.7, 1.0)  # Light blue for negative
            else:
                self.color = (1.0, 1.0, 1.0)  # White for neutral
        else:
            self.color = color

    def update(self, dt):
        E = self.field.get_electric_field(self.position)
        B = self.field.get_magnetic_field(self.position)
        if self.interaction:
            E += self.particle_field.get_electric_field(self.position)
            B += self.particle_field.get_magnetic_field(self.position)
        self.update_boris(E, B, dt)
        # Handle boundaries if enabled
        if self.bound:
            self.bounce_boundaries()
        
    def update_boris(self, E, B, dt):
        # Half acceleration from electric field
        v_minus = self.velocity + self.q_m* E * dt/2

        # Magnetic field rotation
        t = self.q_m* B * dt/2
        t_squared = np.dot(t, t)

        
        s = 2*t/(1 + t_squared)
        v_prime = v_minus + np.cross(v_minus, t)
        v_plus = v_minus + np.cross(v_prime, s)
        #else:
            #v_plus = v_minus

        # Half acceleration from electric field
        self.velocity = v_plus + self.q_m * E * dt/2
        
        # Update position
        self.position += self.velocity * dt



    def bounce_boundaries(self, boundary_size=10.0):
        for i in range(3):
            if abs(self.position[i]) > boundary_size:
                self.position[i] = np.sign(self.position[i]) * boundary_size
                self.velocity[i] *= -1

    def add_particle_field(self, particle_field):
        """Add a particle field to this particle."""
        self.particle_field = particle_field