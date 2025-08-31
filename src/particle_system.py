from collections import deque
import numpy as np
from .particle import Particle
from .field.emfield import ElectroMagneticField as eElectroMagneticField
from .field.electric import ElectricObject
from .field.magnetic import MagneticObject
from .field.electric import ElectricPoint
from .field.magnetic import MagneticMoment
from .field.magnetic import MovingPointCharge

class ParticleSystem:
    def __init__(self, field, particles, trail_length=25000, interaction = True, charge_scale=1.0):
        self.field = field
        self.particles = particles
        self.num_particles = len(particles)
        self.interaction = interaction
        self.charge_scale = charge_scale
        if(interaction):
            self.particle_field = eElectroMagneticField()
            for particle in self.particles:
                particle.add_particle_field(self.particle_field)
        
        self.trails = [deque(maxlen=trail_length) if p.traced else None 
                      for p in self.particles]
        self.update_arrays()
    
    def update_arrays(self):
        self.positions = np.array([p.position for p in self.particles])
        self.velocities = np.array([p.velocity for p in self.particles])
        
        for i, (particle, trail) in enumerate(zip(self.particles, self.trails)):
            if particle.traced and trail is not None:
                trail.append(particle.position.copy())

    def update(self, dt):
        if(self.interaction):
            self.update_particle_field()
        """Update particles using CPU."""
        for particle in self.particles:
            particle.update(dt)  # Only pass dt, particle gets E/B itself
        self.update_arrays()

    def get_particle_electric_field(self, index):
        """Get the electric field at the particle's position."""
        if 0 <= index < self.num_particles:
            return self.field.get_electric_field(self.particles[index].position)
        else:
            raise IndexError("Particle index out of range.")
        
    def update_particle_field(self):
        """Update the particle field with the current state of particles."""
        self.particle_field.electric_objects.clear()
        self.particle_field.magnetic_objects.clear()
        
        for particle in self.particles:
            if particle.interaction:
                # Add a point charge for each traced particle
                self.particle_field.add_electric_object(
                    ElectricPoint(
                        position=particle.position,
                        charge=particle.charge * self.charge_scale
                    )
                )
                self.particle_field.add_magnetic_object(
                    MovingPointCharge(
                        position=particle.position,
                        velocity=particle.velocity,
                        charge=particle.charge * self.charge_scale
                    )
                )
        
        # Update magnetic field if needed
        # This can be extended to include magnetic effects from particles
        # For now, we assume no magnetic contributions from particles