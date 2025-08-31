import numpy as np
from .particle_system import ParticleSystem
from .field.emfield import ElectroMagneticField
from .particle import Particle

class Universe:
    """Manages the physics simulation including particles and fields."""
    
    def __init__(self, config=None):
        try:
            self.config = config or {}
            self.running = True
            self.field = ElectroMagneticField()
            self.time = 0.0
            
            # Get physics configuration
            physics_config = self.config.get('physics', {})
            self.time_step = physics_config.get('time_step', 1e-12)  # Default 1 picosecond
            self.charge_scale = physics_config.get('charge_scale', 1.0)
            # Create particle system
            particle_list = self.config.get('particle', [])
            display_config = self.config.get('display', {})
            trail_length = display_config.get('trail_length', 25000)

            # Create Particle objects from the list of dicts
            particles = [
                Particle(
                    field=self.field,
                    position=p.get('position'),
                    velocity=p.get('velocity'),
                    mass=p.get('mass', 9.1093837015e-31),
                    charge=p.get('charge', -1.602176634e-19),
                    bound=p.get('bound', True),
                    classical=p.get('classical', True),
                    traced=p.get('traced', True),
                    color=p.get('color', None)
                )
                for p in particle_list
            ]

            self.particle_system = ParticleSystem(
                field=self.field,
                particles=particles,
                trail_length=trail_length,
                charge_scale=self.charge_scale
            )
            
        except Exception as e:
            print(f"Error initializing Universe: {e}")
            raise

    def add_field_source(self, source, field_type='both'):
        """Add electromagnetic field sources to the universe."""
        if field_type in ['electric', 'both']:
            self.field.add_electric_object(source)
        if field_type in ['magnetic', 'both']:
            self.field.add_magnetic_object(source)
    
    def update(self):
        """Update the state of the universe for one time step."""
        if not self.running:
            return
        self.time += self.time_step
        self.particle_system.update(self.time_step)

    def toggle_simulation(self):
        """Pause/resume the simulation."""
        self.running = not self.running
    
    def set_time_step(self, dt):
        """Change the simulation time step."""
        self.time_step = dt
    
    def reset(self):
        """Reset the simulation to initial conditions."""
        self.time = 0.0
        # Re-create particles from config
        particle_list = self.config.get('particle', [])
        display_config = self.config.get('display', {})
        trail_length = display_config.get('trail_length', 25000)
        particles = [
            Particle(
                field=self.field,
                position=p.get('position'),
                velocity=p.get('velocity'),
                mass=p.get('mass', 9.1093837015e-31),
                charge=p.get('charge', -1.602176634e-19),
                bound=p.get('bound', True),
                classical=p.get('classical', True),
                traced=p.get('traced', True),
                color=p.get('color', None)
            )
            for p in particle_list
        ]
        self.particle_system = ParticleSystem(
            field=self.field,
            particles=particles,
            trail_length=trail_length
        )
        self.running = True