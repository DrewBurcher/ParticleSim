import sys
import numpy as np
from PyQt5.QtWidgets import QApplication
from src.field.electric import ConstantElectricField, ElectricPlate
from src.field.magnetic import ConstantMagneticField, CurrentRing
from src.gui.main_window import MainWindow
import src.field

def setup_simulation():
    # Physical constants
    ELECTRON_MASS = 9.1093837015e-31  # kg
    ELECTRON_CHARGE = -1.602176634e-19  # C

    # Initial conditions for particle
    initial_position = np.array([0, 0.0, 0])
    initial_velocity = np.array([0, 0, 0])  # 1000 km/s in x-y plane

    N = 10  # Number of positive and negative particles
    np.random.seed(103)  # For reproducibility

    particles = []
    for i in range(N):
        # Positive particle
        particles.append({
            'mass': ELECTRON_MASS,
            'charge': ELECTRON_CHARGE,
            'bound': False,
            'classical': True,
            'position': (4e-2*(np.array([(np.random.rand()-0.5)*4,(np.random.rand()-0.5)*4,0])+np.random.uniform(-1, 1, 3))).tolist(),
            'velocity': (np.random.uniform(-1e6, 1e6, 3)).tolist(),
            'traced': True
        })
        # Negative particle
        particles.append({
            'mass': ELECTRON_MASS,
            'charge': -ELECTRON_CHARGE,
            'bound': False,
            'classical': True,
            'position': (4e-2*(np.array([(np.random.rand()-0.5)*4,(np.random.rand()-0.5)*4,0])+np.random.uniform(-1, 1, 3))).tolist(),
            'velocity': (np.random.uniform(-1e6, 1e6, 3)),#(np.random.uniform(-1e4, 1e4, 3)).tolist(),
            'traced': True
        })

    config = {
        'particle': particles,
        'display': {
            'trail_length': 1000,
            'point_size': 5.0,
            'line_width': 2.0,
            'window_size': (1600, 1200),
            'show_time': True,
            'zoom': 1e2,
            'sim_render_ratio': 10
        },
        'physics': {
            'time_step': 1e-9,  # 1 picosecond
            'boundary_size': 1.0,
            'charge_scale': 8e6  # Scale for particle charge in field calculations
        }
    }

    field_sources = {
        'magnetic': [
            CurrentRing(
                center=np.array([0,0,0]),
                current=500,  # -1 kA
                radius=5e-2,  # 1 cm radius
                direction = np.array([0,0,1]),  # Current in z direction
                num_points=20  # Number of points in the ring
            )
        ],
        'electric': [
            #ConstantElectricField([1e3, 0, 0]),  # Uniform field in z direction
            """
            #ConstantElectricField([0, 0, 0]),  # Uniform field in z direction
            ElectricPlate(
                position=[-1e-2, 0, 0.0],
                voltage=1e3,
                Dx=1e-1,  # Plate width
                Dy=1e-2,  # Plate height
                Nx=10,  # Number of points in x direction
                Ny=10,  # Number of points in y direction
                normal=[1, 0, 0],  # Normal vector pointing up
                rotation=0  # No rotation
            ),
            ElectricPlate(
                position=[1e-2, 0, 0.0],
                voltage=1e3, 
                Dx=1e-1,  # Plate width
                Dy=1e-2,  # Plate height
                Nx=10,  # Number of points in x direction
                Ny=10,  # Number of points in y direction
                normal=[1, 0, 0],  # Normal vector pointing up
                rotation=0  # No rotation
            ),
            ElectricPlate(
                position=[0, 1e-2, 0.0],
                voltage=-1e3,  # 1 MV
                Dx=1e-2,  # Plate width
                Dy=1e-1,  # Plate height
                Nx=10,  # Number of points in x direction
                Ny=10,  # Number of points in y direction
                normal=[0, 1, 0],  # Normal vector pointing up
                rotation=0  # No rotation
             ),
            ElectricPlate(
                position=[0, -1e-2, 0.0],
                voltage=1e3,  # 1 MV
                Dx=1e-2,  # Plate width
                Dy=1e-1,  # Plate height
                Nx=10,  # Number of points in x direction
                Ny=10,  # Number of points in y direction
                normal=[0, 1, 0],  # Normal vector pointing up
                rotation=0  # No rotation
             ),
             """
        ]
    }

    return config, field_sources

if __name__ == '__main__':
    # Get configuration and field sources
    config, field_sources = setup_simulation()
    
    # Create application
    app = QApplication(sys.argv)
    
    # Create main window with configuration
    window = MainWindow(config=config, field_sources=field_sources)
    window.show()
    
    sys.exit(app.exec_())
