from PyQt5.QtOpenGL import QGLWidget, QGLFormat
from PyQt5.QtCore import Qt, QTimer
from OpenGL.GL import *
from OpenGL.GLU import *
import numpy as np
import time

from src.field.electric import ElectricPlate
from src.field.magnetic import CurrentRing
from ..universe import Universe

class GLWidget(QGLWidget):
    def __init__(self, parent=None, config=None, field_sources=None):
        try:
            # Create GL format first
            fmt = QGLFormat()
            fmt.setDepthBufferSize(24)
            fmt.setDoubleBuffer(True)
            fmt.setSampleBuffers(True)
            fmt.setSwapInterval(0)  # Disable VSync
            
            # Pass format to parent constructor
            super().__init__(fmt, parent)
            
            # Store config and get display settings
            self.config = config or {}
            display_config = self.config.get('display', {})
            self.scale = display_config.get('zoom', 1.0)
            self.sim_render_ratio = display_config.get('sim_render_ratio', 10)  # Add this line
            
            # Set window size from config
            window_size = display_config.get('window_size', (800, 600))
            self.setMinimumSize(*window_size)
            
            # Initialize other properties
            self.rotation = [30, 30, 0]  # Start with a better view angle
            self.last_pos = None
            self.universe = None
            self.field_sources = field_sources
            self.show_individual_charges = False
            
            # Initialize FPS tracking
            self.fps = 0
            self.frame_count = 0
            self.last_fps_time = time.time()
            
            # Use QTimer for simulation updates
            self.sim_timer = QTimer(self)
            self.sim_timer.timeout.connect(self.update_simulation)
            self.sim_timer.start(0)  # Run as fast as possible
            
            # Separate timer for display updates
            self.display_timer = QTimer(self)
            self.display_timer.timeout.connect(self.updateGL)
            self.display_timer.start(16)  # ~60 FPS for display
            
        except Exception as e:
            print(f"Error in GLWidget initialization: {e}")
            raise

    def initializeGL(self):
        try:
            # Create universe here when GL context is ready
            self.universe = Universe(self.config)
            
            # Add field sources if provided
            if self.field_sources:
                for source_type, sources in self.field_sources.items():
                    for source in sources:
                        self.universe.add_field_source(source, source_type)

            # Set up OpenGL state
            glEnable(GL_DEPTH_TEST)
            glEnable(GL_POINT_SMOOTH)
            glEnable(GL_LINE_SMOOTH)
            glEnable(GL_BLEND)
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
            glClearColor(0.0, 0.0, 0.0, 1.0)
            
            # Set display parameters
            display_config = self.config.get('display', {})
            glPointSize(display_config.get('point_size', 5.0))
            glLineWidth(display_config.get('line_width', 2.0))
            
        except Exception as e:
            print(f"OpenGL initialization error: {e}")
            raise

    def update_simulation(self):
        """Update simulation state multiple times per render"""
        if self.universe and self.universe.running:
            # Run multiple simulation steps per render
            for _ in range(self.sim_render_ratio):
                self.universe.update()
        
            # Update FPS counter
            self.frame_count += 1
            current_time = time.time()
            if current_time - self.last_fps_time > 1.0:
                self.fps = self.frame_count
                self.frame_count = 0
                self.last_fps_time = current_time

    def draw_axes(self):
        """Draws X, Y, and Z axes at fixed scale."""
        axis_length = 10.0  # Fixed length for reference
        neg_axis_brightness = 0.25  # Dimmed axes for better visibility
        
        glBegin(GL_LINES)
        # X-axis in Red
        glColor3f(1.0, 0.0, 0.0)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(axis_length, 0.0, 0.0)
        glColor3f(1.0 * neg_axis_brightness, 0.0, 0.0)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(-axis_length, 0.0, 0.0)
        
        # Y-axis in Green
        glColor3f(0.0, 1.0, 0.0)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(0.0, axis_length, 0.0)
        glColor3f(0.0, 1.0 * neg_axis_brightness, 0.0)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(0.0, -axis_length, 0.0)
        
        # Z-axis in Blue
        glColor3f(0.0, 0.0, 1.0)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(0.0, 0.0, axis_length)
        glColor3f(0.0, 0.0, 1.0 * neg_axis_brightness)
        glVertex3f(0.0, 0.0, 0.0)
        glVertex3f(0.0, 0.0, -axis_length)
        glEnd()

    def draw_tracer(self):
        glColor3f(0.5, 0.5, 0.5)
        glBegin(GL_LINE_STRIP)
        for pos in self.universe.particle_system.trail:
            glVertex3fv(pos)
        glEnd()

    def draw_field_sources(self):
        """Draws all electromagnetic field sources."""
        # Draw magnetic sources
        for source in self.universe.field.magnetic_objects:
            if isinstance(source, CurrentRing):
                glColor3f(1.0, 1.0, 0.0)  # Yellow
                glBegin(GL_LINE_LOOP)
                for point in source.points:
                    scaled_point = point * self.scale
                    glVertex3fv(scaled_point)
                glEnd()

        # Draw electric sources
        for source in self.universe.field.electric_objects:
            if isinstance(source, ElectricPlate):
                self.draw_electric_plate(source)

    def paintGL(self):
        if not self.universe:
            return
        
        try:
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            glLoadIdentity()

            glTranslatef(0.0, 0.0, -40.0)
            glRotatef(self.rotation[0], 1.0, 0.0, 0.0)
            glRotatef(self.rotation[1], 0.0, 1.0, 0.0)
            
            self.draw_axes()
            self.draw_field_sources()
            
            # Draw particles with their colors
            glBegin(GL_POINTS)
            for i, pos in enumerate(self.universe.particle_system.positions):
                particle = self.universe.particle_system.particles[i]
                glColor3fv(particle.color)
                scaled_pos = pos * self.scale
                glVertex3fv(scaled_pos)
            glEnd()

            # Draw trails for traced particles
            for i, trail in enumerate(self.universe.particle_system.trails):
                if trail and len(trail) > 1:  # Check if trail exists and has points
                    particle = self.universe.particle_system.particles[i]
                    glColor3f(*[c * 0.5 for c in particle.color])  # Dimmer trail color
                    glBegin(GL_LINE_STRIP)
                    for pos in trail:
                        scaled_pos = pos * self.scale
                        glVertex3fv(scaled_pos)
                    glEnd()

        except Exception as e:
            print(f"Error in paintGL: {e}")
            self.universe = None

    def resizeGL(self, w, h):
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        # Handle case where h is zero to avoid division by zero
        if h == 0:
            h = 1
        gluPerspective(45.0, float(w)/float(h), 0.1, 100.0)
        glMatrixMode(GL_MODELVIEW)

    def mousePressEvent(self, event):
        """Handle mouse press events."""
        self.last_pos = event.pos()

    def mouseMoveEvent(self, event):
        """Handle mouse movement for rotation."""
        if event.buttons() & Qt.LeftButton and self.last_pos is not None:
            dx = event.x() - self.last_pos.x()
            dy = event.y() - self.last_pos.y()
            
            # Update rotation with sensitivity adjustment
            self.rotation[0] = (self.rotation[0] + dy/2) % 360
            self.rotation[1] = (self.rotation[1] + dx/2) % 360
            
            self.last_pos = event.pos()
            self.updateGL()

    def draw_current_ring(self, ring):
        """Draws a current ring in yellow."""
        glColor3f(1.0, 1.0, 0.0)  # Yellow
        glBegin(GL_LINE_LOOP)
        for point in ring.points:
            glVertex3fv(point)
        glEnd()

    def draw_electric_plate(self, plate):
        """Draws an electric plate with color based on voltage and optional individual charges."""
        # Set color based on voltage
        if plate.voltage > 0:
            plate_color = (1.0, 0.2, 0.2, 0.3)  # Red for positive
        else:
            plate_color = (0.2, 0.7, 1.0, 0.3)  # Light blue for negative

        # Apply rotation if specified
        if plate.rotation != 0:
            # Convert rotation angle to radians
            angle_rad = np.radians(plate.rotation)
            rot_z = np.array([
                [np.cos(angle_rad), -np.sin(angle_rad), 0],
                [np.sin(angle_rad), np.cos(angle_rad), 0],
                [0, 0, 1]
            ])
            # Apply rotation to up vector before calculating corners
            up = rot_z @ np.array([0.0, 1.0, 0.0])
        else:
            up = np.array([0.0, 1.0, 0.0])

        # Calculate plate corners with rotation
        right = np.cross(plate.normal, up)
        if np.allclose(right, 0):
            right = np.array([1.0, 0.0, 0.0])
        right = right / np.linalg.norm(right)
        up = np.cross(right, plate.normal)  # Recalculate up to ensure orthogonality
        
        # Calculate corners using plate's actual dimensions
        corners = [
            plate.position + right * plate.Dx/2 + up * plate.Dy/2,
            plate.position + right * plate.Dx/2 - up * plate.Dy/2,
            plate.position - right * plate.Dx/2 - up * plate.Dy/2,
            plate.position - right * plate.Dx/2 + up * plate.Dy/2
        ]

        # Draw plate surface with proper scaling and blending
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        
        # Draw plate surface
        glColor4f(*plate_color)
        glBegin(GL_QUADS)
        for corner in corners:
            scaled_corner = corner * self.scale
            glVertex3fv(scaled_corner)
        glEnd()

        # Draw plate border
        glColor4f(*plate_color[:3], 1.0)
        glBegin(GL_LINE_LOOP)
        for corner in corners:
            scaled_corner = corner * self.scale
            glVertex3fv(scaled_corner)
        glEnd()
        
        # Draw normal vector
        glColor3f(*plate_color[:3])
        glBegin(GL_LINES)
        normal_length = min(plate.Dx, plate.Dy) * 0.2
        base_pos = plate.position * self.scale
        normal_end = (plate.position + plate.normal * normal_length) * self.scale
        glVertex3fv(base_pos)
        glVertex3fv(normal_end)
        glEnd()
        
        # Draw individual charges if enabled
        if self.show_individual_charges and hasattr(plate, 'points'):
            glPointSize(3.0)
            glBegin(GL_POINTS)
            for point in plate.points:
                scaled_pos = point.position * self.scale
                if point.charge > 0:
                    glColor3f(1.0, 0.0, 0.0)
                else:
                    glColor3f(0.0, 0.7, 1.0)
                glVertex3fv(scaled_pos)
            glEnd()
    
        glDisable(GL_BLEND)

    def keyPressEvent(self, event):
        """Handle key press events."""
        if event.key() == Qt.Key_C:  # 'C' key toggles charge visualization
            self.show_individual_charges = not self.show_individual_charges
            self.update()
