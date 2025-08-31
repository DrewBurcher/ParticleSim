from PyQt5.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, 
                          QHBoxLayout, QPushButton, QLabel)
from PyQt5.QtCore import QTimer

from src.gui.gl_widget import GLWidget

class MainWindow(QMainWindow):
    def __init__(self, config=None, field_sources=None):
        super(MainWindow, self).__init__()
        self.setWindowTitle('3D Particle Simulation')

        # Create central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Create info panel
        info_layout = QHBoxLayout()
        self.time_label = QLabel('Time: 0.000000 ps')
        self.fps_label = QLabel('FPS: 0')
        info_layout.addWidget(self.time_label)
        info_layout.addStretch()
        info_layout.addWidget(self.fps_label)
        layout.addLayout(info_layout)

        # Create GL widget with configuration
        self.gl_widget = GLWidget(config=config, field_sources=field_sources)
        layout.addWidget(self.gl_widget)

        # Add control buttons
        button_layout = QHBoxLayout()
        
        pause_button = QPushButton('Pause/Resume')
        pause_button.clicked.connect(self.toggle_simulation)
        
        reset_button = QPushButton('Reset')
        reset_button.clicked.connect(self.reset_simulation)
        
        button_layout.addWidget(pause_button)
        button_layout.addWidget(reset_button)
        layout.addLayout(button_layout)

        # Setup timer for updating display
        self.update_timer = QTimer()
        self.update_timer.timeout.connect(self.update_display)
        self.update_timer.start(16)  # Update every 16ms (roughly 60 FPS)

    def update_display(self):
        if hasattr(self.gl_widget, 'universe') and self.gl_widget.universe:
            time_ps = self.gl_widget.universe.time
            self.time_label.setText(f'Time: {time_ps:.12f} s')
            self.fps_label.setText(f'FPS: {self.gl_widget.fps}')

    def toggle_simulation(self):
        self.gl_widget.universe.toggle_simulation()

    def reset_simulation(self):
        self.gl_widget.universe.reset()