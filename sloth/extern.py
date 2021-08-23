# importing external repos
import os
import sys
Diagnostics_path=f'{os.path.dirname(__file__)}/../extern/ana_parflow-diagnostics_pythonheat'
sys.path.append(Diagnostics_path)
import Diagnostics
import IO as diagIO
TrafoRotPole_path=f'{os.path.dirname(__file__)}/../extern/transform_rotated_pole'
sys.path.append(TrafoRotPole_path)
from rotate import undo_grid_rotation
from rotate import rotate_grid

