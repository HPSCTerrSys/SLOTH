# importing external repos
import sys
Diagnostics_path='../extern/ana_parflow-diagnostics_pythonheat'
sys.path.append(Diagnostics_path)
import Diagnostics
import IO as diagIO
TrafoRotPole_path='../extern/transform_rotated_pole'
sys.path.append(TrafoRotPole_path)
from rotate import undo_grid_rotation

