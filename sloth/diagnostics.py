# importing external repo ana_parflow-diagnostics_pythonheat
import os
import sys
Diagnostics_path=f'{os.path.dirname(__file__)}/../extern/ana_parflow-diagnostics_pythonheat'
sys.path.append(Diagnostics_path)
import Diagnostics
import IO as diagIO
