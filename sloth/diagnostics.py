# importing external repo ParFlowDiagnostics
import os
import sys
Diagnostics_path=f'{os.path.dirname(__file__)}/../extern/ParFlowDiagnostics'
sys.path.append(Diagnostics_path)
import Diagnostics
import IO as diagIO
