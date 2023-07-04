# importing external repo colormaps
# https://github.com/pratiman-91/colormaps
import os
import sys
colormaps_path=f'{os.path.dirname(__file__)}/../extern/colormaps'
sys.path.append(colormaps_path)
import colormaps
