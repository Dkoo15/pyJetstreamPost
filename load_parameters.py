from component import Component
''' usage: python load_parameters.py

This script is for setting any global level variables for the
postprocess scripts to use later on.
'''

IMGWIDTH = 2400
EXPORT_PLOT = True

folder = 'crm/'

files = ['surfaceCPCF_001.dat', 'surfaceCPCF_358.dat']
labels = ['initial', 'optimized']

# Reference area for lift / drag calculation
sref = 3.407
cg = [1.2077, 0, 0.6452]

# Spanwise Sections- y/b locations
span_stations = [0.0235, 0.267, 0.557, 0.695, 0.828, 0.944]

# Indicate which y/b is where the winglet starts. Set >1.0 for no winglet
dihedral_stn = 1.1

# Indicate which y/b is where the wing starts
wing_offset = 0.0

wing = Component("Wing")
components = [wing]
