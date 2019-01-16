import numpy as np
import tecplot
import jtstrm_plot3d
import sys
import os

''' usage: python convert_plot3d.py

    This script loads the calls funtions in jtstrm_plot3d to load a plot3d
    flow solution. It extracts the surface data for the desired variables
    and outputs them to a new tecplot surface dataset. On startup, the user
    is prompted with a folder containing a grid file (.g) and a plot3d solution
    file (.q)
'''


def create_surface_dataset(extracted_surface_list, label, variables):
    '''Create a tecplot dataset with all the sections

    arguments:
        extracted_section_list - (list of SectionData objects)
        label - (str) dataset name
        variables - (list of str)
    '''
    new_dataset = tecplot.active_frame().create_dataset(label, variables, True)
    for boundary in extracted_surface_list:
        srf_name = 'S0%d%d' % (boundary.block+1, boundary.face)
        dimensions = (boundary.imax, boundary.jmax, 1)
        z = new_dataset.add_ordered_zone(name=srf_name, shape=dimensions)

        for i, v in enumerate(variables):
            z.values(v)[:] = boundary.data[:, i]

    return new_dataset

# == Main script begins
if sys.version_info > (3, 0):  # Python3
    folder = input('Folder containing files: ')
else:  # Python2
    folder = raw_input('Folder containing files: ')

listing = os.listdir(folder)

q_files = [f for f in listing if f.endswith('.q')]
g_files = [f for f in listing if f.endswith('.g')]
con_files = [f for f in listing if f.endswith('.con')]
if (len(g_files) != 1 or len(q_files) != 1):
    print("Error in folder %s, skipping..." % folder)
    raise SystemExit

grid = folder + '/' + g_files[0]
sol = folder + '/' + q_files[0]
conn = folder + '/' + con_files[0]

jtstrm_plot3d.load_solution(grid, sol)
coordinates = ('X', 'Y', 'Z')
coefficients = ('Cp', 'Cf', 'Cfx', 'Cfy', 'Cfz')
surface_list = jtstrm_plot3d.get_aero_surfaces(coordinates,
                                               coefficients, conn)

surface_dataset = create_surface_dataset(surface_list, 'Surface CPCF',
                                         variables=coordinates + coefficients)
print('Exporting surface data file: surfaceCPCF_new.dat ...')
tecplot.data.save_tecplot_ascii(folder + '/surfaceCPCF_new.dat',
                                dataset=surface_dataset, use_point_format=True)
