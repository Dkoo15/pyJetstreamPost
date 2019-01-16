import numpy as np
import tecplot

from aero_sections import *
from load_parameters import *
from component import Component

''' usage: python postprocess_sections.py

    Input files and settings are specified in load_parameters.py
    This script loads the files and calls functions in aero_sections to
    get the x/c, z/c, cp, and cf data for specified spanwise sections. Data
    is written in Tecplot form and can be plotted here.

    functions:
        create_section_dataset - creates dataset from section list
'''


def create_section_dataset(extracted_section_list, variables):
    '''Create a tecplot dataset with all the sections

    Arguments:
        extracted_section_list - (list of SectionData objects
        variables - (tuple of str)
    '''
    new_dataset = tecplot.active_frame().create_dataset('Surface CPCF',
                                                        variables, True)
    for section in extracted_section_list:
        sec_name = '%s_y/b=%.3f' % (section.label, section.location)
        z = new_dataset.add_fe_zone(tecplot.constant.ZoneType.FELineSeg,
                                    name=sec_name,
                                    num_points=len(section.data),
                                    num_elements=len(section.nodemap))
        for i, v in enumerate(variables):
            z.values(v)[:] = section.data[:, i]
        z.nodemap[:] = section.nodemap

    return new_dataset

# == Main script begins

colors = [tecplot.constant.Color.Red,
          tecplot.constant.Color.Blue,
          tecplot.constant.Color.Green,
          tecplot.constant.Color.Purple]

section_list = list()

for i, surface in enumerate(files):
    init_surface(folder + surface)
    cpcf_sections(labels[i], span_stations, ('Cp', 'Cf'), section_list,
                  dihedral_stn)

dataset = create_section_dataset(section_list, ('x/c', 'z/c', 'Cp', 'Cf'))

print('Exporting section data file: cp_sections.dat ...')
tecplot.data.save_tecplot_ascii(folder+'cp_sections.dat', dataset=dataset,
                                use_point_format=True)

print('-----Completed Extraction of Slices-----')

if not EXPORT_PLOT:
    raise SystemExit

# Create Line Plots from extracted section zones
framelist = []
page = tecplot.active_page()
main = tecplot.active_frame()

for stn in span_stations:
    frame = page.add_frame()
    frame.name = 'y/b='+str(stn)
    frame.plot_type = tecplot.constant.PlotType.XYLine
    frame.show_border = False
    frame.show_header = False
    frame.plot().delete_linemaps()
    frame.add_text('%.2f %% span' % (100*stn), (37, 90), size=10)
    framelist.append(frame)

    plot = frame.plot()
    plot.axes.axis_mode = tecplot.constant.AxisMode.Independent
    plot.axes.viewport.top = 90
    plot.axes.viewport.bottom = 10
    plot.axes.viewport.left = 12

    for i, name in enumerate(labels):
        section_name = '%s_y/b=%.3f' % (labels[i], stn)
        zone = dataset.zone(section_name)

        af_map = plot.add_linemap(name=name, zone=zone,
                                  x=dataset.variable('x/c'),
                                  y=dataset.variable('z/c')
                                  )
        af_map.line.color = colors[i]
        af_map.line.line_thickness = 0.2
        af_map.y_axis_index = 0
        af_map.y_axis.max = 0.8
        af_map.y_axis.min = -0.12
        af_map.y_axis.title.position = 20
        af_map.show_in_legend = tecplot.constant.LegendShow.Never
        af_map.y_axis.title.offset = 4
        af_map.y_axis.tick_labels.format.remove_leading_zeros = True
        af_map.y_axis.ticks.length = 1

        cp_map = plot.add_linemap(name=name, zone=zone,
                                  x=dataset.variable('x/c'),
                                  y=dataset.variable('Cp')
                                  )
        cp_map.line.color = colors[i]
        cp_map.line.line_thickness = 0.2
        cp_map.y_axis_index = 1
        cp_map.y_axis.reverse = True
        cp_map.y_axis.max = 1.2
        cp_map.y_axis.min = -1.25
        cp_map.x_axis.max = 1.05
        cp_map.x_axis.min = -0.05
        cp_map.show_in_legend = tecplot.constant.LegendShow.Always
        cp_map.y_axis.tick_labels.format.show_decimals_on_whole_numbers = True

    plot.legend.show = True
    plot.legend.box.box_type = tecplot.constant.TextBox.None_
    plot.legend.box.margin = 0.05
    plot.legend.font.size = 2.0
    plot.legend.position = (70, 40)

main.move_to_bottom()
framesPerRow = (len(framelist)+1)//2
frameHght = page.paper.dimensions[1]/2
frameWdth = page.paper.dimensions[0]/framesPerRow
for i, f in enumerate(framelist):
    frameX = (i % framesPerRow)*frameWdth
    frameY = (i//framesPerRow)*frameHght
    f.position = (frameX, frameY)
    f.height = frameHght
    f.width = frameWdth

print('Exporting section image: section_and_cp_distro.png ...')
tecplot.export.save_png(folder+'section_and_cp_distro.png', IMGWIDTH,
                        region=tecplot.constant.ExportRegion.AllFrames)
