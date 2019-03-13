import numpy as np
import tecplot

from aero_sections import *
from load_parameters import *
from component import Component

''' usage: python postprocess_geometry.py

    Input files and settings are specified in load_parameters.py
    This script loads the files and calls functions in aero_sections to
    get the twist and thickness info for specified spanwise sections. Data
    is written in Tecplot form and can be plotted here.

    functions:
        create_spanwise_dataset - creates dataset from section list
'''


def create_spanwise_dataset(extracted_section_list, labels, tc_locations):
    '''Create a tecplot dataset with the twist and thickness distribution

    arguments:
        extracted_section_list - (list of SectionData objects)
        labels - (list of str)
        tc_locations - (list of float)
    '''
    variables = ['y/b', 'twist', 'tcmax']
    for x in tc_locations:
        variables.append('tc %4.2f' % x)

    new_dataset = (tecplot.active_frame().create_dataset(
                   'Spanwise twist and thickness distribution',
                   variables, True)
                   )

    for name in labels:
        sub_section_list = [s for s in extracted_section_list
                            if s.label == name]

        n = len(sub_section_list)
        yb = np.zeros(n)
        twist = np.zeros(n)
        tcmax = np.zeros(n)
        thick = np.zeros(n)
        for i, section in enumerate(sub_section_list):
            yb[i] = section.location
            twist[i] = section.twist/np.pi*180.0
            tcmax[i] = section.max_thickness()

        zone = new_dataset.add_ordered_zone(name=name, shape=(n, 1, 1))
        zone.values('y/b')[:] = yb[:]
        zone.values('twist')[:] = twist[:]
        zone.values('tcmax')[:] = tcmax[:]

        for x in tc_locations:
            for i, section in enumerate(sub_section_list):
                thick[i] = section.thickness(x)
            zone.values('tc %4.2f' % x)[:] = thick[:]

    return new_dataset

# == Main script begins

colors = [tecplot.constant.Color.Red,
          tecplot.constant.Color.Blue,
          tecplot.constant.Color.Green,
          tecplot.constant.Color.Purple]

symbols = [tecplot.constant.GeomShape.Del,
           tecplot.constant.GeomShape.Diamond,
           tecplot.constant.GeomShape.Circle]


section_list = []
tc_locations = [0.15, 0.90]

for i, surface in enumerate(files):
    init_surface(folder + surface, box=wing.box)
    geometry_sections(labels[i], section_list, dihedral_stn=dihedral_stn)

dataset = create_spanwise_dataset(section_list, labels, tc_locations)

# ================ Create Line Plots ================#
print('Exporting section data file: twist_thick_distro.dat ...')
tecplot.data.save_tecplot_ascii(folder + 'twist_thick_distro.dat',
                                dataset=dataset,
                                use_point_format=True)

if not EXPORT_PLOT:
    raise SystemExit

page = tecplot.active_page()
tframe = tecplot.active_frame()
frameheight = page.paper.dimensions[0]/2
framewidth = page.paper.dimensions[0]/2

# Twist plot first:
tframe.height = frameheight
tframe.width = framewidth
tframe.plot_type = tecplot.constant.PlotType.XYLine
tframe.show_border = False
tframe.show_header = False
tframe.position = (0, 0.25)
plot = tframe.plot()
plot.delete_linemaps()
plot.axes.axis_mode = tecplot.constant.AxisMode.Independent
plot.axes.viewport.left = 12
plot.axes.viewport.right = 95
plot.axes.viewport.top = 95
plot.axes.viewport.bottom = 10

for i, name in enumerate(labels):
    lmap = plot.add_linemap(name=name + ' twist', zone=dataset.zone(name),
                            x=dataset.variable('y/b'),
                            y=dataset.variable('twist')
                            )
    lmap.line.color = colors[i]
    lmap.line.line_thickness = 0.5
    lmap.show_in_legend = tecplot.constant.LegendShow.Always
    lmap.symbols.show = True
    lmap.symbols.symbol().shape = tecplot.constant.GeomShape.Circle
    lmap.symbols.size = 1
    lmap.symbols.fill_mode = tecplot.constant.FillMode.UseSpecificColor
    lmap.symbols.fill_color = colors[i]
    lmap.symbols.color = colors[i]

lmap.x_axis.max = 1.0
lmap.x_axis.min = 0.0
lmap.y_axis.max = 7.
lmap.y_axis.min = -7.
plot.legend.show = True
plot.legend.box.box_type = tecplot.constant.TextBox.None_
plot.legend.box.margin = 0.05
plot.legend.font.size = 2.5
plot.legend.position = (85, 85)
plot.show_symbols = True

# Thickness plot:
frame = page.add_frame()
frame.height = frameheight
frame.width = framewidth
frame.plot_type = tecplot.constant.PlotType.XYLine
frame.show_border = False
frame.show_header = False
frame.position = (framewidth, 0.25)
plot = frame.plot()
plot.delete_linemaps()
plot.axes.axis_mode = tecplot.constant.AxisMode.Independent
plot.axes.viewport.left = 12
plot.axes.viewport.right = 95
plot.axes.viewport.top = 95
plot.axes.viewport.bottom = 10
plot.show_symbols = True

for i, name in enumerate(labels):
    lmap = plot.add_linemap(name=name + ' max', zone=dataset.zone(name),
                            x=dataset.variable('y/b'),
                            y=dataset.variable('tcmax')
                            )
    lmap.line.color = colors[i]
    lmap.line.line_thickness = 0.2
    lmap.show_in_legend = tecplot.constant.LegendShow.Always
    lmap.symbols.show = False
    lmap.line.line_pattern = tecplot.constant.LinePattern.Dashed

    for k, x in enumerate(tc_locations):
        'tc %4.2f' % x
        var_name = 'tc'+str(tc_locations[k])
        lmap = plot.add_linemap(name=name + ' x/c=%d%%' % (100*x),
                                zone=dataset.zone(name),
                                x=dataset.variable('y/b'),
                                y=dataset.variable('tc %4.2f' % x)
                                )
        lmap.line.color = colors[i]
        lmap.symbols.show = True
        lmap.symbols.symbol().shape = symbols[k]
        lmap.symbols.size = 1
        lmap.symbols.fill_mode = tecplot.constant.FillMode.UseSpecificColor
        lmap.symbols.fill_color = colors[i]
        lmap.symbols.color = colors[i]
        lmap.line.line_thickness = 0.2
        lmap.show_in_legend = tecplot.constant.LegendShow.Always

    lmap.x_axis.max = 1.0
    lmap.x_axis.min = 0.0
    lmap.y_axis.max = 0.22
    lmap.y_axis.min = 0.01
    lmap.y_axis.title.text = 'Thickness-to-chord ratio'
    lmap.y_axis.title.font.size = 3
    lmap.y_axis.title.title_mode = tecplot.constant.AxisTitleMode.UseText
    lmap.y_axis.tick_labels.format.remove_leading_zeros = True
    plot.legend.show = True
    plot.legend.box.box_type = tecplot.constant.TextBox.None_
    plot.legend.box.margin = 0.05
    plot.legend.font.size = 2.5
    plot.legend.position = (85, 85)

print('Exporting image: twist_and_thickness_distro.png ...')
tecplot.export.save_png(folder+'twist_and_thickness_distro.png', IMGWIDTH,
                        region=tecplot.constant.ExportRegion.AllFrames)
