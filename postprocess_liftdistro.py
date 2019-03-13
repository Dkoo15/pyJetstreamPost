import numpy as np
import tecplot

from aero_sections import *
from load_parameters import *

''' usage: python postprocess_liftdistro.py

    Input files and settings are specified in load_parameters.py
    This script loads the files and calls functions in aero_sections to
    get the twist and thickness info for specified spanwise sections. Data
    is written in Tecplot form and can be plotted here.

    functions:
        create_lift_dataset - creates dataset from section list
'''


def create_lift_dataset(extracted_section_list, labels):
    '''Create a tecplot dataset with the twist and thickness distribution

    Arguments:
        extracted_section_list - (list of SectionData objects)
        labels - (list of str)
    '''

    new_dataset = tecplot.active_frame().create_dataset(
                                                        'Spanwise lift',
                                                        ['y/b', 'Fn'],
                                                        True)

    for name in labels:
        sub_section_list = [s for s in extracted_section_list
                            if s.label == name]

        n = len(sub_section_list)
        yb = np.zeros(n)
        lift = np.zeros(n)

        for i, section in enumerate(sub_section_list):
            yb[i] = section.location
            lift[i] = section.section_lift()

        zone = new_dataset.add_ordered_zone(name=name, shape=(n, 1, 1))
        zone.values('y/b')[:] = yb[:]
        zone.values('Fn')[:] = lift[:]

    return new_dataset


def elliptical_lift(dataset, cl, span, npts=100, wing_offset=0.0):
    lift = np.zeros(npts)
    y = np.linspace(0, 1, npts)
    dy = 1/npts

    den = 0
    for i in range(npts-1):
        den = den + np.sqrt(1-(y[i])**2)*dy
    lroot = cl/den
    for i in range(npts):
        lift[i] = lroot*np.sqrt(1-(y[i])**2)

    y = y*(1-wing_offset) + wing_offset

    elliptical = dataset.add_ordered_zone(name='lift_elliptical',
                                          shape=(npts, 1, 1)
                                          )
    elliptical.values('Fn')[:] = lift[:]/abs(span)
    elliptical.values('y/b')[:] = y[:]

    return elliptical


# == Main script begins
colors = [tecplot.constant.Color.Red,
          tecplot.constant.Color.Blue,
          tecplot.constant.Color.Green,
          tecplot.constant.Color.Purple]

section_list = []

for i, surface in enumerate(files):
    init_surface(folder + surface, box=wing.box)
    lift, span = force_sections(labels[i], section_list, cg)

dataset = create_lift_dataset(section_list, labels)
elliptical_zone = elliptical_lift(dataset, lift, span)

zonelist = list(dataset.zones())
main = tecplot.active_frame()
main.plot_type = tecplot.constant.PlotType.XYLine
main.show_border = False
main.show_header = False
main.plot().delete_linemaps()
pltxy = main.plot()
pltxy.axes.axis_mode = tecplot.constant.AxisMode.Independent
pltxy.axes.viewport.left = 10
pltxy.axes.viewport.right = 90
spanwise = dataset.variable('y/b')
lift = dataset.variable('Fn')
for i, zone in enumerate(zonelist):
    lmap = pltxy.add_linemap(name=zone.name, zone=zone, x=spanwise, y=lift)
    lmap.line.color = colors[i]
    lmap.line.line_thickness = 0.2
    if 'elliptical' in zone.name:
        lmap.line.line_pattern = tecplot.constant.LinePattern.Dashed
    lmap.x_axis.max = 1.0
    lmap.x_axis.min = 0.0
    lmap.y_axis.max = 0.7
    lmap.y_axis.min = 0.0
    lmap.show_in_legend = tecplot.constant.LegendShow.Always
    lmap.y_axis.title.text = 'Lift'
    lmap.y_axis.title.title_mode = tecplot.constant.AxisTitleMode.UseText
    lmap.y_axis.ticks.length = 1
    lmap.y_axis.ticks.minor_length = 0.5
    lmap.y_axis.tick_labels.format.remove_leading_zeros = True
    pltxy.legend.show = True
    pltxy.legend.box.box_type = tecplot.constant.TextBox.None_
    pltxy.legend.box.margin = 0.05
    pltxy.legend.font.size = 2.5
    pltxy.legend.position = (85, 85)

print('Exporting image: lift_distro.png ...')
tecplot.export.save_png(folder+'lift_distro.png', IMGWIDTH)
tecplot.data.save_tecplot_ascii(folder + 'lift_distro.dat', dataset=dataset,
                                use_point_format=True)
