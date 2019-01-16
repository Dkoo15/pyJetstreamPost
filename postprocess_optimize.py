import tecplot
import numpy as np
import sys

''' usage: python postprocess_optimize.py

    This script loads optimization history files and generates plots for
    Optimality, Feasibility, Merit Function etc. This is mainly customized
    for lift-constrainted drag minimization of aerodynamic geometries.

    functions:
        open_results_check
        optimize_his
'''


def open_results_check(filename):
    '''Open .ohis file format with checks'''
    try:
        fp = open(filename, 'r')
        return fp
    except IOError:
        return None


def optimize_his(filename):
    '''open optimize.his and parse file with checks'''
    try:
        f = open(filename, 'r')
    except IOError:
        print("No optimize.his file found!!")
        raise SystemExit

    raw = f.read()
    f.close()
    raw = raw.replace('(', '')
    raw = raw.replace(')', '')
    data = raw.splitlines()

    iteration_lines = []
    for line in data:
        if len(line) > 0 and line != ' ':
            words = line.split()
            if words[0].isdigit():
                if int(words[0]) % 100 == 0 and len(words) < 8:
                    continue
                iteration_lines.append(line.split())
    # fill data
    nfun = len(iteration_lines)
    optim = np.zeros(nfun)
    feasi = np.zeros(nfun)
    merit = np.zeros(nfun)
    feval = [0]*nfun

    for i in range(nfun):
        nums = iteration_lines[i]
        j = 3
        if nums[0] == '0':
            j = 2
        feval[i] = int(nums[j])
        feasi[i] = float(nums[j+1])
        optim[i] = float(nums[j+2])
        merit[i] = float(nums[j+3])

    return nfun, feval, feasi, optim, merit


# == Main script begins

colors = [tecplot.constant.Color.Red,
          tecplot.constant.Color.Blue,
          tecplot.constant.Color.Green,
          tecplot.constant.Color.Purple,
          tecplot.constant.Color.Cyan]

if sys.version_info > (3, 0):  # Python3
    folder = input('Folder containing files: ')
    sref = float(input('Enter reference area: '))
else:  # Python2
    folder = raw_input('Folder containing files: ')
    sref = float(raw_input('Enter reference area: '))

nfun, feval, feasi, optim, merit = optimize_his(folder + '/optimize.his')
dataset = tecplot.active_frame().create_dataset('Optimization History',
                                                ['Iteration',
                                                 'Feasibility',
                                                 'Optimality'],
                                                True)
zone = dataset.add_ordered_zone(name=folder, shape=(nfun, 1))
zone.values('Iteration')[:] = range(nfun+1)
zone.values('Feasibility')[:] = feasi
zone.values('Optimality')[:] = optim

page = tecplot.active_page()
width = page.paper.dimensions[0]/2
height = width

# ========Feasibility & Optimality Plot==========#
frame = page.add_frame()
frame.name = 'opfeas_frame'
frame.plot_type = tecplot.constant.PlotType.XYLine
frame.show_border = False
frame.show_header = False
frame.height = height
frame.width = width
frame.position = (0.0, 2.5)

plot = frame.plot()
plot.delete_linemaps()
plot.axes.axis_mode = tecplot.constant.AxisMode.Independent
plot.axes.viewport.left = 12
plot.axes.viewport.right = 95
plot.axes.viewport.top = 95
plot.axes.viewport.bottom = 10

lmap = plot.add_linemap(name='feasibility', zone=zone,
                        x=dataset.variable('Iteration'),
                        y=dataset.variable('Feasibility')
                        )
lmap.line.color = tecplot.constant.Color.Blue
lmap.line.line_thickness = 0.8
lmap.x_axis.fit_range()
lmap.show_in_legend = tecplot.constant.LegendShow.Always

lmap = plot.add_linemap(name='optimality', zone=zone,
                        x=dataset.variable('Iteration'),
                        y=dataset.variable('Optimality')
                        )
lmap.line.color = tecplot.constant.Color.Green
lmap.line.line_thickness = 0.8
lmap.y_axis.log_scale = True
lmap.y_axis.fit_range()
lmap.y_axis.max = 1e-1
lmap.y_axis.min = 1e-7
lmap.y_axis.title.text = 'Feasibility, Optimality'
lmap.y_axis.title.font.size = 3
lmap.y_axis.title.title_mode = tecplot.constant.AxisTitleMode.UseText
lmap.y_axis.title.offset = 7
lmap.show_in_legend = tecplot.constant.LegendShow.Always

plot.legend.show = True
plot.legend.box.box_type = tecplot.constant.TextBox.None_
plot.legend.box.margin = 0.05
plot.legend.font.size = 2.5
plot.legend.position = (85, 85)

# ========Merit Function (Drag) Plot=========#
dataset.add_variable('cd1')
drag = np.zeros(nfun)

ohis = open_results_check(folder + '/results-001.ohis')
multipoint = ohis is not None

if not multipoint:
    ohis = open_results_check(folder+'/results.ohis')
    if ohis is None:
        print('No results.ohis filename, plotting merit function')
        drag = merit/sref*1e4
        map_name = 'drag'
    else:
        data = ohis.read().splitlines()
        ohis.close()
        del data[0]
        iterations = [data[i] for i in feval]
        lift = np.zeros(nfun)
        for i in range(nfun):
            it = iterations[i].split()
            lift[i] = float(it[4])/sref
            drag[i] = float(it[5])/sref*1e4

        dataset.add_variable('cl')
        zone.values('cl')[:] = lift
        mach = float(it[1])
        map_name = 'M=%.2f, CL=%.2f' % (mach, lift[nfun-1])

else:
    print('Loading multipoint results-001.ohis filename')
    data = ohis.read().splitlines()
    del data[0]
    ohis.close()
    iterations = [data[i] for i in feval]
    for i in range(nfun):
        it = iterations[i].split()
        drag[i] = float(it[5])/sref*1e4

    mach = float(it[1])
    cl = float(it[4])/sref
    map_name = 'M=%.2f, CL=%.2f' % (mach, cl)

zone.values('cd1')[:] = drag

frame = page.add_frame()
frame.name = 'cd_frame'
frame.plot_type = tecplot.constant.PlotType.XYLine
frame.show_border = False
frame.show_header = False
frame.height = height
frame.width = width
frame.position = (width, 2.5)

plot = frame.plot()
plot.delete_linemaps()
plot.axes.axis_mode = tecplot.constant.AxisMode.Independent
plot.axes.viewport.left = 12
plot.axes.viewport.right = 95
plot.axes.viewport.top = 95
plot.axes.viewport.bottom = 10
plot.axes.viewport.right = 80

cdmap = plot.add_linemap(name=map_name, zone=zone,
                         x=dataset.variable('Iteration'),
                         y=dataset.variable('cd1')
                         )
cdmap.line.color = tecplot.constant.Color.Red
cdmap.line.line_thickness = 0.8
cdmap.x_axis.fit_range()
cdmap.x_axis.min = -1
cdmap.y_axis.ticks.length = 1
cdmap.y_axis.ticks.minor_length = 0.5
cdmap.y_axis.title.font.size = 3
cdmap.y_axis.title.offset = 7
cdmap.y_axis.max = drag.max()
cdmap.y_axis.min = drag.min()*0.99
cdmap.y_axis_index = 0
cdmap.y_axis.title.text = 'Drag (counts)'
cdmap.y_axis.title.title_mode = tecplot.constant.AxisTitleMode.UseText

if singlepoint:
    cdmap.y_axis.line.color = tecplot.constant.Color.Red
    cdmap.y_axis.tick_labels.color = tecplot.constant.Color.Red
    cdmap.y_axis.title.color = tecplot.constant.Color.Red
    clmap = plot.add_linemap(name='cl_map', zone=zone,
                             x=dataset.variable('Iteration'),
                             y=dataset.variable('cl')
                             )
    clmap.line.color = tecplot.constant.Color.Blue
    clmap.line.line_thickness = 0.8
    clmap.y_axis_index = 1
    clmap.y_axis.title.font.size = 3
    clmap.y_axis.title.offset = 7
    clmap.y_axis.tick_labels.format.remove_leading_zeros = True
    clmap.y_axis.min = lift.max()*1.02
    clmap.y_axis.min = lift.min()*0.98
    clmap.y_axis.line.color = tecplot.constant.Color.Blue
    clmap.y_axis.tick_labels.color = tecplot.constant.Color.Blue
    clmap.y_axis.title.color = tecplot.constant.Color.Blue
    clmap.y_axis.title.text = 'Lift Coefficient'
    clmap.y_axis.title.title_mode = tecplot.constant.AxisTitleMode.UseText

if multipoint:
    j = 2
    cdmap.y_axis.max = 310
    cdmap.y_axis.min = 220
    ohis = open_results_check(folder + '/results-%03d.ohis' % j)
    while ohis is not None:
        print('Loaded multipoint results-%03d.ohis filename' % j)
        data = ohis.read().splitlines()
        del data[0]
        ohis.close()
        iterations = [data[i] for i in feval]
        for i in range(nfun):
            it = iterations[i].split()
            drag[i] = float(it[5])/sref*1e4  # drag is already allocated

        mach = float(it[1])
        cl = float(it[4])/sref
        map_name = 'M=%.2f, CL=%.2f' % (mach, cl)

        var_name = 'cd%02d' % j
        dataset.add_variable(var_name)
        zone.values(var_name)[:] = drag
        cdmap = plot.add_linemap(name=map_name, zone=zone,
                                 x=dataset.variable('Iteration'),
                                 y=dataset.variable(var_name)
                                 )
        cdmap.line.color = colors[j-1]
        cdmap.line.line_thickness = 0.8
        cdmap.y_axis_index = 0
        j = j + 1

    plot.legend.show = True
    plot.legend.box.box_type = tecplot.constant.TextBox.None_
    plot.legend.box.margin = 0.05
    plot.legend.font.size = 2.5
    plot.legend.position = (85, 85)

page.delete_frame(main)
print('Exporting image: optimization_history.png ...')
tecplot.export.save_png(folder+'/optimization_history.png', 2000,
                        region=tecplot.constant.ExportRegion.AllFrames)
