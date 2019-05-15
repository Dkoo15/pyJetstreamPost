import tecplot
import numpy as np
import jtstrm_surface

from load_parameters import *

''' usage: python postprocess_contour.py

    Loads 2 surfaces and plots a CP contour with overhead view, one surface
    mirrored This is different from the other post process scripts as it
    calls the jtstrm_surface module directly.
'''

if len(files) != 2:
    print("Requires 2 surfaces as input")
    raise SystemExit

jtstrm_surface.load(folder + files[0])
jtstrm_surface.mirror()
dataset = (tecplot.data.load_tecplot(
           filenames=folder + files[1],
           read_data_option=tecplot.constant.ReadDataOption.Append)
           )
print('-> Loaded File: %s.' % (folder + files[1]))
main = tecplot.active_frame()
main.add_text(labels[0], (20, 90), size=14)
main.add_text(labels[1], (70, 90), size=14)

main.plot_type = tecplot.constant.PlotType.Cartesian3D
main.plot().view.rotate_to_angles(0, 0, -90)
main.plot().axes.orientation_axis.show = False
main.height = 8.
main.width = 10.

main.plot().show_contour = True
contour = main.plot().contour(0)
contour.levels.reset_levels(np.linspace(-1, 1, 41))
for fmap in main.plot().fieldmaps():
    fmap.contour.contour_type = tecplot.constant.ContourType.Overlay
legend = contour.legend
legend.show = True
legend.vertical = False
legend.label_step = 4
legend.position = (95,  12)
legend.header_font.size = 1.5
legend.number_font.size = 1.5
legend.box.box_type = tecplot.constant.TextBox.None_
main.plot().view.fit()
main.plot().use_lighting_effect = False
main.plot().view.center()
print('Exporting image: cp_contour.png ...')
tecplot.export.save_png(folder+'cp_contour.png', IMGWIDTH)
tecplot.export.save_ps(folder+'cp_contour.ps')
