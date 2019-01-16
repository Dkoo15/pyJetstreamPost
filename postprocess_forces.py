import numpy as np
import jtstrm_surface

from load_parameters import *
from component import Component

''' usage: python postprocess_forces.py

    Loads the surfaces and computes the lift, drag, and momment forces
    Components can be provided for a more detailed breakdown.
'''

for i, surface in enumerate(files):
    jtstrm_surface.load(folder + surface)

    for comp in components:
        comp.zonelist = jtstrm_surface.filter(comp.box, False)
    other = jtstrm_surface.excluded_component(components, new_name='Fuselage')
    if other is not None:
        components.append(other)

    for comp in components:
        comp.clear()
        _ = jtstrm_surface.surface_force(component=comp,
                                         force_variables=None, cg=cg)

    aoa = jtstrm_surface.angle_of_attack()*np.pi/180
    lvec = np.array([-np.sin(aoa), 0, np.cos(aoa)])
    dvec = np.array([np.cos(aoa), 0, np.sin(aoa)])

    print('== Force Analysis: %s ==' % labels[i])
    print('Reference Area = %.3f' % sref)

    total_area = sum([c.area for c in components])
    print('Surface Area = %.3f' % total_area)

    cl = sum([c.total_force(lvec) for c in components])
    print('CL  = %.4f' % (cl/sref))
    cd = sum([c.total_force(dvec) for c in components])
    print('CD  = %.1f (counts)' % (cd/sref*10000))
    cm = sum([c.moment for c in components])
    print('CM  = %.4f' % (cm/sref))
    cy = sum([c.total_force([0, -1, 0]) for c in components])
    print('CY  = %.4f' % (cy/sref))
    print('L/D = %.3f' % (cl/cd))
    print(' ')
    print('== Lift Breakdown==')
    for c in components:
        cl_breakdown = c.total_force(lvec)
        print('CL_%s = %.5f   (%5.2f%%)' % (c.name, cl_breakdown/sref,
                                            cl_breakdown/cl*100))

    print('== Drag Breakdown==')
    for c in components:
        cd_breakdown = c.total_force(dvec)
        print('CD_%s = %.5f   (%5.2f%%)' % (c.name, cd_breakdown/sref,
                                            cd_breakdown/cd*100))

    print('== Moment about %.4f,%.4f ==' % (cg[0], cg[2]))
    for c in components:
        print('CM_%s = %.5f   (%5.2f%%)' % (c.name, c.moment/sref,
                                            c.moment/cm*100))

    print('== Side Force Breakdown ==')
    for c in components:
        cy_breakdown = c.total_force((0, -1, 0))
        print('CY_%s = %.5f   (%5.2f%%)' % (c.name, cy_breakdown/sref,
                                            cy_breakdown/cy*100))

    print('== Drag Pressure/Friction Split ==')
    cdpres = sum([c.pressure_force(dvec) for c in components])
    cdfric = sum([c.friction_force(dvec) for c in components])
    print('CD_pres = %.5f   (%5.2f%%)' % (cdpres/sref, cdpres/cd*100))
    print('CD_fric = %.5f   (%5.2f%%)' % (cdfric/sref, cdfric/cd*100))
    print('   ')
