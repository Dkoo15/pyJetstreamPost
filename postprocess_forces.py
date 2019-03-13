import numpy as np
import jtstrm_surface

from load_parameters import *
from component import Component

''' usage: python postprocess_forces.py

    Loads the surfaces and computes the lift, drag, and momment forces
    Components can be provided for a more detailed breakdown.
'''
body = None

for i, surface in enumerate(files):
    jtstrm_surface.load(folder + surface)

    if body in components:
        components.remove(body)

    for comp in components:
        comp.zonelist = jtstrm_surface.filter(comp.box, False)
 
    body = jtstrm_surface.excluded_component(components, new_name='Fuselage')
    if body is not None:
        components.append(body)

    for comp in components:
        comp.clear()
        _ = jtstrm_surface.surface_force(component=comp,
                                         force_variables=None, cg=cg)
     
    aoa = jtstrm_surface.angle_of_attack()*np.pi/180
    lvec = np.array([-np.sin(aoa), 0, np.cos(aoa)])
    dvec = np.array([np.cos(aoa), 0, np.sin(aoa)])

    print('== Force Analysis: %s ==' % labels[i])
    print('Angle of Attack \t= %.3f' % (aoa*180/np.pi))   
    print('Reference Area \t= %.3f' % sref)

    total_area = sum([c.area for c in components])
    print('Surface Area \t= %.3f' % total_area)

    cl = sum([c.total_force(lvec) for c in components])
    print('CL  \t= %.4f' % (cl/sref))
    cd = sum([c.total_force(dvec) for c in components])
    print('CD  \t= %.1f (counts)' % (cd/sref*10000))
    cm = sum([c.moment for c in components])
    print('CM  \t= %.4f' % (cm/sref))
    cy = sum([c.total_force([0, -1, 0]) for c in components])
    print('CY  \t= %.4f' % (cy/sref))
    print('L/D \t= %.3f' % (cl/cd))
    print(' ')
    print('== Lift Breakdown==')
    for c in components:
        cl_breakdown = c.total_force(lvec)
        print('CL_%s \t= %.5f   (%5.2f%%)' % (c.name, cl_breakdown/sref,
                                              cl_breakdown/cl*100))

    print('== Drag Breakdown==')
    for c in components:
        cd_breakdown = c.total_force(dvec)
        print('CD_%s \t= %.5f   (%5.2f%%)' % (c.name, cd_breakdown/sref,
                                              cd_breakdown/cd*100))

    print('== Moment about %.4f,%.4f ==' % (cg[0], cg[2]))
    for c in components:
        print('CM_%s \t= %.5f   (%5.2f%%)' % (c.name, c.moment/sref,
                                              c.moment/cm*100))

    print('== Side Force Breakdown ==')
    for c in components:
        cy_breakdown = c.total_force((0, -1, 0))
        print('CY_%s \t= %.5f   (%5.2f%%)' % (c.name, cy_breakdown/sref,
                                              cy_breakdown/cy*100))

    print('== Drag Pressure/Friction Split ==')
    cdpres = sum([c.pressure_force(dvec) for c in components])
    cdfric = sum([c.friction_force(dvec) for c in components])
    print('CD_pres \t= %.5f   (%5.2f%%)' % (cdpres/sref, cdpres/cd*100))
    print('CD_fric \t= %.5f   (%5.2f%%)' % (cdfric/sref, cdfric/cd*100))
    print('   ')
