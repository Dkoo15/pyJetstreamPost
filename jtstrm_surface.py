import tecplot
import numpy as np
import block_utils
from component import Component

''' This module is for loading analyzing surface data files from Jetstream or
    created from convert_plot3d.py

    functions:
        load - loads a file into the dataset
        filter - look for all zones within a specified box
        angle_of_attack - returns reference angle of attack
        span - returns y-extents of the dataset
        dihedral_section - general arbitrary cross section
        cross_section - takes a section across the surface
        zonal_forces - loops through and sums forces for a zone
        surface_forces - computes and sums forces on this surface and places
                         the result into the component class
        mirror - flips surface across x-z plane
        excluded component - creates a final component, one which is not
                             included in any other provided component
                             (often the fuselage).
'''

_dataset = None


def load(filename):
    ''' Load tecplot surface file into the workspace'''
    global _dataset
    try:
        _dataset = (tecplot.data.load_tecplot(
                    filenames=filename,
                    read_data_option=tecplot.constant.ReadDataOption.Replace)
                    )
    except (tecplot.exception.TecplotSystemError,
            tecplot.exception.TecplotTypeError
            ):
        print('-> Error loading file %s' % filename)
        raise SystemExit

    tecplot.active_frame().plot_type = tecplot.constant.PlotType.Cartesian3D
    print('-> Loaded File: %s.' % filename)


def filter(box, activate=True):
    ''' Filter surface zones within a given box

    arguments:
        box - (numpy array of floats) format:[[x1, y1, z1], [x2, y2, z2]]
    optional arguments
        activate - (boolean) set only these active?
    '''
    global _dataset
    zonelist = list(_dataset.zones())
    if box is None:
        return zonelist

    xmin = box[0][0]
    xmax = box[1][0]
    ymin = box[0][1]
    ymax = box[1][1]
    zmin = box[0][2]
    zmax = box[1][2]
    subzones = []

    for zone in zonelist:
        if (zone.values('Y').min() > ymin and
                zone.values('Y').max() < ymax):
                if (zone.values('X').min() > xmin and
                        zone.values('X').max() < xmax):
                        if (zone.values('Z').min() > zmin and
                                zone.values('Z').max() < zmax):
                                subzones.append(zone)

    if (activate):
        _, = tecplot.active_frame().active_zones(subzones)

    return subzones


def angle_of_attack(default=0.0):
    ''' Load the angle of attack from the aux_data'''
    global _dataset
    try:
        alpha = float(_dataset.aux_data['Common.AngleOfAttack'])

    except (tecplot.exception.TecplotKeyError):
        alpha = default

    return alpha


def span():
    ''' Compute the y-coordinate of the wing root, wing tip, and the span'''
    global _dataset
    ymin = _dataset.variable('Y').min()
    ymax = _dataset.variable('Y').max()
    if ymax < 1e-6:  # Some files have negative span
        tip = ymin + 0.05
        root = ymax + 0.05*tip
    else:
        tip = ymax - 0.05
        root = ymin + 0.05*tip

    span = abs(tip-root)
    return tip, root, span


def dihedral_section(y, dihedral=None, variables=None, get_nodemap=False,
                     untwist=True):
    '''Compute a cross section with dihedral

    arguments:
        y - (float)
    optional arguments:
        variables - (list of str)
        get_nodemap - (boolean)
    '''
    global _dataset
    # take section slice at y+dy and y-dy
    dy = 0.01
    section_pos = (tecplot.data.extract.extract_slice(
                    origin=(0, y+dy, 0), normal=(0, 1, 0),
                    source=tecplot.constant.SliceSource.SurfaceZones,
                    dataset=_dataset)
                   )
    section_neg = (tecplot.data.extract.extract_slice(
                    origin=(0, y-dy, 0), normal=(0, 1, 0),
                    source=tecplot.constant.SliceSource.SurfaceZones,
                    dataset=_dataset)
                   )

    x = section_pos.values('X').as_numpy_array()
    ip = x.argmin()  # leading edge at y+dy
    xp = x[ip]
    zp = section_pos.values('Z')[ip]
    x = section_neg.values('X').as_numpy_array()
    im = x.argmin()  # leading edge at y-dy
    xm = x[im]
    zm = section_neg.values('Z')[im]

    dy = 2*dy
    dz = zp-zm
    d = np.sqrt(dy**2 + dz**2)
    normal = (0, dy/d, dz/d)      # normal of plane to extract slice
    gamma = np.arctan(abs(dz/dy))  # this is used later to rotate back
    origin = ((xm+xp)/2, y, (zm+zp)/2)  # origin of rotation

    _dataset.delete_zones([section_neg, section_pos])  # clean up

    if dihedral is not None:  # Override with given dihedral
        gamma = dihedral*np.pi/180.0
        normal = (0, np.cos(gamma), np.sin(gamma))

    af, nmap, twst, gamma, chrd, = cross_section(origin, normal, gamma,
                                                 variables, get_nodemap,
                                                 untwist)

    return af, nmap, twst, gamma, chrd


def cross_section(origin, normal=(0, 1, 0), gamma=0.0,
                  variables=None, get_nodemap=False, untwist=True):
    '''Take a section on the surface, getting any variables requested.
    Returns the airfoil coordinates, nodemap and other variables
    (if requested), twist, dihedral, and chord of the section

    arguments:
       origin - (3-tuple of floats)
    optional arguments:
       normal - (3-tuple of floats)
       gamma - (float)
       variables - (list of str)
       get_nodemap - (boolean)
       untwist - (boolean)
    '''
    global _dataset
    try:
        section = (tecplot.data.extract.extract_slice(origin=origin,
                   normal=normal,
                   source=tecplot.constant.SliceSource.SurfaceZones,
                   dataset=_dataset)
                   )
    except (tecplot.exception.TecplotLogicError):
        return None, None, 0.0, 0.0, 0.0

    x = section.values('X').as_numpy_array()
    z = section.values('Z').as_numpy_array()

    i = x.argmin()  # leading Edge
    j = x.argmax()  # trailing Edge
    xc = (x - x[i])
    zc = (z - z[i])

    if (gamma > np.pi/180):  # if dihedral is significant, transform
        y = section.values('Y').as_numpy_array()
        yc = y - y[i]
        ztemp = zc*np.cos(gamma) + yc*np.sin(gamma)
        zc = ztemp

    dx = xc[i] - xc[j]
    dz = zc[i] - zc[j]
    chord = np.sqrt(dx*dx + dz*dz)  # compute and Scale by chord
    xc = xc/chord
    zc = zc/chord
    twist = np.arctan(dz/dx)  # compute Twist

    if (untwist):   # untwist the section
        xtemp = xc*np.cos(twist) + zc*np.sin(twist)
        ztemp = -xc*np.sin(twist) + zc*np.cos(twist)
        xc = xtemp
        zc = ztemp

    if variables is None:
        m = 2
    else:
        m = len(variables) + 2

    n = section.num_points
    data = np.zeros((n, m), dtype=np.float32)
    data[:, 0] = xc[:]
    data[:, 1] = zc[:]
    if variables is not None:
        for j, v in enumerate(variables, 2):
            data[:, j] = section.values(v).as_numpy_array()

    if get_nodemap is False:
        nodemap = None
    else:
        nodemap = section.nodemap[:]  # Change this

    _dataset.delete_zones(section)  # clean up
    twist = -twist
    dihedral = gamma

    return data, nodemap, twist, dihedral, chord


def zonal_forces(zone, cg):
    ''' Compute the forces for a surface zone. Returns the pressure force,
    friction forces, moment about cg, and area, summed for this entire zone.

    arguments:
        zone - (tecplot zone)
        cg  - (list of float) (1-by-3)
        store - (boolean)
    '''
    n = zone.num_points
    imax, jmax, _, = zone.dimensions

    xyz = np.zeros((n, 3))
    xyz[:, 0] = zone.values('X').as_numpy_array()
    xyz[:, 1] = zone.values('Y').as_numpy_array()
    xyz[:, 2] = zone.values('Z').as_numpy_array()
    cp = zone.values('Cp').as_numpy_array()
    cfx = zone.values('Cfx').as_numpy_array()
    cfy = zone.values('Cfy').as_numpy_array()
    cfz = zone.values('Cfz').as_numpy_array()

    fx = np.zeros(n, dtype=np.float32)
    fy = np.zeros(n, dtype=np.float32)
    fz = np.zeros(n, dtype=np.float32)

    xi = block_utils.surface_metrics(imax, jmax, xyz)
    pres = np.zeros(3)
    fric = np.zeros(3)
    area = 0
    mmnt = 0

    ind = 0
    for j in range(jmax):
        if (j == 0 or j == jmax - 1):
            deta = 0.5
        else:
            deta = 1.

        for i in range(imax):
            if (i == 0 or i == imax - 1):
                dxi = 0.5
            else:
                dxi = 1.

            nrm = np.sqrt(xi[ind, 0]**2 + xi[ind, 1]**2 + xi[ind, 2]**2)
            fp = cp[ind]*xi[ind, :]
            ff = nrm*np.array([cfx[ind], cfy[ind], cfz[ind]])
            frc = fp + ff
            fx[ind] = frc[0]/nrm
            fy[ind] = frc[1]/nrm
            fz[ind] = frc[2]/nrm
            arm_x = -(xyz[ind, 0] - cg[0])
            arm_z = (xyz[ind, 2] - cg[2])

            area = area + dxi*deta*nrm
            pres[:] = pres[:] + fp*dxi*deta
            fric[:] = fric[:] + ff*dxi*deta
            mmnt = mmnt + (frc[2]*arm_x + frc[0]*arm_z)*dxi*deta

            ind = ind + 1

    if 'Fx' in _dataset.variable_names:
        zone.values('Fx')[:] = fx
    if 'Fy' in _dataset.variable_names:
        zone.values('Fy')[:] = fy
    if 'Fz' in _dataset.variable_names:
        zone.values('Fz')[:] = fz

    return pres, fric, mmnt, area


def surface_force(component=None, force_variables=None,
                  cg=(0, 0, 0)):
    ''' Compute the forces for a surface. Returns the forces, and
    moment about cg, and area for each component.

    optional arguments:
        components - (Component)
        force_variables - (list of str)
        cg - (list of float) (1-by-3)
    '''
    global _dataset

    if component is not None:
        assert isinstance(component, Component)
        zonelist = component.zonelist
    else:
        component = Component('All')
        zonelist = list(_dataset.zones())

    if force_variables is not None:
        for v in force_variables:
            _dataset.add_variable(v)

    for zone in zonelist:
        pres, fric, mmnt, area = zonal_forces(zone, cg)
        component.force_pressure[:] = component.force_pressure[:] + pres[:]
        component.force_friction[:] = component.force_friction[:] + fric[:]
        component.area = component.area + area
        component.moment = component.moment + mmnt

    return component


def mirror():
    ''' Mirror surface on X-Y plane'''
    global _dataset
    for zone in _dataset.zones():
        y = zone.values('Y').as_numpy_array()
        y = -1*y  # eflect the y axis
        zone.values('Y')[:] = y[:]


def excluded_component(components, new_name="Other"):
    ''' Filter surface zones within a given box

    arguments:
        components - (list of Component)
    optional arguments:
        new_name - (str)
    '''
    global _dataset
    zonelist = list(_dataset.zones())
    left_out_zones = []

    for zone in zonelist:
        found = False
        for comp in components:
            if zone in comp.zonelist:
                found = True
                exit
        if not found:
            left_out_zones.append(zone)

    if len(left_out_zones) == 0:
        return None
    else:
        new_component = Component(new_name)
        new_component.zonelist[:] = left_out_zones[:]
        return new_component
