import numpy as np
import jtstrm_surface

'''Module for extracting sections from Jetstream surface file

classes:
    SurfaceSection - Container for section data

functions:
    init_surface - Wrapper for functions in jtstrm_surface
    cpcf_sections - Loops over sections, extracting CP profiles
    geometry_sections - Extracts spanwise profile for twist / thickness
    force_sections - Extracts spanwise lift distribution
'''


class SurfaceSection(object):
    '''Class for data and calculations associated with a surface section'''

    def __init__(self, surf_label, y,
                 gamma=0.0, chord=1.0, twist=0.0):
        '''Constructor for section.

        arguments:
            surf_label - (str) name of the dataset
            y - (float) spanwise location
        optional arguments:
            gamma - (float) dihedral angle
            chord - (float) chord length
            twist - (float) aerodynamic twist
        '''
        self.label = surf_label
        self.location = y
        self.dihedral = gamma
        self.chord = chord
        self.twist = twist
        self.alpha = 0

        self.data = None
        self.nvar = 0
        self.nodemap = None

    def assign_data(self, data, nodemap=None):
        '''Assign coordinate data and nodemaps

        arguments:
            data - (numpy array) coordinates and data
        optional arguments:
            nodemap - (list) nodemap from Tecplot
        '''
        self.data = data
        self.nodemap = nodemap

    def thickness(self, loc):
        '''Compute thickness at the specified x/c location.

        Arguments:
            loc - (float) x/c where to compute thickness
        '''
        xc = self.data[:, 0]
        zc = self.data[:, 1]
        xdist = np.abs(xc - loc)
        pts = xdist.argsort()[:4]
        zloc = zc[pts]
        thick = zloc.max() - zloc.min()
        return thick

    def max_thickness(self):
        '''Compute maximum thickness of this section.'''
        locs = np.linspace(0.1, 0.6, 21)
        thicknesses = np.zeros(21)
        for i, x in enumerate(locs):
            thicknesses[i] = self.thickness(x)

        return thicknesses.max()

    def section_lift(self):
        '''Compute the lift force integrated over the section'''
        lift = 0.0
        lvec = np.array([-np.sin(self.alpha),
                         np.sin(self.dihedral)*np.cos(self.alpha),
                         np.cos(self.dihedral)*np.cos(self.alpha)])

        if self.nodemap is None or self.data.shape[1] < 5:
            return lift

        x = self.data[:, 0]*self.chord
        z = self.data[:, 1]*self.chord
        for arc in self.nodemap:
            j = arc[0]
            k = arc[1]
            ds = np.sqrt((x[j] - x[k])**2 + (z[j] - z[k])**2)
            frc = 0.5*self.data[j, 2:5] + 0.5*self.data[k, 2:5]
            lift = lift + np.dot(frc, lvec)*ds

        return lift


def init_surface(filename, box=None):
    '''Wrapper to initialize jtstrm_surface module'''
    jtstrm_surface.load(filename)
    jtstrm_surface.filter(box)


def cpcf_sections(label, stations, coeffs, section_list,
                  dihedral, winglet_y=1.1, untwist=True):
    '''Loop over the stations and extract slices for CP/CF plots.
    init_surface() must be run.

    arguments:
        label - (str) dataset name
        stations - (list of float) spanwise stations
        coeffs - (list of str) coefficients to extract
        section_list - (list of SurfaceSection)
    optional arguments
        winglet_y - (float) winglet location
    '''
    tip, _, _, = jtstrm_surface.span()
    for stn in stations:
        y = stn*tip

        if stn > winglet_y:
            data, nmap, twst, gamm, chrd = jtstrm_surface.dihedral_section(
                y, variables=coeffs, get_nodemap=True, untwist=untwist)
        elif dihedral > 1e-4:
            data, nmap, twst, gamm, chrd = jtstrm_surface.dihedral_section(
                y, dihedral=dihedral, variables=coeffs, get_nodemap=True,
                untwist=untwist)
        else:
            data, nmap, twst, gamm, chrd = jtstrm_surface.cross_section(
                (0, y, 0), variables=coeffs, get_nodemap=True, untwist=untwist)

        extracted = SurfaceSection(label, stn, gamm, chrd, twst)
        extracted.assign_data(data, nmap)
        section_list.append(extracted)


def geometry_sections(label, section_list,
                      num_sections=21, winglet_y=1.1):
    '''Loop over the stations and extract slices for twist / thickness
    init_surface() must be run.

    arguments:
        label - (str) dataset name
        section_list - (list of SurfaceSection)
    optional arguments:
        num_sections - (int)
        winglet_y  - (float) winglet location
    '''
    tip, root, _, = jtstrm_surface.span()
    stations = np.linspace(root, tip, num_sections)

    for y in stations:
        if y/tip > winglet_y:
            airfoil, _, twst, gamm, chrd = jtstrm_surface.dihedral_section(y)
        else:
            airfoil, _, twst, gamm, chrd = jtstrm_surface.cross_section((0, y,
                                                                         0))

        extracted = SurfaceSection(label, y/tip, gamm, chrd, twst)
        extracted.assign_data(airfoil)
        section_list.append(extracted)


def force_sections(label, section_list, cg,
                   dihedral=0.0, winglet_y=1.1,
                   component=None, y_locations=None):
    '''Loop over the stations and extract slices for twist / thickness
    init_surface() must be run

    arguments:
        label - (str) dataset name
        section_list - (list of SurfaceSection)
        cg - (list of float) centre of gravity
    optional arguments:
        num_sections - (int)
        winglet_y  - (float) winglet location
    '''
    from component import Component

    forces = ('Fx', 'Fy', 'Fz')
    component = jtstrm_surface.surface_force(component=component,
                                             force_variables=forces, cg=cg)
    tip, root, _, = jtstrm_surface.span()

    if y_locations is None:
        stations = np.linspace(root, tip, 51)
    else:
        stations = np.copy(y_locations)

    aoa = jtstrm_surface.angle_of_attack()*np.pi/180
    lvec = np.array([-np.sin(aoa), 0, np.cos(aoa)])

    for y in stations:
        if y/tip > winglet_y:
            data, nmap, twst, gamm, chrd = jtstrm_surface.dihedral_section(
                y, variables=forces, get_nodemap=True, untwist=False)
        elif dihedral > 1e-4:
            data, nmap, twst, gamm, chrd = jtstrm_surface.dihedral_section(
                y, dihedral=dihedral, variables=forces, get_nodemap=True,
                untwist=False)
        else:
            data, nmap, twst, gamm, chrd = jtstrm_surface.cross_section(
                (0, y, 0), variables=forces, get_nodemap=True,
                untwist=False)

        extracted = SurfaceSection(label, y, gamm, chrd, twst)
        extracted.assign_data(data, nmap)
        extracted.alpha = aoa
        section_list.append(extracted)

    lift = component.total_force(lvec)

    return lift, stations
