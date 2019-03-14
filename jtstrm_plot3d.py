import tecplot
import numpy as np
import block_utils

''' This module is for loading Jetstream's plot3d files

classes:
    Boundary - Container for boundary data

functions:
    find_boundaries - Parse grid.con file for surfaces
    extract_surface - Subsamples volume zones for surface data
    get_aero_surfaces - Overall wrapper function for surface extraction
    compute_cpcf - Compute CP/CF for a plot3d block
'''
_dataset = None


class Boundary(object):
    '''Class for data and calculations associated with a boundary'''

    def __init__(self, blk, face):
        '''Constructor for Boundary.

        arguments:
            block - (int)
            face - (int)
        '''
        self.block = blk
        self.face = face
        self.data = None
        self.nvar = 0
        self.imax = 0
        self.jmax = 0

    def assign_data(self, imax, jmax, data, num_variables):
        '''Assign coordinates and other data to Boundary

        arguments:
            imax - (int)
            jmax - (int)
            data - (numpy array of float) (n-by-m)
            num_variables - (int)
        '''
        self.imax = imax
        self.jmax = jmax
        self.data = data
        self.nvar = num_variables


def load_solution(g_file, q_file):
    '''Load the Plot3D data into _dataset module variable

    arguments:
        g_file - (str) .g file path
        q_file - (str) .q file path
    '''
    global _dataset
    try:
        _dataset = tecplot.data.load_plot3d(g_file, q_file,
                                            append=False,
                                            include_boundaries=False)
    except (tecplot.exception.TecplotSystemError,
            tecplot.exception.TecplotTypeError):
        print('-> Error loading file %s, or %s' % (g_file, q_file))
        raise SystemExit

    tecplot.active_frame().plot_type = tecplot.constant.PlotType.Cartesian3D
    print('Loaded %s and %s' % (g_file, q_file))


def find_boundaries(btype='1', con_file='grid.con'):
    '''Read Jetstream's connectivity file, establishing boundaries

    optional arguments:
        btype  - (str) actually an integer, the surface code
        con_file - (str) grid.con file path
    '''
    f = open(con_file, 'r')
    for i in range(9):
        f.readline()  # data starts in line 9
    connectivity = f.read().splitlines()
    f.close()

    boundaries = []
    for face in connectivity:
        numbers = face.split()
        if(len(numbers) < 10):
            continue  # not a boundary
        if(numbers[1] != btype):
            continue  # not a surface we want
        bdy = Boundary(int(numbers[2])-1, int(numbers[3]))
        boundaries.append(bdy)

    return boundaries


def extract_surface(boundary, variables):
    ''' Take the values from the boundary of a volume zone and place in Boundary

    arguments:
        boundary : (Boundary)
        variables : (list of str)
    '''
    global _dataset
    assert(isinstance(boundary, Boundary))

    volzone = _dataset.zone(boundary.block)
    bface = boundary.face
    imax, jmax, kmax = volzone.dimensions
    nvar = len(variables)

    if bface == 1:  # not correct for metrics, need transpose!
        m = jmax*kmax
        surfvar = np.zeros((m, nvar))
        for n, v in enumerate(variables):
            var = volzone.values(v).as_numpy_array()
            surfvar[:, n] = var[0::imax]
        imax = jmax
        jmax = kmax

    elif bface == 2:
        m = jmax*kmax
        surfvar = np.zeros((m, nvar))
        for n, v in enumerate(variables):
            var = volzone.values(v).as_numpy_array()
            surfvar[:, n] = var[imax-1::imax]
        imax = jmax
        jmax = kmax

    elif bface == 3:  # not correct for metrics, need transpose!
        m = imax*kmax
        surfvar = np.zeros((m, nvar))
        for n, v in enumerate(variables):
            var = volzone.values(v).as_numpy_array()
            for k in range(kmax):
                surfvar[k*imax:(k+1)*imax, n] = var[k*jmax*imax:(k*jmax+1)*imax]
        jmax = kmax

    elif bface == 4:
        m = imax*kmax
        surfvar = np.zeros((m, nvar))
        for n, v in enumerate(variables):
            var = volzone.values(v).as_numpy_array()
            for k in range(kmax):
                surfvar[k*imax:(k+1)*imax, n] = var[((k+1)*jmax-1)*imax:(k+1)*jmax*imax]
        jmax = kmax

    elif bface == 5:  # Transpose data for metric calculation
        m = imax*jmax
        surfvar = np.zeros((m, nvar), dtype=np.float64)
        for n, v in enumerate(variables):
            var = volzone.values(v).as_numpy_array()
            for i in range(imax):
                surfvar[i*jmax:(i+1)*jmax, n] = var[i:m:imax]

        tmp = jmax
        jmax = imax
        imax = tmp  # swap coordinates

    elif bface==6:
        m = imax*jmax
        surfvar = np.zeros((m, nvar), dtype=np.float64)
        for n, v in enumerate(variables):
            var = volzone.values(v).as_numpy_array()
            surfvar[:, n] = var[-m:]
    else:
        print("Error, incorrect bface!")
        raise SystemExit

    boundary.assign_data(imax, jmax, surfvar, nvar)


def get_aero_surfaces(coordinates=('X', 'Y', 'Z'),
                      coefficients=('Cp', 'Cf', 'Cfx', 'Cfy', 'Cfz'),
                      con_file='grid.con'):
    '''Converts the current _dataset into a surface data file.

    optional arguments:
        coordinates : (tuple of str)
        coefficients : (tuple of str)
        con_file : (str) grid.con file path
    '''
    global _dataset
    if _dataset is None:
        print("Error, dataset not yet loaded!")
        raise SystemExit

    variables = coordinates + coefficients

    mach = float(_dataset.aux_data['Common.ReferenceMachNumber'])
    alpha = float(_dataset.aux_data['Common.AngleOfAttack'])
    re = float(_dataset.aux_data['Common.ReynoldsNumber'])

    boundaries = find_boundaries(con_file=con_file)

    for c in coefficients:
        _dataset.add_variable(c)

    for b in boundaries:
        compute_cpcf(b, re, mach)
        extract_surface(b, variables)

    return boundaries


def compute_cpcf(boundary, re, mach):
    '''Compute friction and pressure coefficients for the block of this boundary

    arguments:
        boundary - Boundary
        re - (float)
        mach - (float)
    '''
    global _dataset
    assert(isinstance(boundary, Boundary))

    volzone = _dataset.zone(boundary.block)
    bface = boundary.face

    imax, jmax, kmax = volzone.dimensions
    n = volzone.num_points
    offset = imax*jmax
    cf = np.zeros((n, 3))
    cf_mag = np.zeros(n)
    cp = np.zeros(n)

    if bface == 5:
        wallidx = 0
        sgn = -1
    else:
        wallidx = n - offset
        sgn = 1

    xyz = np.zeros((n, 3))
    xyz[:, 0] = volzone.values('X').as_numpy_array()
    xyz[:, 1] = volzone.values('Y').as_numpy_array()
    xyz[:, 2] = volzone.values('Z').as_numpy_array()

    xi = block_utils.partial_block_metrics(xyz, imax, jmax, n, wallidx, sgn)

    # mu = mach/Re #farfield viscosity
    ST = 198.6/460.0  # Sutherland's Law value
    cfac = -2/(mach*mach)  # factor

    rho = volzone.values(3).as_numpy_array()
    u = volzone.values(4).as_numpy_array()/rho
    v = volzone.values(5).as_numpy_array()/rho
    w = volzone.values(6).as_numpy_array()/rho
    e = volzone.values(7).as_numpy_array()

    v1 = np.zeros(3)
    v2 = np.zeros(3)
    v3 = np.zeros(3)

    ind = wallidx
    for j in range(jmax):
        for i in range(imax):
            # differences in i direction
            if i == 0:
                u_xi = -1.5*u[ind] + 2*u[ind+1] - 0.5*u[ind+2]
                v_xi = -1.5*v[ind] + 2*v[ind+1] - 0.5*v[ind+2]
                w_xi = -1.5*w[ind] + 2*w[ind+1] - 0.5*w[ind+2]
            elif i == imax-1:
                u_xi = 1.5*u[ind] - 2*u[ind-1] - 0.5*u[ind-2]
                v_xi = 1.5*v[ind] - 2*v[ind-1] - 0.5*v[ind-2]
                w_xi = 1.5*w[ind] - 2*w[ind-1] - 0.5*w[ind-2]
            else:
                u_xi = 0.5*(u[ind+1] - u[ind-1])
                v_xi = 0.5*(v[ind+1] - v[ind-1])
                w_xi = 0.5*(w[ind+1] - w[ind-1])
            # differences in j direction
            if j == 0:
                u_eta = -1.5*u[ind] + 2*u[ind+imax] - 0.5*u[ind+imax+imax]
                v_eta = -1.5*v[ind] + 2*v[ind+imax] - 0.5*v[ind+imax+imax]
                w_eta = -1.5*w[ind] + 2*w[ind+imax] - 0.5*w[ind+imax+imax]
            elif j == jmax-1:
                u_eta = 1.5*u[ind] - 2*u[ind-imax] - 0.5*u[ind-imax-imax]
                v_eta = 1.5*v[ind] - 2*v[ind-imax] - 0.5*v[ind-imax-imax]
                w_eta = 1.5*w[ind] - 2*w[ind-imax] - 0.5*w[ind-imax-imax]
            else:
                u_eta = 0.5*(u[ind+imax] - u[ind-imax])
                v_eta = 0.5*(v[ind+imax] - v[ind-imax])
                w_eta = 0.5*(w[ind+imax] - w[ind-imax])
            # differences in k direction
            if bface == 5:  # k==0
                u_zeta = -1.5*u[ind] + 2*u[ind+offset] - 0.5*u[ind+2*offset]
                v_zeta = -1.5*v[ind] + 2*v[ind+offset] - 0.5*v[ind+2*offset]
                w_zeta = -1.5*w[ind] + 2*w[ind+offset] - 0.5*w[ind+2*offset]
            else:  # k==kmax
                u_zeta = 1.5*u[ind] - 2*u[ind-offset] + 0.5*u[ind-2*offset]
                v_zeta = 1.5*v[ind] - 2*v[ind-offset] + 0.5*v[ind-2*offset]
                w_zeta = 1.5*w[ind] - 2*w[ind-offset] + 0.5*w[ind-2*offset]

            J = block_utils.jac(xi[ind, :, :])
            xi_x = xi[ind, 0, 0]
            xi_y = xi[ind, 0, 1]
            xi_z = xi[ind, 0, 2]
            eta_x = xi[ind, 1, 0]
            eta_y = xi[ind, 1, 1]
            eta_z = xi[ind, 1, 2]
            zeta_x = xi[ind, 2, 0]
            zeta_y = xi[ind, 2, 1]
            zeta_z = xi[ind, 2, 2]

            p = 0.4*(e[ind] - 0.5*rho[ind]*(u[ind]**2 + v[ind]**2 + w[ind]**2))
            a2 = 1.4*p/rho[ind]

            mu = np.sqrt(a2)*a2*(ST+1)/(a2+ST) * (mach/re)  # local mu

            tauxy = mu * (xi_y*u_xi + eta_y*u_eta + zeta_y*u_zeta -
                          (xi_x*v_xi + eta_x*v_eta + zeta_x*v_zeta)) * J
            tauxz = -mu * (xi_z*u_xi + eta_z*u_eta + zeta_z*u_zeta -
                           (xi_x*w_xi + eta_x*w_eta + zeta_x*w_zeta)) * J
            tauyz = mu * (xi_z*v_xi + eta_z*v_eta + zeta_z*v_zeta -
                          (xi_y*w_xi + eta_y*w_eta + zeta_y*w_zeta)) * J

            # di=2, it1=0, it2=1
            nrm = xi[ind, 2, :]
            mag = np.sqrt(nrm[0]*nrm[0] + nrm[1]*nrm[1] + nrm[2]*nrm[2])
            v1[2] = nrm[0]
            v1[0] = -nrm[2]
            v1[1] = 0.
            v2[2] = -nrm[1]
            v2[0] = 0.
            v2[1] = nrm[2]
            v3[2] = 0.
            v3[0] = nrm[1]
            v3[1] = -nrm[0]

            cf[ind, :] = (tauxz*v1[:] + tauyz*v2[:] + tauxy*v3[:])*cfac/mag*sgn
            cf_mag[ind] = np.sqrt(cf[ind, 0]**2 + cf[ind, 1]**2 +
                                  cf[ind, 2]**2)
            cp[ind] = 2*(1.4*p-1) / (1.4*mach**2)
            ind = ind + 1

    if 'Cp' in _dataset.variable_names:
        volzone.values('Cp')[:] = cp[:]
    if 'Cf' in _dataset.variable_names:
        volzone.values('Cf')[:] = cf_mag[:]
    if 'Cfx' in _dataset.variable_names:
        volzone.values('Cfx')[:] = cf[:, 0]
    if 'Cfy' in _dataset.variable_names:
        volzone.values('Cfy')[:] = cf[:, 1]
    if 'Cfz' in _dataset.variable_names:
        volzone.values('Cfz')[:] = cf[:, 2]

    return
