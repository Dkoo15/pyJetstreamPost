import numpy as np

''' This module holds functions for computing sbp operators and grid metrics
functions:
    _sbp - Summation-by-parts operator (internal use)
    surface_metrics - metric calculation for surface drag only
    partial_block_metrics - metric calculation for surface-adjacent
                            terms necessary for friction/shear
    jac - Jacobian calculation
'''


def _sbp(n):
    '''Compute the SBP operator coefficients and pointer arrays'''
    ja = np.zeros(6, dtype=np.int16)
    a = np.zeros(6, dtype=np.float64)
    ibeg = np.zeros(n, dtype=np.int16)
    iend = np.zeros(n, dtype=np.int16)

    ptr = 0
    ibeg[0] = ptr
    ja[ptr] = 0
    a[ptr] = -1.
    ptr = ptr + 1
    ja[ptr] = 1
    a[ptr] = 1.
    iend[0] = ptr

    ptr = ptr + 1
    ibeg[1:n-1] = ptr
    ja[ptr] = -1
    a[ptr] = -0.5
    ptr = ptr + 1
    ja[ptr] = 1
    a[ptr] = 0.5
    iend[1:n-1] = ptr

    ptr = ptr + 1
    ibeg[n-1] = ptr
    ja[ptr] = -1
    a[ptr] = -1.
    ptr = ptr + 1
    ja[ptr] = 0
    a[ptr] = 1.
    iend[n-1] = ptr

    return ja, a, ibeg, iend


def surface_metrics(imax, jmax, xyz):
    '''Compute the metric terms on the surface, in the direction normal.

    arguments:
        imax - (int)
        jmax - (int)
        xyz : (numpy array of float) (n-by-3)
    '''
    ia, sbpi, ibeg, iend = _sbp(imax)
    ja, sbpj, jbeg, jend = _sbp(jmax)
    n = len(xyz)

    dxdxi_i = np.zeros((n, 3))  # dxdxi in i direction
    ind = 0
    for j in range(jmax):
        for i in range(0, imax):
            for p in [ibeg[i], iend[i]]:
                sbp = sbpi[p]
                indp = ind + ia[p]
                for di in [0, 1, 2]:
                    dxdxi_i[ind, di] = dxdxi_i[ind, di] + sbp*xyz[indp, di]
            ind = ind + 1

    # dxdxi in j direction
    dxdxi_j = np.zeros((n, 3))  # dxdxi in j direction
    indo = 0
    for j in range(jmax):
        for p in [jbeg[j], jend[j]]:
            jp = ja[p]*imax
            sbp = sbpj[p]
            for i in range(imax):
                ind = indo + i
                indp = indo + i + jp
                for di in [0, 1, 2]:
                    dxdxi_j[ind, di] = dxdxi_j[ind, di] + sbp*xyz[indp, di]
        indo = indo + imax

    dxidx = np.zeros((n, 3))  # di = i direction
    ind = 0
    for j in range(jmax):
        for i in range(imax):
            for p in [ibeg[i], iend[i]]:
                sbp = 0.5*sbpi[p]
                indp = ind + ia[p]
                dxidx[ind, 0] = (dxidx[ind, 0] +
                                 sbp*(xyz[indp, 1]*dxdxi_j[indp, 2] -
                                      xyz[indp, 2]*dxdxi_j[indp, 1]))
                dxidx[ind, 1] = (dxidx[ind, 1] +
                                 sbp*(xyz[indp, 2]*dxdxi_j[indp, 0] -
                                      xyz[indp, 0]*dxdxi_j[indp, 2]))
                dxidx[ind, 2] = (dxidx[ind, 2] +
                                 sbp*(xyz[indp, 0]*dxdxi_j[indp, 1] -
                                      xyz[indp, 1]*dxdxi_j[indp, 0]))
            ind = ind + 1

    indo = 0  # di = j direction
    for j in range(jmax):
        for p in [jbeg[j], jend[j]]:
            jp = ja[p]*imax
            sbp = 0.5*sbpj[p]
            for i in range(imax):
                ind = indo + i
                indp = indo + i + jp
                dxidx[ind, 0] = (dxidx[ind, 0] +
                                 sbp*(xyz[indp, 2]*dxdxi_i[indp, 1] -
                                      xyz[indp, 1]*dxdxi_i[indp, 2]))
                dxidx[ind, 1] = (dxidx[ind, 1] +
                                 sbp*(xyz[indp, 0]*dxdxi_i[indp, 2] -
                                      xyz[indp, 2]*dxdxi_i[indp, 0]))
                dxidx[ind, 2] = (dxidx[ind, 2] +
                                 sbp*(xyz[indp, 1]*dxdxi_i[indp, 0] -
                                      xyz[indp, 0]*dxdxi_i[indp, 1]))
        indo = indo + imax

    return dxidx


def partial_block_metrics(xyz, imax, jmax, n, wallindx, sgn):
    ''' Compute the grid metrics we need to get friction

    arguments
        xyz - (numpy array of floats) (n-by-3)
        imax - (int)
        jmax - (int)
        n - (int) number of points
        wallindx - (int)
        sgn - (float) -1 or +1
    '''

    dxdxi = np.zeros((n, 3, 3))
    for kplane in [wallindx, wallindx-sgn*imax*jmax]:  # need 2 levels of dxdxi
        # dxdxi in i direction
        ind = kplane
        ia, sbpi, ibeg, iend = _sbp(imax)  # define sbp operators in i and j
        for j in range(jmax):
            for i in range(imax):
                for p in [ibeg[i], iend[i]]:
                    sbp = sbpi[p]
                    indp = ind + ia[p]
                    dxdxi[ind, 0, :] = dxdxi[ind, 0, :] + sbp*xyz[indp, :]
                ind = ind + 1

        # dxdxi in j direction
        indo = kplane
        ja, sbpj, jbeg, jend = _sbp(jmax)
        for j in range(0, jmax):
            for p in [jbeg[j], jend[j]]:
                jp = ja[p]*imax
                sbp = sbpj[p]
                for i in range(0, imax):
                    ind = indo + i
                    indp = indo + i + jp
                    dxdxi[ind, 1, :] = dxdxi[ind, 1, :] + sbp*xyz[indp, :]
            indo = indo + imax

    # dxdxi in the k direction
    ind = wallindx
    sbpk = [sgn, -sgn]
    ka = [0, -sgn*imax*jmax]
    for j in range(jmax):
        for i in range(imax):
            for p in [0, 1]:
                indp = ind + ka[p]
                sbp = sbpk[p]
                dxdxi[ind, 2, :] = dxdxi[ind, 2, :] + sbp*xyz[indp, :]
            ind = ind + 1

    # di in i direction
    dxidx = np.zeros((n, 3, 3))
    ind = wallindx
    for j in range(jmax):
        for i in range(imax):
            for p in [ibeg[i], iend[i]]:
                sbp = 0.5*sbpi[p]
                indp = ind + ia[p]
                # di = 0, it1 = 1
                dxidx[ind, 1, 0] = (dxidx[ind, 1, 0] +
                                    sbp*(xyz[indp, 2]*dxdxi[indp, 2, 1] -
                                         xyz[indp, 1]*dxdxi[indp, 2, 2]))
                dxidx[ind, 1, 1] = (dxidx[ind, 1, 1] +
                                    sbp*(xyz[indp, 0]*dxdxi[indp, 2, 2] -
                                         xyz[indp, 2]*dxdxi[indp, 2, 0]))
                dxidx[ind, 1, 2] = (dxidx[ind, 1, 2] +
                                    sbp*(xyz[indp, 1]*dxdxi[indp, 2, 0] -
                                         xyz[indp, 0]*dxdxi[indp, 2, 1]))
                # di = 0, it2 = 2
                dxidx[ind, 2, 0] = (dxidx[ind, 2, 0] +
                                    sbp*(xyz[indp, 1]*dxdxi[indp, 1, 2] -
                                         xyz[indp, 2]*dxdxi[indp, 1, 1]))
                dxidx[ind, 2, 1] = (dxidx[ind, 2, 1] +
                                    sbp*(xyz[indp, 2]*dxdxi[indp, 1, 0] -
                                         xyz[indp, 0]*dxdxi[indp, 1, 2]))
                dxidx[ind, 2, 2] = (dxidx[ind, 2, 2] +
                                    sbp*(xyz[indp, 0]*dxdxi[indp, 1, 1] -
                                         xyz[indp, 1]*dxdxi[indp, 1, 0]))
            ind = ind + 1

    # di in j direction
    indo = wallindx
    for j in range(jmax):
        for p in [jbeg[j], jend[j]]:
            jp = ja[p]*imax
            sbp = 0.5*sbpj[p]
            for i in range(imax):
                ind = indo + i
                indp = indo + i + jp
                # di = 1, it1 = 2
                dxidx[ind, 2, 0] = (dxidx[ind, 2, 0] +
                                    sbp*(xyz[indp, 2]*dxdxi[indp, 0, 1] -
                                         xyz[indp, 1]*dxdxi[indp, 0, 2]))
                dxidx[ind, 2, 1] = (dxidx[ind, 2, 1] +
                                    sbp*(xyz[indp, 0]*dxdxi[indp, 0, 2] -
                                         xyz[indp, 2]*dxdxi[indp, 0, 0]))
                dxidx[ind, 2, 2] = (dxidx[ind, 2, 2] +
                                    sbp*(xyz[indp, 1]*dxdxi[indp, 0, 0] -
                                         xyz[indp, 0]*dxdxi[indp, 0, 1]))
                # di = 1, it2 = 0
                dxidx[ind, 0, 0] = (dxidx[ind, 0, 0] +
                                    sbp*(xyz[indp, 1]*dxdxi[indp, 2, 2] -
                                         xyz[indp, 2]*dxdxi[indp, 2, 1]))
                dxidx[ind, 0, 1] = (dxidx[ind, 0, 1] +
                                    sbp*(xyz[indp, 2]*dxdxi[indp, 2, 0] -
                                         xyz[indp, 0]*dxdxi[indp, 2, 2]))
                dxidx[ind, 0, 2] = (dxidx[ind, 0, 2] +
                                    sbp*(xyz[indp, 0]*dxdxi[indp, 2, 1] -
                                         xyz[indp, 1]*dxdxi[indp, 2, 0]))
        indo = indo + imax

    # di in k direction
    ind = wallindx
    for j in range(jmax):
        for i in range(imax):
            for p in [0, 1]:
                indp = ind + ka[p]
                sbp = 0.5*sbpk[p]
                # di = 2, it1 = 0
                dxidx[ind, 0, 0] = (dxidx[ind, 0, 0] +
                                    sbp*(xyz[indp, 2]*dxdxi[indp, 1, 1] -
                                         xyz[indp, 1]*dxdxi[indp, 1, 2]))
                dxidx[ind, 0, 1] = (dxidx[ind, 0, 1] +
                                    sbp*(xyz[indp, 0]*dxdxi[indp, 1, 2] -
                                         xyz[indp, 2]*dxdxi[indp, 1, 0]))
                dxidx[ind, 0, 2] = (dxidx[ind, 0, 2] +
                                    sbp*(xyz[indp, 1]*dxdxi[indp, 1, 0] -
                                         xyz[indp, 0]*dxdxi[indp, 1, 1]))
                # di = 2, it2 = 1
                dxidx[ind, 1, 0] = (dxidx[ind, 1, 0] +
                                    sbp*(xyz[indp, 1]*dxdxi[indp, 0, 2] -
                                         xyz[indp, 2]*dxdxi[indp, 0, 1]))
                dxidx[ind, 1, 1] = (dxidx[ind, 1, 1] +
                                    sbp*(xyz[indp, 2]*dxdxi[indp, 0, 0] -
                                         xyz[indp, 0]*dxdxi[indp, 0, 2]))
                dxidx[ind, 1, 2] = (dxidx[ind, 1, 2] +
                                    sbp*(xyz[indp, 0]*dxdxi[indp, 0, 1] -
                                         xyz[indp, 1]*dxdxi[indp, 0, 0]))

            ind = ind + 1

    return dxidx


def jac(xi):
    '''Compute the jacobian from the metric terms'''
    # xi is 3x3
    j = (xi[0, 0]*xi[1, 1]*xi[2, 2] +
         xi[0, 1]*xi[1, 2]*xi[2, 0] +
         xi[0, 2]*xi[1, 0]*xi[2, 1] -
         xi[0, 0]*xi[1, 2]*xi[2, 1] -
         xi[0, 1]*xi[1, 0]*xi[2, 2] -
         xi[0, 2]*xi[1, 1]*xi[2, 0])
    j = np.sqrt(1/j)
    return j
