import numpy as np
import scipy.spatial as spspatial
import scipy.sparse as spsparse

def simpvol(p, t):
    """Signed volumes of the simplex elements in the mesh."""
    dim = p.shape[1]
    if dim == 1:
        d01 = p[t[:,1]]-p[t[:,0]]
        return d01
    elif dim == 2:
        d01 = p[t[:,1]]-p[t[:,0]]
        d02 = p[t[:,2]]-p[t[:,0]]
        return (d01[:,0]*d02[:,1]-d01[:,1]*d02[:,0])/2
    else:
        raise NotImplementedError
        
def simpqual(p, t):
    """Simplex quality.

    Usage
    -----
    q = simpqual(p, t)

    Parameters
    ----------
    p : array, shape (np, dim)
        nodes
    t : array, shape (nt, dim+1)
        triangulation

    Returns
    -------
    q : array, shape (nt, )
        qualities
    """
    assert (p.ndim == 2
            and t.ndim == 2
            and p.shape[1]+1 == t.shape[1])

    dim = p.shape[1]
    if dim == 1:
        return np.ones((1, nt))

    elif dim == 2:
        length = lambda p1: np.sqrt((p1**2).sum(1))
        a = length(p[t[:,1]]-p[t[:,0]])
        b = length(p[t[:,2]]-p[t[:,0]])
        c = length(p[t[:,2]]-p[t[:,1]])
        r = 0.5*np.sqrt((b+c-a)*(c+a-b)*(a+b-c)/(a+b+c))
        R = a*b*c/np.sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c))
        return 2*r/R

    else:
        raise NotImplementedError
        


def fixmesh(p, t, ptol=2e-13):
    """Remove duplicated/unused nodes and fix element orientation.

    Parameters
    ----------
    p : array, shape (np, dim)
    t : array, shape (nt, nf)

    Usage
    -----
    p, t = fixmesh(p, t, ptol)
    """
    snap = (p.max(0)-p.min(0)).max()*ptol
    _, ix, jx = unique_rows(np.round(p/snap)*snap, True, True)

    p = p[ix]
    t = jx[t]

    flip = simpvol(p,t)<0
    t[flip, :2] = t[flip, 1::-1]

    return p, t


def dense(I, J, S, shape=None, dtype=None):
    """
    Similar to MATLAB's SPARSE(I, J, S, ...), but instead returning a
    dense array.

    Usage
    -----
    >>> shape = (m, n)
    >>> A = dense(I, J, S, shape, dtype)
    """

    # Advanced usage: allow J and S to be scalars.
    if np.isscalar(J):
        x = J
        J = np.empty(I.shape, dtype=int)
        J.fill(x)
    if np.isscalar(S):
        x = S
        S = np.empty(I.shape)
        S.fill(x)

    # Turn these into 1-d arrays for processing.
    S = S.flat; I = I.flat; J = J.flat
    return spsparse.coo_matrix((S, (I, J)), shape, dtype).toarray()
    


def setdiff_rows(A, B, return_index=False):
    """
    Similar to MATLAB's setdiff(A, B, 'rows'), this returns C, I
    where C are the row of A that are not in B and I satisfies
    C = A[I,:].

    Returns I if return_index is True.
    """
    A = np.require(A, requirements='C')
    B = np.require(B, requirements='C')

    assert A.ndim == 2, "array must be 2-dim'l"
    assert B.ndim == 2, "array must be 2-dim'l"
    assert A.shape[1] == B.shape[1], \
           "arrays must have the same number of columns"
    assert A.dtype == B.dtype, \
           "arrays must have the same data type"

    # NumPy provides setdiff1d, which operates only on one dimensional
    # arrays. To make the array one-dimensional, we interpret each row
    # as being a string of characters of the appropriate length.
    orig_dtype = A.dtype
    ncolumns = A.shape[1]
    dtype = np.dtype((np.character, orig_dtype.itemsize*ncolumns))
    C = np.setdiff1d(A.view(dtype), B.view(dtype)) \
        .view(A.dtype) \
        .reshape((-1, ncolumns), order='C')
    if return_index:
        raise NotImplementedError
    else:
        return C

def unique_rows(A, return_index=False, return_inverse=False):
    """
    Similar to MATLAB's unique(A, 'rows'), this returns B, I, J
    where B is the unique rows of A and I and J satisfy
    A = B[J,:] and B = A[I,:]

    Returns I if return_index is True
    Returns J if return_inverse is True
    """
    A = np.require(A, requirements='C')
    assert A.ndim == 2, "array must be 2-dim'l"

    orig_dtype = A.dtype
    ncolumns = A.shape[1]
    dtype = np.dtype((np.character, orig_dtype.itemsize*ncolumns))
    B, I, J = np.unique(A.view(dtype),
                        return_index=True,
                        return_inverse=True)

    B = B.view(orig_dtype).reshape((-1, ncolumns), order='C')

    # There must be a better way to do this:
    if (return_index):
        if (return_inverse):
            return B, I, J
        else:
            return B, I
    else:
        if (return_inverse):
            return B, J
        else:
            return B


def distmesh2d(fd, fh, h0, bbox, pfix=None, fig='gcf'):
    """
    distmesh2d: 2-D Mesh Generator using Distance Functions.

    Usage
    -----
    >>> p, t = distmesh2d(fd, fh, h0, bbox, pfix)

    Parameters
    ----------
    fd:        Distance function d(x,y)
    fh:        Scaled edge length function h(x,y)
    h0:        Initial edge length
    bbox:      Bounding box, (xmin, ymin, xmax, ymax)
    pfix:      Fixed node positions, shape (nfix, 2)
    fig:       Figure to use for plotting, or None to disable plotting.

    Returns
    -------
    p:         Node positions (Nx2)
    t:         Triangle indices (NTx3)

    Example: (Uniform Mesh on Unit Circle)
    >>> fd = lambda p: sqrt((p**2).sum(1))-1.0
    >>> p, t = distmesh2d(fd, huniform, 2, (-1,-1,1,1))

    Example: (Rectangle with circular hole, refined at circle boundary)
    >>> fd = lambda p: ddiff(drectangle(p,-1,1,-1,1), dcircle(p,0,0,0.5))
    >>> fh = lambda p: 0.05+0.3*dcircle(p,0,0,0.5)
    >>> p, t = distmesh2d(fd, fh, 0.05, (-1,-1,1,1),
                          [(-1,-1), (-1,1), (1,-1), (1,1)])

    Example: (Polygon)
    >>> pv=[(-0.4, -0.5), (0.4, -0.2), (0.4, -0.7), (1.5, -0.4), (0.9, 0.1),
            (1.6, 0.8), (0.5, 0.5), (0.2, 1.0), (0.1, 0.4), (-0.7, 0.7),
            (-0.4, -0.5)]
    >>> fd = lambda p: dpoly(p, pv)
    >>> p, t = distmesh2d(fd, huniform, 0.1, (-1,-1, 2,1), pv)

    Example: (Ellipse)
    >>> fd = lambda p: p[:,0]**2/2**2 + p[:,1]**2/1**2 - 1
    >>> p, t = dm.distmesh2d(fd, dm.huniform, 0.2, (-2,-1, 2,1))

    Example: (Square, with size function point and line sources)
    >>> fd = lambda p: dm.drectangle(p,0,1,0,1)
    >>> fh = lambda p: np.minimum(np.minimum(
            0.01+0.3*abs(dm.dcircle(p,0,0,0)),
            0.025+0.3*abs(dm.dpoly(p,[(0.3,0.7),(0.7,0.5)]))), 0.15)
    >>> p, t = dm.distmesh2d(fd, fh, 0.01, (0,0,1,1), [(0,0),(1,0),(0,1),(1,1)])

    Example: (NACA0012 airfoil)
    >>> hlead=0.01; htrail=0.04; hmax=2; circx=2; circr=4
    >>> a=.12/.2*np.array([0.2969,-0.1260,-0.3516,0.2843,-0.1036])
    >>> a0=a[0]; a1=np.hstack((a[5:0:-1], 0.0))
    >>> fd = lambda p: dm.ddiff(
        dm.dcircle(p,circx,0,circr),
        (abs(p[:,1])-np.polyval(a1, p[:,0]))**2-a0**2*p[:,0])
    >>> fh = lambda p: np.minimum(np.minimum(
            hlead+0.3*dm.dcircle(p,0,0,0),
            htrail+0.3*dm.dcircle(p,1,0,0)), hmax)

    >>> fixx = 1.0-htrail*np.cumsum(1.3**np.arange(5))
    >>> fixy = a0*np.sqrt(fixx)+np.polyval(a1, fixx)
    >>> fix = np.vstack((
            np.array([(circx-circr,0),(circx+circr,0),
                      (circx,-circr),(circx,circr),
                      (0,0),(1,0)]),
            np.vstack((fixx, fixy)).T,
            np.vstack((fixx, -fixy)).T))
    >>> box = (circx-circr,-circr, circx+circr,circr)
    >>> h0 = min(hlead, htrail, hmax)
    >>> p, t = dm.distmesh2d(fd, fh, h0, box, fix)
    """

    dptol=.001; ttol=.1; Fscale=1.2; deltat=.2; geps=.001*h0;
    deps=np.sqrt(np.finfo(np.double).eps)*h0;
    densityctrlfreq=30;

    # Extract bounding box
    xmin, ymin, xmax, ymax = bbox
    if pfix is not None:
        pfix = np.array(pfix, dtype='d')

    # 1. Create initial distribution in bounding box (equilateral triangles)
    x, y = np.mgrid[xmin:(xmax+h0):h0,
                    ymin:(ymax+h0*np.sqrt(3)/2):h0*np.sqrt(3)/2]
    x[:, 1::2] += h0/2                               # Shift even rows
    p = np.vstack((x.flat, y.flat)).T                # List of node coordinates

    # 2. Remove points outside the region, apply the rejection method
    p = p[fd(p)<geps]                                # Keep only d<0 points
    r0 = 1/fh(p)**2                                  # Probability to keep point
    p = p[np.random.random(p.shape[0])<r0/r0.max()]  # Rejection method
    if pfix is not None:
        p = setdiff_rows(p, pfix)                 # Remove duplicated nodes
        pfix = unique_rows(pfix); nfix = pfix.shape[0]
        p = np.vstack((pfix, p))                     # Prepend fix points
    else:
        nfix = 0
    N = p.shape[0]                                   # Number of points N

    count = 0
    pold = float('inf')                              # For first iteration

    while True:
        count += 1

        # 3. Retriangulation by the Delaunay algorithm
        dist = lambda p1, p2: np.sqrt(((p1-p2)**2).sum(1))
        if (dist(p, pold)/h0).max() > ttol:          # Any large movement?
            pold = p.copy()                          # Save current positions
            t = spspatial.Delaunay(p).vertices       # List of triangles
            pmid = p[t].sum(1)/3                     # Compute centroids
            t = t[fd(pmid) < -geps]                  # Keep interior triangles
            # 4. Describe each bar by a unique pair of nodes
            bars = np.vstack((t[:, [0,1]],
                              t[:, [1,2]],
                              t[:, [2,0]]))          # Interior bars duplicated
            bars.sort(axis=1)
            bars = unique_rows(bars)              # Bars as node pairs
            # 5. Graphical output of the current mesh

        # 6. Move mesh points based on bar lengths L and forces F
        barvec = p[bars[:,0]] - p[bars[:,1]]         # List of bar vectors
        L = np.sqrt((barvec**2).sum(1))              # L = Bar lengths
        hbars = fh(p[bars].sum(1)/2)
        L0 = (hbars*Fscale
              *np.sqrt((L**2).sum()/(hbars**2).sum()))  # L0 = Desired lengths

        # Density control - remove points that are too close
        if (count % densityctrlfreq) == 0 and (L0 > 2*L).any():
            ixdel = np.setdiff1d(bars[L0 > 2*L].reshape(-1), np.arange(nfix))
            p = p[np.setdiff1d(np.arange(N), ixdel)]
            N = p.shape[0]; pold = float('inf')
            continue

        F = L0-L; F[F<0] = 0                         # Bar forces (scalars)
        Fvec = F[:,None]/L[:,None].dot([[1,1]])*barvec # Bar forces (x,y components)
        Ftot = dense(bars[:,[0,0,1,1]],
                        np.repeat([[0,1,0,1]], len(F), axis=0),
                        np.hstack((Fvec, -Fvec)),
                        shape=(N, 2))
        Ftot[:nfix] = 0                              # Force = 0 at fixed points
        p += deltat*Ftot                             # Update node positions

        # 7. Bring outside points back to the boundary
        d = fd(p); ix = d>0                          # Find points outside (d>0)
        if ix.any():
            dgradx = (fd(p[ix]+[deps,0])-d[ix])/deps # Numerical
            dgrady = (fd(p[ix]+[0,deps])-d[ix])/deps # gradient
            dgrad2 = dgradx**2 + dgrady**2
            p[ix] -= (d[ix]*np.vstack((dgradx, dgrady))/dgrad2).T # Project

        # 8. Termination criterion: All interior nodes move less than dptol (scaled)
        if (np.sqrt((deltat*Ftot[d<-geps]**2).sum(1))/h0).max() < dptol:
            break

    # Clean up and plot final mesh
    p, t = fixmesh(p, t)

    return p, t

def uniform_mesh_on_unit_circle(h0):
    """Uniform Mesh on Unit Circle"""
    fd = lambda p: np.sqrt((p**2).sum(1))-1.0
    def huniform(p):
        """Implements the trivial uniform mesh size function h=1."""
        return np.ones(p.shape[0])
    return distmesh2d(fd, huniform, h0, (-1,-1,1,1),pfix=[[0,0]])
    #return distmesh2d(fd, huniform, h0, (-1,-1,1,1))


import math
def convert_K_to_RGB(colour_temperature):
    """
    Converts from K to RGB, algorithm courtesy of 
    http://www.tannerhelland.com/4435/convert-temperature-rgb-algorithm-code/
    """
    #range check
    if colour_temperature < 1000: 
        colour_temperature = 1000
    elif colour_temperature > 40000:
        colour_temperature = 40000
    
    tmp_internal = colour_temperature / 100.0
    
    # red 
    if tmp_internal <= 66:
        red = 255
    else:
        tmp_red = 329.698727446 * math.pow(tmp_internal - 60, -0.1332047592)
        if tmp_red < 0:
            red = 0
        elif tmp_red > 255:
            red = 255
        else:
            red = tmp_red
    
    # green
    if tmp_internal <=66:
        tmp_green = 99.4708025861 * math.log(tmp_internal) - 161.1195681661
        if tmp_green < 0:
            green = 0
        elif tmp_green > 255:
            green = 255
        else:
            green = tmp_green
    else:
        tmp_green = 288.1221695283 * math.pow(tmp_internal - 60, -0.0755148492)
        if tmp_green < 0:
            green = 0
        elif tmp_green > 255:
            green = 255
        else:
            green = tmp_green

        
