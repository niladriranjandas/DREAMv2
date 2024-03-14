"""Module with utilities for density estimation.

When using the kernel density estimation capabilities of this module, please
cite:

Bermejo, G.A., Clore, G.M., and Schwieters, C.D. (2012). Smooth statistical
torsion angle potential derived from a large conformational database via
adaptive kernel density estimation improves the quality of NMR protein
structures. Protein Sci. 21, 1824-1836.

"""
# Written by Guillermo A. Bermejo.


import math
import itertools

import cdsMatrix
import cdsVector

import pdfCore


class KDE:
    """General kernel density estimate with a fixed- or variable-width kernel.

    Subclass and implement the method evaluate according to a specific kernel
    function (e.g., class NormalKDE(KDE) evaluates the estimates using a
    multivariate normal kernel).  Optionally, the kernel-specific method
    optimal_width, that calculates an optimal width assuming an underlying
    probability density function (pdf), can also be implemented.
    
    The general formula used for density estimation is that given by (Wand and
    Jones, 1993), modified to include local widths as in (Silverman, 1986):

    f(x) = n^-1 SUM_i{h^-d |C|^-0.5 lambdai^-d K[(h lambdai)^-1 C^-0.5 (x-Xi)]}
 
    where x and Xi are d-dimensional vectors (Xi represents sample data point
    i), n is the total number of points, C is a symmetric positive definite
    d-by-d matrix with determinant |C|, and h and lambdai are, respectively, the
    global and local widths of kernel K at Xi.  The following discussion will
    assume kernel K is a symmetric pdf (e.g., multivariate standard normal), and
    refer to the expression within the sum as the "effective kernel" (which
    might not be symmetric).

    The constructor requires a dataset, a <m cdsMatrix>.CDSMatrix_type(d, n)
    (type=int, double) that contains the sample points from which the density is
    estimated.  For example, for two 3D points with coordinates (1.0, 2.0, 3.0)
    and (4.0, 5.0, 6.0):

        dataset = cdsMatrix.CDSMatrix_double(3, 2).fromList([[1.0, 4.0],
                                                             [2.0, 5.0],
                                                             [3.0, 6.0]])

    The constructor call is something like:
    
        kde = NormalKDE(dataset)
        
    Optionally, the arguments width, localwidths, symmetric, cov, and wrapped
    can also be specified.  width is the global kernel width h; if the method
    optimal_width is implemented, an optimal width is automatically calculated
    when width is not supplied.  localwidths is a sequence of length n that
    specifies a local width for each kernel used in the density estimation:
    localwidths[i]=lambdai; if not supplied lambdai=1 for all i (i.e., a fixed-
    width effective kernel is used).  Different combinations of values for
    arguments symmetric and cov determine matrix C:
    
    symmetric=True => C is the identity matrix;
        
    symmetric=False and cov=None => C is the covariance matrix of the data
    (the Fukunaga method; see Silverman 1986, Sec. 4.2.1, p. 77);
        
    symmetric=False and cov=A => C is the arbitrary symmetric positive definite
    matrix A, input as either a <m cdsMatrix>.CDSMatrix_type(d, d) (type=int,
    double) or a sequence of diagonal elements for a diagonal A.

    Thus, symmetric=True yields a radially symmetric effective kernel, while
    options with symmetric=False yield an asymmetric one (except for trivial
    cases such as setting cov=A=constant*identity matrix). 
    The remaining argument, wrapped, specifies whether the input dataset has
    been augmented by surrounding it with one copy of it (e.g., by using
    function wraparound in this module).  

    The evaluate method in a KDE subclass should be coded so that density
    estimates at arbitrary points are obtained via:

        vals = kde.evaluate(pts)
        
    pts is a <m cdsMatrix>.CDSMatrix_type(d, m) (type=int, double), where m is
    the total number of points.  For example, if kde is to be evaluated at
    (1.5, 2.5, 3.5), (4.5, 5.5, 6.5), and (7.5, 8.5, 9.5):
    
        pts = cdsMatrix.CDSMatrix_double(3, 3).fromList([[1.5, 4.5, 7.5],
                                                         [2,5, 5.5, 8.5],
                                                         [3.5, 6.5, 9.5]])
                                                         
    The returned values, vals, is a <m cdsVector>.CDSVector_double of size m,
    where vals[i] is the density at the point in column i of pts (e.g., vals[1]
    is the density at (4.5, 5.5, 6.5) in the above example).
    
    References:
    Silverman, B.W. (1986).  Density Estimation for Statistics and Data
    Analysis.  Chapman & Hall.
    Wand, M.P. and Jones, M.C. (1993).  Comparison of Smoothing Parametrizations
    in Bivariate Kernel Density-Estimation. J. Amer. Statistical Assoc. 88(422):
    520-528.

    """
    def __init__(self, dataset, width=None, localwidths=None, symmetric=True,
                 cov=None, wrapped=False):
        self.dataset = cdsMatrix.CDSMatrix_double(dataset)  # ensure float
        self.d, self.n = (self.dataset.rows(), self.dataset.cols())
                
        self.symmetric = symmetric
        self.wrapped = wrapped
    
        # Set symmetric positive definite matrix cov (C in formula).
        if symmetric == True:  # Unit matrix.
            unitmatrix = cdsMatrix.CDSMatrix_double(self.d, self.d, 0.0)
            unitmatrix.setDiag(1.0)
            self.cov = unitmatrix
        elif cov == None:      # Sample covariace matrix of the data.
            self.cov = self.get_samplecov() 
        else:                  # Input matrix.
            if type(cov) not in (cdsMatrix.CDSMatrix_double,# input cov is a seq
                                 cdsMatrix.CDSMatrix_int):
                if len(cov) != self.d:
                    msg = 'Mismatch between lengh of cov sequence (%i) and \
dimensionality of data points (%i)'
                    raise Exception(msg % (len(cov), self.d))
                diag = cdsVector.CDSVector_double(cov)
                cov = cdsMatrix.CDSMatrix_double(self.d, self.d, 0.0)
                cov.setDiag(diag)
                self.cov = cov
            else:  # input cov is a matrix
                if cov.rows() != self.d or cov.cols() != self.d:
                    msg = 'Mismatch between dimensionality of cov matrix and \
data points'
                    raise Exception(msg)
                self.cov = cdsMatrix.CDSMatrix_double(cov)  # ensure float

        self.invcov = cdsMatrix.inverse(self.cov) # inverse of cov matrix

        # Set global kernel width.
        if width != None:
            self.width = float(width)   # subjective input width
        else:
            self.width = self.optimal_width()    # optimal width

        # Set local widths.
        if localwidths == None:
            self.localwidths = cdsVector.CDSVector_double(self.n, 1.0)
        else:
            self.localwidths = cdsVector.CDSVector_double(localwidths)

        if len(self.localwidths) != self.n:
            msg = 'Mismatch in number of localwidths (%i) and data points (%i)'
            raise Exception(msg % (len(self.localwidths), self.n))

    def get_samplecov(self):
        """Return a <m cdsMatrix>.CDSMatrix_double with the sample covariance matrix.

        """
        average = []   # average vector
        for dim in range(self.d):
            ave = 0.0  # average value in dimension dim
            for point in range(self.n):
                ave += self.dataset[dim, point]
            ave = ave/float(self.n)
            average.append(ave)

        diff = cdsMatrix.CDSMatrix_double(self.dataset)
        for dim in range(self.d):
            for point in range(self.n):
                diff[dim, point] = diff[dim, point] - average[dim]

        samplecov = diff * cdsMatrix.transpose(diff)
        samplecov.scale(1.0/float(self.n-1))
        return samplecov

    def get_widths(self):
        """Return a <m cdsVector>.CDSVector_double of localwidth*width at each point.
        """
        result = cdsVector.CDSVector_double(self.localwidths)
        result.scale(self.width)
        return result

    def optimal_width(self):
        """Return kernel-specific optimal width assumming an underlying pdf.

        The implementaton of this method in a subclass is not needed if width is
        to be set only subjectively (i.e., specified upon intance creation).
        Raise ValueError if this method is called while not set in the subclass.

        """
        msg = 'Unspecified method to set optimal width; enter subjective value'
        raise ValueError(msg)

    def ndim(self):
        """Return the dimensionality of the estimated density function."""
        return self.d

    def npts(self):
        """Return the number of sample points used for density estimation."""
        return self.n

    def evaluate(self, points):
        """Return kernel-specific density estimates at arbitrary points."""
        msg = 'Unspecified method to evaluate kernel density estimates'
        raise ValueError(msg)



class NormalKDE(KDE):
    """Kernel density estimate with a multivariate normal kernel.

    It uses kernel

        K(x) = (2 pi)^-(d/2) exp(-0.5 x^T x)

    where x is a d-dimensional vector with transpose x^T.  It estimates the
    density using the general formula given in the KDE class.  See KDE class for
    a detailed description of the interface.

    """
    def optimal_width(self):
        """Return optimal width assuming the underlying pdf is normal.

        If symmetric=False it returns hopt = [4/n(2d+1)]^[1/(d+4)].
        If symmetric=True it returns sqrt(Tr{S}/d)*hopt, where Tr{S} is the
        trace of the covariance matrix of the data (Silverman, 1986, Sec. 4.3.2,
        pp.86, 87).

        Reference: Silverman, B.W. (1986). Density Estimation for Statistics
        and Data Analysis. Chapman & Hall.
        
        """
        n = float(self.n)
        d = float(self.d)
        hopt = (4.0 / (n * (2.0*d + 1.0)))**(1.0 / (d + 4.0))
        if self.symmetric == True:
            sigma = math.sqrt(cdsMatrix.trace(self.get_samplecov())/d)
            return sigma*hopt
        else:
            return hopt
        
    def evaluate(self, points):
        """Return a <m cdsVector>.CDSVector_double of density estimates at points.

        points is a <m cdsMatrix>.CDSMatrix_type(d, m) (type=int, double), where
        m is the total number of points.  For example, a NormalKDE instance is
        evaluated at (1.5, 2.5, 3.5), (4.5, 5.5, 6.5), and (7.5, 8.5, 9.5) as
        follows:
        
        points = cdsMatrix.CDSMatrix_double(3, 3).fromList([[1.5, 4.5, 7.5],
                                                            [2.5, 5.5, 8.5],
                                                            [3.5, 6.5, 9.5]])
                                                             
        The returned values vector, vals, has size m, and vals[i] is the density
        estimate at the point in column i of the input points matrix (e.g.,
        vals[1] is the density at (4.5, 5.5, 6.5) in the above example).

        """
        points = cdsMatrix.CDSMatrix_double(points)  # ensure float
        result = pdfCore.pdfEvaluate(points, self.dataset, self.invcov,
                                     self.width, self.localwidths)
        if self.wrapped == True:
            n = float(self.n/(3**self.d)) #number of points in unwrapped dataset
        else:
            n = float(self.n)
        d = float(self.d)

        # (Note that division by sqrt(detC) was already done inside pdfCore.)
        factor = 1.0/(n * self.width**d * (2.0*math.pi)**(d/2.0))
        result.scale(factor)
        return result

    def marginalize(self, dims):
        """Return a NormalKDE instance with the marginal density over dims.

        dims argument is either an integer or a sequence of intergers specifying
        the dimension(s) over which to marginalize.  
        
        """
        if type(dims) not in [list, tuple]: dims = [dims]

        dims = set(dims)  # side effect: removes repetitions
        jointdims = set(range(self.d))
        
        if not dims.issubset(jointdims):
            raise ValueError('Input dimensions are not a subset of the joint')

        keep = list(jointdims.difference(dims))  # dimensions to keep
        keep.sort() # offsets will be used to label new dimensions

        dataset = cdsMatrix.CDSMatrix_double(len(keep), self.n)
        if self.symmetric == True:
            cov = None  # set to unit on instantiation with symmetric=True
        else:  # to be calculated below
            cov = cdsMatrix.CDSMatrix_double(len(keep), len(keep))
            
        for i in keep:
            # New dataset.
            for point in range(self.n):
                dataset[keep.index(i), point] = self.dataset[i, point]
            # New cov matrix.
            if self.symmetric == False:
                for j in keep:
                    cov[keep.index(i), keep.index(j)] = self.cov[i, j]

        return NormalKDE(dataset, width=self.width,localwidths=self.localwidths,
                         symmetric=self.symmetric,cov=cov, wrapped=self.wrapped)
                
                

### Classes for other kernels should go here ###


def adaptive_kde(dataset, width=None, symmetric=True, cov=None, wrapped=False,
                 alpha=0.5, pilot=NormalKDE, adaptive=NormalKDE):
    """Return a class 'adaptive' instance: the adaptive kernel density estimate.

    Arguments dataset, width, symmetric, cov, and wrapped are those required for
    KDE instantiation; alpha is the sensitivity parameter used in the local
    width calculation from the pilot density estimate:

        lambdai = [pilot_f(Xi)/g]^-alpha

    with lambdai the local width of the kernel centered at vector Xi, where the
    pilot desity is pilot_f(Xi), and log[g] = n^-1 SUM_i{log[pilot_f(Xi)]} (eq.
    5.7 in Silverman 1986, Sec. 5.3, p.101).  The kernels used in the pilot and
    adaptive stages are respectively specified with arguments pilot and
    adaptive, which expect the names of the corresponding classes; the global
    width of both kernels is the same, specified with the width argument.  The
    overall procedure is described by (Silverman, 1986, Sec. 5.3); for more
    details see KDE docstring.

    Reference: Silverman, B.W. (1986).  Density Estimation for Statistics and
    Data Analysis.  Chapman & Hall.

    """
    # Pilot kernel estimate.
    pilot = pilot(dataset, width=width, symmetric=symmetric, cov=cov,
                  wrapped=wrapped)

    # Calculate local widths from pilot density estimate.
    localwidths = pilot.evaluate(dataset)
    log_g = cdsVector.sum(cdsVector.log(localwidths))/float(pilot.n)
    localwidths.scale(1.0/math.exp(log_g))
    localwidths = localwidths**(-float(alpha))

    # Adaptive kernel estimate.
    adaptive = adaptive(dataset, width=width, localwidths=localwidths,
                        symmetric=symmetric, cov=cov, wrapped=wrapped)
    return adaptive
    


def wraparound(dataset, spectral_widths):
    """Return an augmented dataset by adding shifted copies of it.

    Impose a periodic or "wrap around" boundary condition on the dataset by
    augmenting it via the addition of shifted copies of it (Silverman, 1986,
    Sec. 2.10, p. 31, 32).

    dataset is <m cdsMatrix>.CDSMatrix_type(d, n) (type=int, double), where d is
    the dimensionality and n the number of points (like in KDE instantiation).
    For example, two 3D points with coordinates (1.0, 2.0, 3.0) and
    (4.0, 5.0, 6.0) yield:
    
        dataset = cdsMatrix.CDSMatrix_double(3, 2).fromList([[1.0, 4.0],
                                                             [2.0, 5.0],
                                                             [3.0, 6.0]])
                                                             
    After wrapping around the number of points is (3^d)n.
    spectral_widths is a sequence with the size of the domain in each dimension.
    For example, if dataset contains directions or angles defined within the
    range [-pi, pi], the associated size in spectral_widths is 2pi.  Sizes in
    spectral_widths and dimensions in dataset are matched by offset: row i in
    dataset corresponds to spectral_widths[i]; ValueError is raised if the
    dimensionality of dataset and spectral_widths do not match.

    Reference: Silverman, B.W. (1986).  Density Estimation for Statistics and
    Data Analysis.  Chapman & Hall.

    """
    n = dataset.cols()       # number of points in dataset
    d = dataset.rows()       # dimensionality of dataset
    m = len(spectral_widths) # dimensionality of spectral_widths
    if m != d:
        msg = 'spectral_widths has %s dimension(s), and dataset %s' % (m, d)
        raise ValueError(msg)

    axes = []
    for dim in range(d):
        axis = [[], [], []]
        for point in range(n):
            axis[0].append(dataset[dim, point])
            axis[1].append(dataset[dim, point] - spectral_widths[dim])
            axis[2].append(dataset[dim, point] + spectral_widths[dim])            
        axes.append(axis)
    
    result = [] 
    for dim in range(d):
        result.append([])
    for points in itertools.product(*axes):
        count = 0
        for axis in points:
            result[count].extend(axis)
            count += 1
    result = cdsMatrix.CDSMatrix_double(d, (3**d)*n).fromList(result)
    return result


def grid(npoints, ranges, periodic=False, axes=False):
    """Return a <m cdsMatrix>.CDSMatrix_double representing a d-dimesional grid.

    npoints is a sequence with the number of points in each dimension.  ranges
    is a sequence of (min, max) pairs, with the min and max values in each
    dimension (the order of values within a pair is irrelevant).  The number of
    dimensions, d, is len(npoints) (or len(ranges); ValueError is raised if
    len(npoints)!=len(ranges)).  For example:

        npoints = (5, 3)
        ranges = ([-10, 10], [-5, 5])

    represents a system with d=2, where the first dimension or axis has values
    in the [-10, 10] range, and the second axis in the [-5, 5].  5 points are
    requested from the first axis and 3 from the second; points are equally
    spaced to make the grid, i.e.

        Axis 1: [-10.0, -5.0, 0.0, 5.0, 10.0]
        Axis 2: [-5.0, 0.0, 5.0]

    The grid can be generated via grid(npoints, ranges), which, following the
    above example, returns a <m cdsMatrix>.CDSMatrix_double(2, 5*3) representing

        [[-10, -10, -10, -5, -5, ..., 10]   # axis 1
         [ -5,   0,   5, -5,  0, ...,  5]]  # axis 2
                                                        
    The optional argument periodic is a sequence with True/False for each axis;
    if True, the last point in the corresponding axis is omitted.  For example,
    grid(npoints, ranges, (True, False)) returns a grid with 

        Axis 1: [-10.0, -5.0, 0.0, 5.0]
        Axis 2: [-5.0, 0.0, 5.0]

    Alternatively, all axes can be given the same True or False periodic value
    by setting periodic=True or periodic=False, respectively, i.e.,

        grid(npoints, ranges, (True, True)) <-> grid(npoints, ranges, True) 
        grid(npoints, ranges, (False, False)) <-> grid(npoints, ranges, False)

    The axes may be obtained by setting argument axes=True, in which case the
    function returns a tuple with the grid (the <m cdsMatrix>.CDSMatrix_double)
    as first element and the axes, packed into a list, as second element.  For
    example, grid(npoints, ranges, True, True) returns

        (grid(npoints, ranges, True), [[-10.0, -5.0, 0.0, 5.0], [-5.0, 0.0]])

    i.e., (the grid, [axis 1, axis 2]).

    """
    d = len(ranges)  # number of dimensions
    m = len(npoints)
    if d != m:
        msg = 'ranges has dimension %i, and npoints %i' % (d, m)
        raise ValueError(msg)
    if periodic != True and periodic != False:
        m = len(periodic)
        if d != m:
            msg = 'ranges has dimension %i, and periodic %i' % (d, m)
            raise ValueError(msg)
    
    if periodic == True: periodic = d * [True]
    if periodic == False: periodic = d * [False]
    # Otherwise leave periodic as is.

    want_axes = axes  # free 'axes' variable 

    grid = []  # for cdsMatrix.CDSMatrix_double instantiation
    axes = []  # to contain the values at each axis
    for dim in range(d):
        grid.append([])
        axes.append([])

    dimensions = list(zip(npoints, ranges))
    dim = 0
    for (npoint, bounds) in dimensions:
        increment = (max(bounds) - min(bounds))/float(npoint-1)
        for point in range(npoint):
            value = min(bounds) + point*increment
            if periodic[dim] and point == npoint-1: break  # skip last point
            axes[dim].append(value)
        dim += 1

    for point in itertools.product(*axes):
        for dim in range(d):
            grid[dim].append(point[dim])
    grid = cdsMatrix.CDSMatrix_double(d, len(grid[0])).fromList(grid)
    
    if want_axes == False:
        return grid
    else:
        return (grid, axes)
    

def gengrid(npoints, ranges, periodic=False):
    """Return a Python itertools.product object representing a d-dimesional grid.

    Note: This is a Python generator version of the grid function; another
    difference with grid is that it lacks the functionality associated with the
    axes argument (there is no such argument here).

    npoints is a sequence with the number of points in each dimension.  ranges
    is a sequence of (min, max) pairs, with the min and max values in each
    dimension (the order of values within a pair is irrelevant).  The number of
    dimensions, d, is len(npoints) (or len(ranges); ValueError is raised if
    len(npoints)!=len(ranges)).  For example:

        npoints = (5, 3)
        ranges = ([-10, 10], [-5, 5])

    represents a system with d=2, where the first dimension or axis has values
    in the [-10, 10] range, and the second axis in the [-5, 5].  5 points are
    requested from the first axis and 3 from the second; points are equally
    spaced to make the grid, i.e.

        Axis 1: [-10.0, -5.0, 0.0, 5.0, 10.0]
        Axis 2: [-5.0, 0.0, 5.0]

    The grid associated with the above axes is the matrix

          -10, -10, -10, -5, -5, ..., 10      # from axis 1
           -5,   0,   5, -5,  0, ...,  5      # from axis 2

    The matrix columns are the multidimensional (in this example two-
    dimensional) points in the grid.  They are returned via a Python
    itertools.product object, a generator that yields each point packed into a
    tuple.
                                                        
    The optional argument periodic is a sequence with True/False for each axis;
    if True, the last point in the corresponding axis is omitted.  For example,
    gengrid(npoints, ranges, (True, False)) is associated with a grid generated
    from axes

        Axis 1: [-10.0, -5.0, 0.0, 5.0]
        Axis 2: [-5.0, 0.0, 5.0]

    Alternatively, all axes can be given the same True or False periodic value
    by setting periodic=True or periodic=False, respectively, i.e.,

    gengrid(npoints, ranges, (True, True)) <-> gengrid(npoints, ranges, True) 
    gengrid(npoints, ranges, (False, False)) <-> gengrid(npoints, ranges, False)

    """
    d = len(ranges)  # number of dimensions
    m = len(npoints)
    if d != m:
        msg = 'ranges has dimension %i, and npoints %i' % (d, m)
        raise ValueError(msg)
    if periodic != True and periodic != False:
        m = len(periodic)
        if d != m:
            msg = 'ranges has dimension %i, and periodic %i' % (d, m)
            raise ValueError(msg)
    
    if periodic == True: periodic = d * [True]
    if periodic == False: periodic = d * [False]
    # Otherwise leave periodic as is.

    axes = []  # to contain the values at each axis
    for dim in range(d):
        axes.append([])

    dimensions = list(zip(npoints, ranges))
    dim = 0
    for (npoint, bounds) in dimensions:
        increment = (max(bounds) - min(bounds))/float(npoint-1)
        for point in range(npoint):
            value = min(bounds) + point*increment
            if periodic[dim] and point == npoint-1: break  # skip last point
            axes[dim].append(value)
        dim += 1

    return itertools.product(*axes)



    
    

