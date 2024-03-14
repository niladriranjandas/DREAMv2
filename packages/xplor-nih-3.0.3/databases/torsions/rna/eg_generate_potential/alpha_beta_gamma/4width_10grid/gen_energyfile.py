
import glob

import densityEstimation
import cdsVector
import cdsMatrix
import trace


# Torsion angles to consider.
torsion_names = ['alpha', 'beta', 'gamma']

# Global width for density estimation. 
width = 4

# Grid spacing. 
spacing = 10


# Non-numeric part of the energy file.
header = """name %s_%s_%s

atom1 resid _RESID-1 and name O3'
atom2 resid _RESID and name P
atom3 resid _RESID and name O5'
atom4 resid _RESID and name C5'

atom1 resid _RESID and name P
atom2 resid _RESID and name O5'
atom3 resid _RESID and name C5'
atom4 resid _RESID and name C4'

atom1 resid _RESID and name O5'
atom2 resid _RESID and name C5'
atom3 resid _RESID and name C4'
atom4 resid _RESID and name C3'

info B-factor<60. 
info Number of data points: %i
info Overall width for density estimation (deg.): %.2f
info Grid spacing for interpolation (deg.): %.2f

""" 


trace.suspend()

infile = open(glob.glob('../../*_allringatoms_epsilonfilter.csv')[0], 'r')
contents = infile.read()
infile.close()


torsion_ranges = [(-180, 180), (-180, 180), (-180, 180)]
torsion_domains = [max(x)-min(x) for x in torsion_ranges]

ndim = len(torsion_ranges)  # dimensionality

points = [[] for x in range(ndim)]
lines = contents.splitlines()
for line in lines:
    if line.startswith('#'):  # line with field names; get offsets
        allfields = [x.strip() for x in line.split('#')[1].split(',')]
        torsion1_idx = allfields.index(torsion_names[0])
        torsion2_idx = allfields.index(torsion_names[1])
        torsion3_idx = allfields.index(torsion_names[2])
    else: # data entry line
        words = line.split(',')
        torsion1 = float(words[torsion1_idx])
        torsion2 = float(words[torsion2_idx])
        torsion3 = float(words[torsion3_idx])

        torsions = [torsion1, torsion2, torsion3]
        for i in range(ndim):
            points[i].append(torsions[i])

ndata = len(points[0]) # number of data points
dataset = cdsMatrix.CDSMatrix_double(ndim, ndata).fromList(points)
dataset = densityEstimation.wraparound(dataset, torsion_domains)

kde = densityEstimation.adaptive_kde(dataset, width=width, wrapped=True)

# Number of points in each dimension of the grid.
npoints = [(x/spacing)+1 for x in torsion_domains if not x%spacing]
if len(npoints) != ndim: raise Exception, 'domain not multiple of spacing'

(grid, axes) = densityEstimation.grid(npoints, torsion_ranges, True, True)

# Add the last point to each grid axis. It is missing because of periodicity
# (i.e., it is missing from the actual grid), but the interpolation function
# requires it when inputting the axes regardless.
for axis in axes:
    axis.append(axis[-1]+spacing)

densities = kde.evaluate(grid)
energies = -cdsVector.log(densities)

filename = '%s_%s_%s_energy.dat' % (torsion_names[0], torsion_names[1],
                                    torsion_names[2])
outfile = open(filename, 'w')

# Write header.
outfile.write(header % (torsion_names[0], torsion_names[1], torsion_names[2],
                        ndata, width, spacing))



linesize = 10  # number of figures per line of numbers (axis1, ..., energy)

axis_count = 0
for axis in axes:
    count = 0
    axis_count += 1
    for value in axis:
        if count == 0:
            outfile.write('axis%i %.3f ' % (axis_count, value))
            count += 1
        elif count > 0 and count < linesize-1:
            outfile.write('%.3f ' % value)
            count += 1
        else:
            outfile.write('%.3f \n' % value)
            count = 0
    outfile.write('\n\n')


count = 0
for value in energies:
    if count == 0:
        outfile.write('energy %.3f ' % value)
        count += 1
    elif count > 0 and count < linesize-1:
        outfile.write('%.3f ' % value)
        count +=1
    else:
        outfile.write('%.3f \n' % value)
        count = 0
outfile.write('\n\n\n')

outfile.close()

trace.resume()



