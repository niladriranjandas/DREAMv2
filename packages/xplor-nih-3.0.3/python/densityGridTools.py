

def trimDensityGrid(dmap, pad=0.5, sel='not pseudo'):
    """Return a subregion of the input atom density map around atom coordinates.

    Given the current atomic coordinates defined by sel (an atom selection
    string) and an atomic density map, dmap (a string with the filename path of
    a ccp4 or mrc file or, alternatively, an <m atomDensity>.DensityGrid
    instance), return a subregion of dmap (an <m atomDensity>.DensityGrid
    instance) around the coordinates.  The subregion is defined by the smallest
    sphere (of radius r) that envelops all the coordinates, with r augmented by
    the fraction specified by the pad argument (i.e., the sphere has a radius
    r + pad * r, centered around the centroid of the coordinates).
    
    """
    # Convert sel to an atomSel.AtomSel instance.
    import selectTools
    sel = selectTools.convertToAtomSel(sel) # need atoms for both centroid and
                                        # distance calculation: convert once
    # Calculate centroid.
    import atomAction
    centroid = atomAction.centerOfMass(sel=sel, useMass=False)

    # Get distance from centroid of the atom farthest away from it.
    import vec3
    distances = [vec3.norm(atom.pos()-centroid) for atom in sel]
    distances.sort()

    # Radius of the smallest sphere that evelopes all atoms.
    radius = distances[-1]

    # Actual radius that determines subregion of map to use.
    radius = radius + pad * radius

    # If input dmap is a filename, create atomDensity.DensityGrid instance. 
    if type(dmap) is str:
        from atomDensity import DensityGrid
        filename = dmap
        dmap = DensityGrid()
        dmap.readCCP4(filename, verbose=True)

    # Total no. of points for each dimension in output map.
    numx = int(round(2*radius/dmap.xdelta))
    numy = int(round(2*radius/dmap.ydelta))
    numz = int(round(2*radius/dmap.zdelta))

    # For each dimension, get point in input map that will be the first point 
    # in output map.
    startx = int(round( (centroid[0]-radius - dmap.xmin)/dmap.xdelta ))
    starty = int(round( (centroid[1]-radius - dmap.ymin)/dmap.ydelta ))
    startz = int(round( (centroid[2]-radius - dmap.zmin)/dmap.zdelta ))

    # Output map.
    result = DensityGrid(numx, numy, numz)
    result.xdelta = dmap.xdelta 
    result.ydelta = dmap.ydelta 
    result.zdelta = dmap.zdelta

    result.xmin = dmap.xmin + startx * dmap.xdelta 
    result.ymin = dmap.ymin + starty * dmap.ydelta 
    result.zmin = dmap.zmin + startz * dmap.zdelta 

    for i in range(numx):
        for j in range(numy):
            for k in range(numz):
                result.setData(i, j, k,
                               dmap.getData(startx+i, starty+j, startz+k))
    return result
        
