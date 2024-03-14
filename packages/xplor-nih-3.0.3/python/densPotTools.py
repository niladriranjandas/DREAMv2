
def create_DensPot(name,selection="not PSEUDO",targetDens=1.1,rSolvent=0):
    """
    create a DensPot including all atoms in selection. targetDens specifies
    the target density value used in energy calculation.

    The target density is given in units of amu/A^3. For rSolvent=0 we have
    found that the density is 1.101 +/- 0.010 for 30 high-resolution
    crystal structures.
    """

    import types
    from atomSel import AtomSel
    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)

    heavyAtoms = AtomSel("(%s) and (not name H*)" % selection.string())

    mass=0
    for atom in heavyAtoms:
        mass += atom.mass()
        # make sure we include hydrogens - they are included in the volume
        # calculation (radiusType='heavyatom')
        for prot in AtomSel("bondedto id %d and name H*" %
                            (atom.index()+1)):
            mass += prot.mass()
            pass
        pass

#    print selection.string(), mass


    pot = DensPot(name,heavyAtoms,targetDens,mass=mass)

    return pot

from pyPot import PyPot
class DensPot(PyPot):
    def __init__(s,name,selection,target=1.1,radiusType='heavyatom',
                 radii=None,mass=None):
        """
        target is the target density in units of amu/A^3
        """

        from selectTools import convertToAtomSel
        selection = convertToAtomSel(selection)

        s.sel = selection

        PyPot.__init__(s,name)
        from surfaceTools import create_SurfaceArea

        s.sa = create_SurfaceArea(rSolvent=0,
                                  sel=selection,
                                  radiusType=radiusType)

        if radii:
            for atom in selection:
                s.sa.setRadius(atom, radii[atom.index()])
                pass
            pass

        #s.sa.setVerbose(4)
        s.sa.setUseConvexHull(False)
        s.sa.setMoveTol(0)
        s.sa.setVerbose(0)
        #s.sa.setUseSphereMetric(1)

        from volume import Volume
        s.vol = Volume(s.sa)

        if mass:
            s.mass=mass
        else:
            s.mass=0
            for atom in selection:
                s.mass += atom.mass()
                pass
            pass

        s.target = target
        
        pass
    def __del__(self,destroy=0):
        PyPot.__del__(self)

    def calcEnergy(s):
        return s.calcEnergyAndDerivList()

    def calcEnergyAndDerivList(s,derivs=None):

        sim = s.sel.simulation()
        if derivs:
            from derivList import DerivList
            dList = DerivList()
            dens = s.density(dList)
            derivs[sim] += dList[sim] * s.scale() * (dens - s.target)
        else:
            dens = s.density()
            pass

            
        return 0.5 * s.scale() * (dens - s.target)**2
        
    def density(s,derivs=None):
        """calculate molecular density based on the surface defined by the
        sa member. If derivs is specified, the gradient of the density
        with respect to each atom will be calculated.
        """

        timing=0
        if timing:
            import time
            start=time.clock() ; print("start...", end=' ')
            pass
        if derivs:
            from vec3 import Vec3
            vol = s.vol.calc(derivs)
        else:
            vol = s.vol.calc()
            pass
        
        if timing:
            print("stop",time.clock()-start)
            pass
        dens = s.mass / vol

        if derivs:
            sim = s.sel.simulation()
            derivs[sim].scale( -s.mass / vol**2)
            pass

        return dens

    pass

# analysis:
#   report density
#   might also report per-residue density

        
