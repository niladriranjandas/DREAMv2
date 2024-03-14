
"""
Tools for molecular symmetry analysis
"""

class SymmetryViolation(Exception):
    def __init__(s,string):
        Exception.__init__(s,string)
        return
    pass

def testC2(sel1,
           sel2,
           rmsdTol=1e-2,
           transTol=1e-2,
           rotTol=0.1):
    """
    test for C2 Symmetry

    First, subtract off center of mass of sel1 and sel2.

    Swap atoms in sel1 and sel2 and rigid-body fit the result to the original.
    The resulting atomic rmsd should be less than rmsdTol.
    The norm of the translation should be less than transTol.
    The resulting rotation matrix should be within rotTol of 180 degrees.

    Atomic coordinates are not modified by this function.

    If any tolerance is violated, a SymmetryViolation exception is thrown.

    A SymmetryResults class is returned, containing:
    rmsd, trans, rot, and rotVec (the rotation axis)    

    """

    from selectTools import convertToAtomSel
    sel1=convertToAtomSel(sel1)
    sel2=convertToAtomSel(sel2)

    import simulation
    if sel1.simulation().name() != sel2.simulation().name():
        raise Exception("sel1 and sel2 must refer to the same Simulation")

    if sel1.size() != sel2.size():
        print("%d != %d" % (sel1.size(), sel2.size()))
        raise Exception("sel1 and sel2 must have the same size")

    
    #clone current Simulation, so atom positions aren't affected
    import simulation
    oldSimulation=simulation.currentSimulation()
    simulation.makeCurrent(sel1.simulation())

    from xplorSimulation import XplorSimulation
    simulation.makeCurrent(sel1.simulation())
    xsim=XplorSimulation(True)
    simulation.makeCurrent(xsim)


    from atomSel import AtomSel
    sel1=AtomSel(sel1.string())
    sel2=AtomSel(sel2.string())


    from atomSel import union
    fullSel = union(sel1,sel2)

    if len(fullSel)<4:
        raise Exception("too few selected atoms (%d)" % len(fullSel))
    


    from atomAction import centerOfMass
    cm=centerOfMass(fullSel)

    # subtract off center of mass
    for atom in fullSel: atom.setPos( atom.pos() - cm )

    compCoords=xsim.atomPosArr()

    #swap coords
    from vec3 import Vec3, norm, unitVec
    for i in range(len(sel1)):
        tmp=Vec3(sel1[i].pos())
        sel1[i].setPos( sel2[i].pos() )
        sel2[i].setPos( tmp )
        pass

    # perform rigid body fit
    from atomSelAction import Fit, RMSD
    fitter = Fit(compCoords,fullSel)
    fullSel.apply( fitter )
    measure=RMSD(compCoords)
    fullSel.apply(measure)

    class SymmetryResults:
        pass
    ret= SymmetryResults()

    ok=True
    from vec3 import Vec3, norm, unitVec
    ret.rmsd=measure.rmsd()
    if ret.rmsd>rmsdTol:
        print("testC2: atomic fit rmsd violated: %7.3f" % ret.rmsd)
        ok=False
        pass
    from vec3 import norm
    ret.trans=fitter.translation()
    fitNorm=norm(ret.trans)
    if fitNorm>transTol:
        print("testC2: fit translation violated: %7.3f" % fitNorm)
        ok=False
        pass
    from mat3 import rotationAmount
    from math import pi
    ret.theta=rotationAmount*180/pi
    if abs(ret.theta-180)>rotTol:
        print("testC2: 180 fit rotation violated: %7.3f" % ret.theta)
        ok=False
        pass
    
    #now determine rotation axis
    from mat3 import eigen
    e=eigen(fitter.rotation())

    for xyz in range(3):
        val=e[xyz].value()
        vec =e[xyz].vector()
        
        if abs(val.imag)<1e-14:
            ret.rotVec = unitVec(Vec3(vec[0].real,vec[1].real,vec[2].real))
            break
        pass
    simulation.makeCurrent(oldSimulation)

    if not ok:
        raise SymmetryViolation("one or more symmetry violations")
    
    return ret

