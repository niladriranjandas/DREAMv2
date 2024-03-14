
"""
Perform some action on an AtomSel. Available actions are the following
functions
   randomizeDomainPos
   randomizeVelocities
   translateFit

 and the following classes
   PrintPos
   SetProperty
   RandomizeVelocities
   TranslateFit

   The classes are used like
   AtomSel("name C").apply(PrintPos())

   or
   AtomSel("resname GLY").apply(SetProperty("residueName","G"))
"""

import random
from math import sqrt
import simulation
from simulationWorld import SimulationWorld_world as world
from selectTools import convertToSelectionIndices
from random import gauss
from atomSel import AtomSel
from atomSelAction import PyAtomSelAction
from atom import Atom

def maxwellDist(temperature,
                factor):
    range = sqrt( factor*temperature );
    return (gauss( 0, range ),
            gauss( 0, range ),
            gauss( 0, range ))

def randomizeVelocities(temperature,
                        sel="known"):
    """
    randomize atomic velocities according to a Maxwell distribution
    corresponding to the given temperature. The default selection
    contains all known atoms. 
    """
    from selectTools import convertToAtomSel
    sel = convertToAtomSel(sel)

    sel.apply( RandomizeVelocities(temperature) )
    pass

def translateFit(fitCoords,
                 atomSel=0):
    """Return array of atom positions which best-fit coordinates
    specified by fitCoords. fitCoords must have size sim.numAtoms()
    This routine does not change the current coordinates.
    the default atomSel is known atoms.
    """
    from selectTools import convertToAtomSel
    atomSel = convertToAtomSel(atomSel)

    trans = TranslateFit(fitCoords,atomSel).trans
    
    sim = atomSel.simulation()
    coord=sim.atomPosArr()
    for i in range( len(coord) ):
        coord[i] +=  trans
        pass
    return coord

#
#
#

class TranslateFit(PyAtomSelAction):
    """translate atom positions to best-fit coordinates specified by fitCoords.
    fitCoords must have size sim.numAtoms()"""
    def __init__(s,fitCoords,fitSel=0):
        """
        fitCoords- coordinates to use for fitting.
        fitSel   - subset of atoms to use in calculating fit (default: known)

        After construction, the fit is calculated and stored in member trans.
        """
        from vec3 import Vec3
        PyAtomSelAction.__init__(s,s)
        s.fitCoords = fitCoords
        if not fitSel: fitSel = AtomSel("known")
        if type(fitSel) == type("string"):
            fitSel = AtomSel(fitSel)
            pass
        s.atomSel = fitSel

        s.trans=Vec3(0,0,0)
        for atom in fitSel:
            s.trans += atom.pos() - fitCoords[atom.index()]
            pass
        s.trans *= -1.0 /len(fitSel)
        return
    def run(s,sim,index):
        sim.setAtomPos(index, sim.atomPos(index) + s.trans)
        return
    pass
    

class RandomizeVelocities(PyAtomSelAction):
    """randomize atomic velocities according to a Maxwell distribution
    corresponding to the given temperature.     """
    def __init__(s,temperature):
        PyAtomSelAction.__init__(s,s)
        s.temp = temperature
        return
    def run(s,sim,index):
        try:
            sim.setAtomVel(index,
                           maxwellDist(s.temp,
                                  world().kBoltzmann()/sim.atomMass(index)))
        except:
            print(index, world().kBoltzmann(), sim.atomMass(index))
            raise
        return
    pass

class PrintPos(PyAtomSelAction):
    """print position of atoms in selection"""
    def __init__(s,otherCoords=0):
        PyAtomSelAction.__init__(s,s)
        s.otherCoords = otherCoords
        return
    def init(s,sel):
        line = "%4s %4s %4s %4s " % ("segid", "resid", "rname", "name")
        line += "        position           "
        if ( s.otherCoords ):
            line += "    comparison position      "
            pass
        print(line)
        print('-'*len(line))
        return
    def run(s,sim,index):
        atom = Atom(sim,index)
        line = "%4s  %4d  %4s %4s " % (atom.segmentName(),
                                     atom.residueNum(),
                                     atom.residueName(),
                                     atom.atomName())
        line += "[%8.3f %8.3f %8.3f]" % tuple(atom.pos())
        if s.otherCoords:
            line += " [%8.3f %8.3f %8.3f]" % tuple(s.otherCoords[ atom.index() ])
        print(line)
        return
    pass

from atomSelAction import SetProperty
class SetPropertyArr(PyAtomSelAction):
    """set named property of all <m atom>.Atoms in a selection
    usage:
    sel = AtomSel("string")
    sel.apply(SetProperty(name,val)
    where
      name is a property of an atom (e.g. pos) and val is the value to set
      that property. If val is an sequence of length len(sel), values will be
      set appropriately.

    <m atomSelAction>.SetProperty is an optimized version of this class which
    allows atom properties top be set to a constant value.
    """
    def __init__(s,name,val):
        PyAtomSelAction.__init__(s,s)
        capName = name[0].upper()
        capName += name[1:]
        s.accessorName = "atom.set" + capName
        s.val = val
        s.count=-2
        import types
        if type(val) == tuple or type(val) == list:
            s.count=-1
            pass
        return
    def run(s,sim,index):
        atom = Atom(sim,index)
        if s.count<-1:
            cmd = s.accessorName + "(s.val)"
        else:
            cmd = s.accessorName + "(s.val[s.count])"
            s.count += 1
            pass
        exec( cmd )
        return
    pass

def getProperty(name,sel):
    """
    get array corresponding to the named property for the given selection.
    """
    if type(sel)==type(""):
        sel = AtomSel(sel)
        pass

    ret=[]
    cmd = "atom.%s()" % name
    for atom in sel:
        ret.append( eval(cmd ) )
        pass
    return ret    
    

def copyAtomCoords(sel1,sel2):
    """ copy atomic coordinates from sel2 to sel1.
    The two selections must contain the same number of atoms. This is most
    useful for copying coordinates between segments or between separate
    <m simulation>.Simulations. The two selections must contain the same number
    of atoms.

    The current implementation requires that identical ordering of the atoms
    in the two selections, but there may be gaps.
    """

    from atomSel import AtomSel
    if type(sel1)==type(""):
        sel1 = AtomSel(sel1)
        pass
    if type(sel2)==type(""):
        sel2 = AtomSel(sel2)
        pass

    if len(sel1) != len(sel2):
        raise Exception("selections are of different size")

    indices1=sel1.indices()
    indices2=sel2.indices()
    sim1 = sel1.simulation()
    sim2 = sel2.simulation()

    for cnt in range(len(indices1)):
        sim1.setAtomPos(indices1[cnt],sim2.atomPos(indices2[cnt]))
        pass
    return

def mass(selection="not PSEUDO"):
    """
    Return the mass of the specified <m atomSel>.AtomSel.
    """
    from selectTools import convertToAtomSel
    sel = convertToAtomSel(selection)
    mass = 0.
    for atom in sel:
        mass += atom.mass()
        pass
    return mass
    

def centerOfMass(sel="not PSEUDO",
                 useMass=True):
    """
    Return the center of mass of the specified atoms.

    The returned object is a <m vec3>.Vec3 instance representing the center of
    mass (if useMass is True) or centroid (if useMass is False) of the atoms
    specified by the argument sel.  sel can be either an atom selection string
    or an <m atomSel>.AtomSel instance.
    
    """
    from selectTools import convertToAtomSel
    sel = convertToAtomSel(sel)

    from atomSel import intersection
    from vec3 import Vec3
    center=Vec3(0,0,0)

    if len(intersection(AtomSel("not known") , sel)) > 0:
        raise Exception('atoms with undefined positions for center of mass/centroid calculation')
    elif len(sel) == 0:
        raise Exception('no atoms selected for center of mass/centroid calculation')
    else:
        mass=0
        for atom in sel:
            center += atom.mass() * atom.pos() if useMass else atom.pos()
            mass += atom.mass()
            pass
        center /= mass if useMass else len(sel)
        pass

    return center
    

from math import pi
def randomizeDomainPos(sel,
                       deltaPos=20,
                       deltaAngle=pi,
                       centerSel=None):
    """
    Randomly rotate and translate the given atom selection, 
    as a rigid-body about its center of mass.

    The deltaPos parameter specifies the maximum displacement to apply.

    The deltaAngle parameter specifies the maximal +/- angle (in radians) to
    rotate. 

    If centerSel is specified, atoms in sel at first centered around the average
    position of centerSel before rotation/translation.
    """

    from selectTools import convertToAtomSel
    sel = convertToAtomSel(sel)

    from vec3 import Vec3
    center = centerOfMass(sel)

    from atomSelAction import Rotate, Translate
    from mat3 import rotVector

    import random
    from math import cos, sin, pi, sqrt
    theta = random.uniform(-pi,pi)
    # not uniform in phi because that distribution is biased toward the poles
    # sinPhi is uniform, however.
    cosPhi = random.uniform(-1,1)
    sinPhi = sqrt(1-cosPhi**2)
    randVec = Vec3(cos(theta)*sinPhi,sin(theta)*sinPhi,cosPhi)
    rotMat = rotVector(randVec,random.uniform(-deltaAngle,deltaAngle))

    if centerSel:
        centerSel = convertToAtomSel(centerSel)
        sel.apply( Translate(-centerOfMass(sel)) )
        sel.apply( Translate(centerOfMass(centerSel)) )
        pass                   

    sel.apply( Rotate(rotMat,center) )

    sel.apply( Translate(Vec3(random.uniform(-deltaPos,deltaPos),
                              random.uniform(-deltaPos,deltaPos),
                              random.uniform(-deltaPos,deltaPos))) )

    return


def genRandomCoords(selection="not pseudo",
                    xyzInterval=10.0):
    """
    Assign atoms in the specified atom selection random x,y,z, coordinates in the
    interval -xyzInterval .. +xyzInterval.
    """
    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)

    from vec3 import Vec3
    for atom in selection:
        atom.setPos( Vec3(random.uniform(-xyzInterval,xyzInterval),
                          random.uniform(-xyzInterval,xyzInterval),
                          random.uniform(-xyzInterval,xyzInterval)) )
        pass
    return
