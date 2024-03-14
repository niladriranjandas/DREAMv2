#
# Two potential terms involving ratios of cos(theta)
#
# this was developed to restrain rdc/csa Da values to be close to oneanother
#  (but not necessarily identical). These potentials have been reimplemented
#  the C++ module cosRatioPot for speed.
#
###################################################
#
# potential term #1:
#          1/2 * scale * (cos(theta1) - ratio*cos(theta2))^2
#
# where theta1 is one bond angle, theta2 is a second bond angle, and
# scale and ratio are constants.
#
# use:
#  from potCosAngleRatio import CosAngleRatioPot
#
#  term = CosAngleRatioPot("name",             
#                         AtomSel("atom 01"),
#                         AtomSel("atom a1"),
#                         AtomSel("atom b1"),
#                         AtomSel("atom 02"),
#                         AtomSel("atom a2"),
#                         AtomSel("atom b2"),
#                         [simulation])
#
# theta1 is given by the angle a1 - 01 - b1
# theta2 is given by the angle a2 - 02 - b2
#  [simulation is an optional simulation parameter - usually not needed.]
#
# to set the scale factor:
#   term.scale = 100
#
# to set the constant ratio:
#   term.ratio = 0.2
#
# remember to add this to the terms to be evaluated:
#   simulation.addPot(term)
#
###################################################
#
# potential term #2:
#       1/2 * scale * (cos(theta1)*cos(theta4) - cos(theta2)*cos(theta3))^2
#
# where theta1, theta2, theta3, and theta4 are bond angles, and
# scale is a constant. Note that the minimum of this potential occurs when
#
#    cos(theta1)/cos(theta2) = cos(theta3)/cos(theta4)
#
# use:
#  from potCosAngleRatio import CosAngleRatioPot
#
#  term = Cos2AngleRatioPot("name",             
#                         AtomSel("atom 01"),
#                         AtomSel("atom a1"),
#                         AtomSel("atom b1"),
#                         AtomSel("atom 02"),
#                         AtomSel("atom a2"),
#                         AtomSel("atom b2"),
#                         AtomSel("atom 03"),
#                         AtomSel("atom a3"),
#                         AtomSel("atom b3"),
#                         AtomSel("atom 04"),
#                         AtomSel("atom a4"),
#                         AtomSel("atom b4"),
#                         [simulation])
#
# theta1 is given by the angle a1 - 01 - b1
# theta2 is given by the angle a2 - 02 - b2
#  [simulation is an optional simulation parameter - usually not needed.]
#
# to set the scale factor:
#   term.scale = 100
#
# to set the constant ratio:
#   term.ratio = 0.2
#
#
from pyPot import PyPot
from math import sqrt
from atomSel import AtomSel
import simulation

def norm(v):
    return sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def dot(v1,v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def unitVec(v):
    return vecScale(1/norm(v) ,v)

def subVec(v1,v2):
    """subtract v2 from v1"""
    return (v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2])

def addVec(v1,v2):
    return (v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2])

def vecScale(a,v):
    return (a*v[0],a*v[1],a*v[2])


class CosAngleRatioPot(PyPot):
    def __init__(s,name,atomSel01,atomSela1,atomSelb1,
                 atomSel02,atomSela2,atomSelb2,sim=0):
        PyPot.__init__(s,name,s)
        s.id01 = atomSel01.indices()[0]
        s.ida1 = atomSela1.indices()[0]
        s.idb1 = atomSelb1.indices()[0]
        s.id02 = atomSel02.indices()[0]
        s.ida2 = atomSela2.indices()[0]
        s.idb2 = atomSelb2.indices()[0]
        if sim==0:
            sim = simulation.currentSimulation()
            pass
        s.sim = sim
        
        s.scale = 1
        s.ratio = 1
        return
    def calcEnergy(s):
        q01 = s.sim.atomPos(s.id01).tuple()
        qa1 = s.sim.atomPos(s.ida1).tuple()
        qb1 = s.sim.atomPos(s.idb1).tuple()
        q02 = s.sim.atomPos(s.id02).tuple()
        qa2 = s.sim.atomPos(s.ida2).tuple()
        qb2 = s.sim.atomPos(s.idb2).tuple()

        s.va1 = unitVec(subVec(qa1,q01))
        s.vb1 = unitVec(subVec(qb1,q01))
        s.na1 = norm(subVec(qa1,q01))
        s.nb1 = norm(subVec(qb1,q01))
        s.cos1 = dot(s.va1,s.vb1)
        
        s.va2 = unitVec(subVec(qa2,q02))
        s.vb2 = unitVec(subVec(qb2,q02))
        s.na2 = norm(subVec(qa2,q02))
        s.nb2 = norm(subVec(qb2,q02))
        s.cos2 = dot(s.va2,s.vb2)

        return 0.5 * s.scale * (s.cos1 - s.ratio*s.cos2)**2
    def calcEnergyAndDerivs(s,derivs):
        energy = s.calcEnergy()

        pref = s.scale * (s.cos1 - s.ratio*s.cos2)

        da1 = vecScale(pref/s.na1, subVec(s.vb1 ,
                                          vecScale(dot(s.va1,s.vb1),s.va1)))
        db1 = vecScale(pref/s.nb1, subVec(s.va1 ,
                                          vecScale(dot(s.va1,s.vb1),s.vb1)))
        derivs[s.ida1] = addVec(derivs[s.ida1], da1)
        derivs[s.idb1] = addVec(derivs[s.idb1], db1)
        derivs[s.id01] = subVec(derivs[s.id01], addVec(da1,db1))
        
        pref = -s.ratio * s.scale * (s.cos1 - s.ratio*s.cos2)

        da2 = vecScale(pref/s.na2, subVec(s.vb2 ,
                                          vecScale(dot(s.va2,s.vb2),s.va2)))
        db2 = vecScale(pref/s.nb2, subVec(s.va2 ,
                                          vecScale(dot(s.va2,s.vb2),s.vb2)))
        derivs[s.ida2] = addVec(derivs[s.ida2], da2)
        derivs[s.idb2] = addVec(derivs[s.idb2], db2)
        derivs[s.id02] = subVec(derivs[s.id02], addVec(da2,db2))
        return energy
    pass

class Cos2AngleRatioPot(PyPot):
    def __init__(s,name,
                 atomSel01,atomSela1,atomSelb1,
                 atomSel02,atomSela2,atomSelb2,
                 atomSel03,atomSela3,atomSelb3,
                 atomSel04,atomSela4,atomSelb4,
                 sim=0):
        PyPot.__init__(s,name,s)
        s.id01 = atomSel01.indices()[0]
        s.ida1 = atomSela1.indices()[0]
        s.idb1 = atomSelb1.indices()[0]
        s.id02 = atomSel02.indices()[0]
        s.ida2 = atomSela2.indices()[0]
        s.idb2 = atomSelb2.indices()[0]
        s.id03 = atomSel03.indices()[0]
        s.ida3 = atomSela3.indices()[0]
        s.idb3 = atomSelb3.indices()[0]
        s.id04 = atomSel04.indices()[0]
        s.ida4 = atomSela4.indices()[0]
        s.idb4 = atomSelb4.indices()[0]
        if sim==0:
            sim = simulation.currentSimulation()
            pass
        s.sim = sim
        
        s.scale = 1
        return
    def calcEnergy(s):
        q01 = s.sim.atomPos(s.id01).tuple()
        qa1 = s.sim.atomPos(s.ida1).tuple()
        qb1 = s.sim.atomPos(s.idb1).tuple()
        q02 = s.sim.atomPos(s.id02).tuple()
        qa2 = s.sim.atomPos(s.ida2).tuple()
        qb2 = s.sim.atomPos(s.idb2).tuple()
        q03 = s.sim.atomPos(s.id03).tuple()
        qa3 = s.sim.atomPos(s.ida3).tuple()
        qb3 = s.sim.atomPos(s.idb3).tuple()
        q04 = s.sim.atomPos(s.id04).tuple()
        qa4 = s.sim.atomPos(s.ida4).tuple()
        qb4 = s.sim.atomPos(s.idb4).tuple()

        s.va1 = unitVec(subVec(qa1,q01))
        s.vb1 = unitVec(subVec(qb1,q01))
        s.na1 = norm(subVec(qa1,q01))
        s.nb1 = norm(subVec(qb1,q01))
        s.cos1 = dot(s.va1,s.vb1)
        
        s.va2 = unitVec(subVec(qa2,q02))
        s.vb2 = unitVec(subVec(qb2,q02))
        s.na2 = norm(subVec(qa2,q02))
        s.nb2 = norm(subVec(qb2,q02))
        s.cos2 = dot(s.va2,s.vb2)

        s.va3 = unitVec(subVec(qa3,q03))
        s.vb3 = unitVec(subVec(qb3,q03))
        s.na3 = norm(subVec(qa3,q03))
        s.nb3 = norm(subVec(qb3,q03))
        s.cos3 = dot(s.va3,s.vb3)

        s.va4 = unitVec(subVec(qa4,q04))
        s.vb4 = unitVec(subVec(qb4,q04))
        s.na4 = norm(subVec(qa4,q04))
        s.nb4 = norm(subVec(qb4,q04))
        s.cos4 = dot(s.va4,s.vb4)

        return 0.5 * s.scale * (s.cos1*s.cos4 - s.cos2*s.cos3)**2
    def calcEnergyAndDerivs(s,derivs):
        energy = s.calcEnergy()

        pref = s.scale * (s.cos1*s.cos4 - s.cos2*s.cos3) * s.cos4

        da1 = vecScale(pref/s.na1, subVec(s.vb1 ,
                                          vecScale(dot(s.va1,s.vb1),s.va1)))
        db1 = vecScale(pref/s.nb1, subVec(s.va1 ,
                                          vecScale(dot(s.va1,s.vb1),s.vb1)))
        derivs[s.ida1] = addVec(derivs[s.ida1], da1)
        derivs[s.idb1] = addVec(derivs[s.idb1], db1)
        derivs[s.id01] = subVec(derivs[s.id01], addVec(da1,db1))
        
        pref = -s.scale * (s.cos1*s.cos4 - s.cos2*s.cos3) * s.cos3

        da2 = vecScale(pref/s.na2, subVec(s.vb2 ,
                                          vecScale(dot(s.va2,s.vb2),s.va2)))
        db2 = vecScale(pref/s.nb2, subVec(s.va2 ,
                                          vecScale(dot(s.va2,s.vb2),s.vb2)))
        derivs[s.ida2] = addVec(derivs[s.ida2], da2)
        derivs[s.idb2] = addVec(derivs[s.idb2], db2)
        derivs[s.id02] = subVec(derivs[s.id02], addVec(da2,db2))
        
        pref = s.scale * (s.cos1*s.cos4 - s.cos2*s.cos3) * s.cos1

        da4 = vecScale(pref/s.na4, subVec(s.vb4 ,
                                          vecScale(dot(s.va4,s.vb4),s.va4)))
        db4 = vecScale(pref/s.nb4, subVec(s.va4 ,
                                          vecScale(dot(s.va4,s.vb4),s.vb4)))
        derivs[s.ida4] = addVec(derivs[s.ida4], da4)
        derivs[s.idb4] = addVec(derivs[s.idb4], db4)
        derivs[s.id04] = subVec(derivs[s.id04], addVec(da4,db4))
        
        pref = -s.scale * (s.cos1*s.cos4 - s.cos2*s.cos3) * s.cos2

        da3 = vecScale(pref/s.na3, subVec(s.vb3 ,
                                          vecScale(dot(s.va3,s.vb3),s.va3)))
        db3 = vecScale(pref/s.nb3, subVec(s.va3 ,
                                          vecScale(dot(s.va3,s.vb3),s.vb3)))
        derivs[s.ida3] = addVec(derivs[s.ida3], da3)
        derivs[s.idb3] = addVec(derivs[s.idb3], db3)
        derivs[s.id03] = subVec(derivs[s.id03], addVec(da3,db3))
        
        return energy
    pass
