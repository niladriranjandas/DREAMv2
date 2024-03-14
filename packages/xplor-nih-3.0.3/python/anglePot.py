            
from pyPot import PyPot
class AnglePot(PyPot):
    def __init__(s,name,
                 angleList):
        """
        Create a new AnglePot term, given the name and a list of atom triples.
        Create a new AnglePot potential with the given instance name and
        an optional simulation argument.
        """
        PyPot.__init__(s,name)
        s.resetPotName('AnglePot')

        s.restraints_=[]
        for (a,b,c,forcec,theta0) in angleList:
            s.restraints_.append(
                Restraint(a,b,c,forcec,theta0,"",s) )
    
        s.setThreshold(5)  # default threshold in degrees from potential minima
        s.simulation=None
        return
    def restraints(s):
        return s.restraints_
    def numRestraints(s):
        return len(s.restraints_)
    def violations(s):
        ret=0
        for r in s.restraints():
            if r.violated():
                ret += 1
                pass
            pass
        return ret

    def rms(s):
        """RMSD in degrees
        """
        s.calcEnergy()
        ret=0.
        for r in s.restraints():
            ret += r.diff()**2
            pass
        from math import sqrt
        if s.numRestraints()>0:
            ret = sqrt(ret/s.numRestraints())
            pass
        return ret
        

    def calcEnergy(s):
        if not s.simulation:
            s.simulation=s.restraints()[0].atom1.simulation()
            pass
        ret=0.
        for r in s.restraints():
            ret += r.energy()
            pass
        ret *= s.scale()
        return ret

    def calcEnergyAndDerivList(s,derivs):
        ret=s.calcEnergy()
        from derivList import DerivList
        dlist = DerivList()
        dlist.init()
        
        for r in s.restraints():
            r.deriv(dlist)
            pass
        dlist[s.simulation].scale(s.scale())
        derivs[s.simulation] += dlist[s.simulation]
        return ret

    def showRestraints(s):
        ret  = "\n     Angle Restraints: term: %s\n\n" % s.instanceName()
        ret += "Index     %20s %20s %20s" % ("Atom 1","Atom 2","Atom 3")
        ret += "      Calcd    Target    Diff\n"
        ret += "-"*79
        ret += '\n'
        cnt=0
        for r in s.restraints():
            ret += "%3d "%cnt
            cnt += 1
            ret += "(%20s) (%20s) (%20s) " % (r.atom1.string(),r.atom2.string(),
                                              r.atom3.string())
            ret +=  "%7.2f %7.2f %7.2f\n" %(RAD2DEG*r.theta,
                                            RAD2DEG*r.target,
                                            r.diff())
            pass
        return ret

    pass

from math import pi
DEG2RAD = pi/180.
RAD2DEG = 180./pi
from vec3 import Vec3 #FIX: why is this hack needed?
from vec3 import dot, cross, norm, unitVec
from bondAngle import BondAngle
import protocol

class Restraint:
    def __init__(s,atom1,atom2,atom3,forcec,target,comment,pot):
        """
        min1 and min2 are input in degrees
        """
        s.pot=pot
        from bondAngle import BondAngle
        from atomSel import AtomSel
        s.atom1 = atom1
        s.atom2 = atom2
        s.atom3 = atom3
        s.bondAngle = BondAngle(atom1,atom2,atom3)
        s.forcec=forcec
        s.target = target * DEG2RAD
        s.comment=comment
        return
    def name(s):
        ret = "%s %s %s" % (s.atom1.string(),
                            s.atom2.string(),
                            s.atom3.string())
        return ret
    def diff(s):
        """Return difference from target value in degrees.
        """
        return RAD2DEG*(s.theta-s.target)
    def violated(s):
        if abs(s.diff()) > s.pot.threshold():
            return True
        else:
            return False
        pass
    def energy(s):
        s.theta=s.bondAngle.value()
#        print s.name(), s.theta
        return s.forcec * (s.theta - s.target)**2
    def deriv(s,dlist):
        """
        Called after energy.
        """
        d1,d2,d3= s.bondAngle.derivs()
        pref = 2* s.forcec * (s.theta - s.target)
        dlist[s.atom1] += pref * d1
        dlist[s.atom2] += pref * d2
        dlist[s.atom3] += pref * d3
        return
    pass
    
    
            

def analyze(potList):
    "perform analysis of AnglePot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'AnglePot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        print(term.showRestraints())

        pass
    
    return ret
#
# print print analysis
import simulationTools
simulationTools.registerTerm(analyze,
                             "AnglePot Analysis","AnglePot",
r"""
For each term report the calculated, target value and difference in angle.
""")
        
