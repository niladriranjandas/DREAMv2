"""
Energy term used to add a repulsive energy term to regions of Ramachadran space
where a violation is detected (as informed by TorsionDB).
"""

#create Pot derived from PyPot
from pyPot import PyPot
class RamaPushPot(PyPot):
    def __init__(s,name,
                 selection="not pseudo",
                 width=60,
                 simulation=None):
        """
        Create a new RamaPushPot potential with the given instance name and
        an optional simulation argument.
        """
        PyPot.__init__(s,name)
        s.resetPotName('RamaPushPot')
        from math import pi
        s.width_=width * pi/180. #radians
        s.simulation=simulation
        if s.simulation==None:
            import simulation
            s.simulation=simulation.currentSimulation()
            pass

        s.restraints_=[]
        from torsionDBPotTools import create_TorsionDBPot
        s.tdb = create_TorsionDBPot('tdb',selection=selection,
                                    threshold='99.8',
                                    database='top8K_bb')
        
        s.setThreshold(0.5)  # default threshold for violation count
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
    def clearRestraints(s):
        s.restraints_=[]
        return
    def ageRestraints(s):
        restraints = s.restraints_
        s.restraints_ = []
        for r in restraints:
            if r.lifetime()>0:
                r.age()
                s.restraints_.append(r)
                pass
            pass
        return
    def addRestraint(s,**kwargs):
        s.restraints_.append( Restraint(s,**kwargs) )
        return
                                        
        return
    def readRestraints(s,restraintString):# XX Fix
        """
        Read in restraint table
        """
        # first deal with {} comments
        import potUtils
        restraintString = potUtils.stripBracketedComments( restraintString )

        # one assignment statement per line
        #
        import re
        for line in restraintString.split('\n'):
            line = line.lstrip()
            m=re.match("assi[a-z]*[^ (]",line,re.IGNORECASE)
            if m:
                s.restraints_.append( s.readRestraint(line[len(m.group()):]) )
                pass
            pass
        return

    def readRestraint(s,line): # XX Fix
        from parseTools import findNested, readFloat

        end = findNested("(",")",0,line,0)
        end += 1
        sel1=line[:end]
        line=line[end:]
        end = findNested("(",")",0,line,0)
        end += 1
        sel2=line[:end]
        line=line[end:]
        obs,line=readFloat(line)
        err,line=readFloat(line)
        line=line.lstrip()
        comment=''
        if line and line[0]=='!':
            comment=line[1:]
            pass


        return Restraint(sel1,sel2,obs,err,comment,s)

    def configure(s,verbose=False,lifetime=0):
        """Configure the specified RamaPushPot based on the currently violated
        residues.
        """
        s.ageRestraints()
        for pots in sorted(s.tdb,key=lambda p: p.instanceName()):
            for pot in sorted(pots,key=lambda p: p.instanceName()):
                for r in pot.restraints():
                    if not r.violated():
                        continue
                    phi = r.angle1
                    psi = r.angle2
                    size = r.energy() - pot.threshold()
                    s.addRestraint(phi=phi,psi=psi,
                                   width=s.width_,size=size,
                                   lifetime=lifetime)
                    if verbose:
                        r = s.restraints()[-1]
                        from math import pi
                        phi = r.phi.value() * 180 / pi
                        psi = r.psi.value() * 180 / pi
                        print("RamaPushPot.configure: adding:", end=' ')
                        print(r.name(), end=' ')
                        print("at phi=%7.1f phi=%7.1f" %(phi,psi))
                    pass
                pass
            pass
        return

    def rms(s): 
        return -1

    def calcEnergy(s):
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

    def showRestraints(s):  # XX Fix
        ret  = "\n     RamaPush Restraints: term: %s\n\n" % s.instanceName()
        ret += "Index %14s   phi     phi0    dphi    psi   psi0     dpsi     energy\n" % (
            "Residue",)
        ret += "-"*79
        ret += '\n'
        cnt=0
        from math import pi
        for r in s.restraints():
            ret += "%3d "%cnt
            cnt += 1
            ret += " %14s " % r.name()
            phi  = r.phi.value() * 180./pi
            phi0 = r.phi0        * 180./pi
            psi  = r.psi.value() * 180./pi
            psi0 = r.psi0        * 180./pi
            ret += "%7.2f %7.2f %7.2f " % (phi,phi0,phi-phi0)
            ret += "%7.2f %7.2f %7.2f " % (psi,psi0,psi-psi0)
            ret += "%8.3f\n" % r.energy()
            pass
        return ret

    pass

class Restraint: # XX Fix
    def __init__(s,pot,phi,psi,width,size,comment='',
                 lifetime=0):
        s.pot=pot
        s.phi = phi
        s.psi = psi
        s.phi0 = phi.value()
        s.psi0 = psi.value()
        s.sigma2 = width**2
        s.size = size
        s.comment=comment
        s.lifetime_=lifetime
        return
    def name(s):
        atom = s.phi.atom1()
        ret = "%4s %4d %4s" % (atom.segmentName(),
                               atom.residueNum(),
                               atom.residueName())
        return ret
    def lifetime(s): return s.lifetime_
    def age(s):
        s.lifetime_ -= 1
        return
    def diff(s):
        from math import sqrt
        if s.dist2<s.sigma2:
            return sqrt(s.dist2) - sqrt(s.sigma2)
        else:
            return 0.
    def violated(s):
        if abs(s.diff()) > s.pot.threshold():
            return True
        else:
            return False
        pass
    def energy(s):
        """For now, circular: width phi = width psi

           (phi-s.phi0)**2 + (psi-s.psi0)**2
        """
        s.dist2 = (s.phi.value()-s.phi0)**2 + (s.psi.value()-s.psi0)**2
        if s.dist2>s.sigma2:
            return 0.
        return s.size * (s.dist2/s.sigma2 -1.)**2
    def deriv(s,dlist):
        """
        Called after energy.
        """
        if s.dist2>s.sigma2:
            return
        pref0 = 4*(s.dist2/s.sigma2 -1.) * s.size / s.sigma2
        pref = pref0 * (s.phi.value()-s.phi0) #phi
        d0,d1,d2,d3 = s.phi.derivs()
        dlist[s.phi.atom0()] += pref * d0
        dlist[s.phi.atom1()] += pref * d1
        dlist[s.phi.atom2()] += pref * d2
        dlist[s.phi.atom3()] += pref * d3

        pref = pref0 * (s.psi.value()-s.psi0) #psi
        d0,d1,d2,d3 = s.psi.derivs()
        dlist[s.psi.atom0()] += pref * d0
        dlist[s.psi.atom1()] += pref * d1
        dlist[s.psi.atom2()] += pref * d2
        dlist[s.psi.atom3()] += pref * d3
        return
    pass

    
            

def analyze(potList): # XX Fix
    "perform analysis of RamaPushPot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'RamaPushPot')

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
                             "RamaPush Analysis","RamaPushPot",
r"""
This term generates no PDB header information.
""")
        
