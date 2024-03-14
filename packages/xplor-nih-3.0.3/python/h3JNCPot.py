
#create Pot derived from PyPot
from pyPot import PyPot
class H3JNCPot(PyPot):
    def __init__(s,name,
                 simulation=None):
        """
        Create a new H3JNC potential with the given instance name and
        an optional simulation argument.
        """
        PyPot.__init__(s,name)
        s.resetPotName('H3JNCPot')
        s.simulation=simulation
        if s.simulation==None:
            import simulation
            s.simulation=simulation.currentSimulation()
            pass

        s.restraints_=[]
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
    def readRestraints(s,restraintString):
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

    def readRestraint(s,line):
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

    def rms(s):
        s.calcEnergy()
        ret=0.
        for r in s.restraints():
            ret += (r.calcd-r.obs)**2
            pass
        from math import sqrt
        if s.numRestraints()>0:
            ret = sqrt(ret/s.numRestraints())
            pass
        return ret
        

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

    def showRestraints(s):
        ret  = "\n     H3JNC Restraints: term: %s\n\n" % s.instanceName()
        ret += "Index     Proton                Oxygen       "
        ret += "      Calcd    Obs    Diff    Err\n"
        ret += "-"*79
        ret += '\n'
        cnt=0
        for r in s.restraints():
            ret += "%3d "%cnt
            cnt += 1
            ret += "(%s) (%s) " % (r.h3JNC.h.string(),
                                   r.h3JNC.o.string())
            ret +=  "%7.2f %7.2f %7.2f %7.2f\n" %(r.calcd,r.obs,
                                                  r.diff(),r.err)
            pass
        return ret

    pass

class Restraint:
    def __init__(s,sel1,sel2,obs,err,comment,pot):
        s.pot=pot
        from h3JNC import H3JNC
        s.h3JNC = H3JNC(sel1,sel2)
        s.obs=obs
        s.err=err
        s.comment=comment
        return
    def name(s):
        ret = s.h3JNC.h.string() + ' ' + s.h3JNC.o.string()
        return ret
    def diff(s):
        return s.calcd - s.obs
    def violated(s):
        if abs(s.diff()) > s.pot.threshold():
            return True
        else:
            return False
        pass
    def calc(s):
        return s.h3JNC.calc()
    def energy(s):
        s.calcd=s.calc()
        return (s.calcd - s.obs)**2 / s.err**2
    def deriv(s,dlist):
        """
        Called after energy.
        """
        dh,do,dc,dn = s.h3JNC.derivs()
        pref = 2 * (s.calcd - s.obs) / s.err**2
        dlist[s.h3JNC.h] += pref * dh
        dlist[s.h3JNC.o] += pref * do
        dlist[s.h3JNC.c] += pref * dc
        dlist[s.h3JNC.n] += pref * dn
        return
    pass
    
    
            

def analyze(potList):
    "perform analysis of H3JNCPOt terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'H3JNCPot')

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
                             "H3JNC Analysis","H3JNCPot",
r"""
This term should generate no PDB header information.
""")
        
