"""
Helper functions for <m distSplinePot>, splined energy function of distance
between two atoms. 
"""

from __future__ import print_function

def create_DistSplinePot(instanceName,
                         fileOrData,
                         rSwitch=1,
                         simulation=None):
    import os
    if os.path.exists(fileOrData):
        lines = open(fileOrData).readlines()
    else:
        lines = fileOrData.split('\n')
        pass
    data = []
    for strippedLine in [line.strip() for line in lines]:
        if not strippedLine or strippedLine.startswith('#'):
            continue
        data.append( [float(v) for v in
                      strippedLine.split()[:2]] )
        pass

    from distSplinePot import DistSplinePot
    pot = DistSplinePot(instanceName,
                        [t[0] for t in data],
                        [t[1] for t in data],
                        rSwitch,
                        simulation)
    return pot

    

class Restraint: # XX Fix
    def __init__(s,pot,atom0,atom1,comment=''):
        s.pot=pot

        from selectTools import convertToAtom
        s.atom0 = convertToAtom( atom0,sim=pot.simulation )
        s.atom1 = convertToAtom( atom1,sim=pot.simulation )
        s.comment=comment
        return
    def name(s):
        ret = "%s %s" % (s.atom0.string(),
                         s.atom1.string())
        return ret
#    def diff(s):
#        from math import sqrt
#        if s.dist2<s.sigma2:
#            return sqrt(s.dist2) - sqrt(s.sigma2)
#        else:
#            return 0.
#    def violated(s):
#        if abs(s.diff()) > s.pot.threshold():
#            return True
#        else:
#            return False
#        pass
    def energy(s):
        from vec3 import norm
        s.dist = norm( s.atom0.pos() - s.atom1.pos() )
        return s.pot.eval( s.dist )
    def deriv(s,dlist):
        """
        Called after energy.
        """
        dE_dr = s.pot.dE_dr(s.dist)
        from vec3 import unitVec
        dE_dq = dE_dr * unitVec( s.atom0.pos() - s.atom1.pos() )

        dlist[s.atom0] += dE_dq
        dlist[s.atom1] -= dE_dq
        
        return
    pass

    
            

def analyze(potList): # XX Fix
    "perform analysis of DistPlotPot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'DistSplinePot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]


        print(term.showViolations())
        print(term.info())

        pass
    
    return ret
#
# print print analysis
import simulationTools
simulationTools.registerTerm(analyze,
                             "DistSplinePot Analysis","DistSplinePot",
r"""
This term generates no PDB header information.
""")

        
def plotCurve(pot,
              rmin=None,
              rmax=None,
              numPoints=100,
              ):
    """
    Make a plot of energy vs r between rmin and rmax with the specified number
    of points.
    """

    import matplotlib
    import pylab

    rmin0=pot.rValues()[0]
    rmax0=pot.rValues()[-1]

    if rmin==None: rmin = rmin0
    if rmax==None: rmax = rmax0
    delta_x=(rmax-rmin)/(numPoints-1)

    xvals=[ rmin+i*delta_x for i in range(numPoints)]
    yvals=[ pot.energy(x) for x in xvals ]

    pylab.plot(xvals,yvals)

    ymin = min(yvals)
    ymax = max(yvals)

    pylab.scatter(pot.rValues(),
                  pot.eValues())

    if rmin<rmin0:
        pylab.plot([rmin0,rmin0],[ymin,ymax],color='k')
        pass
    if rmax>rmax0:
        pylab.plot([rmax0,rmax0],[ymin,ymax],color='k')
        pass

    pylab.xlabel("r($\AA$)")
    pylab.ylabel("E(kcal/mol)")

    

    pylab.show()
    return
