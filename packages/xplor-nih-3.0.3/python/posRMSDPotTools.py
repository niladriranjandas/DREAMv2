"""tools to aid in setup/analysis of potential terms restricting atomic
positions in ensemble calculations.

this module provides functions to simplify the creation and
analysis of <m posRMSDPot>.POSRMSDPot potential terms. 
"""


def create_BFactorPot(name,file=0,restraints="",centerSel="",simulation=0):
    """
    create a <m posRMSDPot>.PosRMSDPot term appropriate for refining
    against crystallographic temperature factors. This potential term is
    only meaningful in the context of <m ensembleSimulation> calculations.
    """

    from posRMSDPot import realPosRMSDPot

    if file: restraints += open(file).read()

    if simulation:
        pot = realPosRMSDPot(name,restraints,simulation)
    else:
        pot = realPosRMSDPot(name,restraints)
        pass

    pot.setRMSDType("bfactor")

    if centerSel:
        pot.setCenterSel(centerSel)
    else:
        pot.setCenterSel("known")
        pass

    pot.setThreshold(0.)

    import potProxy
    return potProxy.PotProxy(pot)

def to_RAPPot(ptr):
    """cast an rc_Pot to a RAPPot
     Note: the tol member is lost."""
#    print str(ptr)
#    if str(ptr).potName() != 'RAPPot': raise "pot is not a RAPPot"
    from posRMSDPot import to_PosRMSDPot
    pot = to_PosRMSDPot(ptr)
    pot.resetPotName("RAPPot")
    return pot
    

from posRMSDPot import realPosRMSDPot
class RAPPot(realPosRMSDPot):
    """
    relative atomic position (RAP) potential term based on a
    <m posRMSDPot>.PosRMSDPot. This potential term is
    only meaningful in the context of <m ensembleSimulation> calculations.
    """
    def __init__(s,name,selection,simulation=0):

        from atomSel import AtomSel
        if type(selection)==type(AtomSel("known")):
            selection=selection.string()
            pass

        s.tol_=0.
        s.constructPot(name,selection,simulation)
        return

    def constructPot(s,name,selection,simulation):
        from posRMSDPot import realPosRMSDPot
        if simulation:
            realPosRMSDPot.__init__(s,name,simulation)
        else:
            realPosRMSDPot.__init__(s,name)
            pass
        
        s.resetPotName("RAPPot")

        s.setRMSDType("rap")

        s.addRestraints("assign (%s) 0 %f" % (selection,s.tol_))
        return
    def tol(s):
        return s.tol_
    def setTol(s,val):
        s.tol_ = val
        for r in s.restraints():
            r.obsErr = s.tol_
            pass
        
        return
    pass

realRAPPot = RAPPot
def RAPPot(*args):
    from potProxy import PotProxy
    return PotProxy( realRAPPot(*args) )


#def create_RAPPot(name,selection,simulation=0):
#    """
#    create a relative atomic position (RAP) potential term using a
#    <m posRMSDPot>.PosRMSDPot. This potential term is
#    only meaningful in the context of <m ensembleSimulation> calculations.
#    """
#
#    from posRMSDPot import PosRMSDPot
#
#    if simulation:
#        pot = PosRMSDPot(name,simulation)
#    else:
#        pot = PosRMSDPot(name)
#        pass
#
#    pot.setRMSDType("rap")
#
#    if type(selection)==type(AtomSel("known")):
#        selection=selection.string()
#        pass
#
#    pot.addRestraints("assign (%s) 0 0" % selection)
#        
#
#    return pot

        
def analyze(potList):
    "perform analysis of PosRMSDPot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,('PosRMSDPot','RAPPot'))

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret += "%-9s  %6s\n" % \
           ( "", "RMS")

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        print(term.showViolations())
        print(term.showPositions())
        print(term.info())



        ret += "%-9s  %6.3f\n" % \
               (name , term.rms() )
        pass
    
    return ret

from simulationTools import registerTerm
registerTerm(analyze,"Positional RMSD Analysis","PosRMSD",
r"""
For each term the root mean square difference between calculated and observed
is reported.
""")
