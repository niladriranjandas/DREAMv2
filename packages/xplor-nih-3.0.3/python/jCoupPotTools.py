"""tools to aid in setup/analysis of dipolar coupling potential term.

this module provides functions to simplify the creation and
analysis of <m jCoupPot>.JCoupPot potential terms.
"""

def create_JCoupPot(name,file=0,
                    A=0,B=0,C=0,phase=-60,
                    sim=0,restraints=""):
    """create an <m jCoupPot>.JCoupPot term with given name,
    and Karplus relationship parameters, A, B, C, and phase.
    If the file argument is specified, restraints are read from the given
    filename.
    sim is an optional Simulation argument- the current default is used if
    it is not specified.
    """
    from jCoupPot import JCoupPot
    from simulation import currentSimulation


    if file: restraints += open(file).read()

    #FIX this when JCoupPot inherits from EnsemblePot
    if sim:
        pot = JCoupPot(name,restraints,sim)
    else:
        pot = JCoupPot(name,restraints)

    if A: pot.setA(A)
    if B: pot.setB(B)
    if C: pot.setC(C)
    if phase: pot.setPhase(phase)

    if pot.A()==0. and pot.B()==0. and pot.C()==0.:
        raise "createJCoupPot: please define the Karplus relation coefficients"
    
    # default threshold value
    pot.setThreshold( 1.0 )

    return pot

def analyze(potList):
    "perform analysis of JCoupPot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'JCoupPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret += "%-9s  %6s  %6s  %6s\n" % \
           ( "", "RMS", "Devia", "Viols")

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        print(term.showViolations())
        print(term.info())



        ret += "%-9s  %6.3f  %6.3f  %6d\n" % \
               (name , term.rms(), term.deviation(), term.violations() )
        pass
    
    return ret

from simulationTools import registerTerm
registerTerm(analyze,"J-Coupling terms","JCoup",
r"""
For each term report::

  RMS           - root mean square deviation between calculated and target
                  RDC values.
  Devia         - the deviation in calculated J-coupling values of different
                  ensemble members.
  Viols         - the number of violated restraints.
  
""")
