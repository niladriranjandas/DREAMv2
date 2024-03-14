"""tools to aid in setup/analysis of <m orderPot>.OrderPot potential terms

"""


def create_OrderPot(name,file=0,restraints="",simulation=0):
    """
    create a <m orderPot>.OrderPot term appropriate for refining against
    S^2 order parameters. This potential term is only meaningful in the context
    of <m ensembleSimulation> calculations.
    """

    from orderPot import OrderPot

    if file: restraints += open(file).read()

    if simulation:
        pot = OrderPot(name,restraints,simulation)
    else:
        pot = OrderPot(name,restraints)
        pass

    pot.setThreshold(0.)

    return pot

def analyze(potList):
    "perform analysis of OrderPot terms and return a nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'OrderPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret += "%-9s  %6s  %6s  %6s  %6s\n" % \
           ( "", "RMS", "ave", "min", "max")

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        print(term.showViolations())
        print(term.info())

        from cdsVector import CDSVector_double, sum
        values = CDSVector_double([r.calcd() for r in term.restraints()])

        ave = sum(values) / len(values)

        ret += "%-9s  %6.3f  %6.3f  %6.3f  %6.3f\n" %\
               (name , term.rms() ,ave, min(values), max(values))
        pass
    
    return ret

from simulationTools import registerTerm
registerTerm(analyze,"Order Parameter","OrderPot",
r"""
For each term, the following is reported:
 RMS    - root mean square deivation between calculated and observed
 ave    - average calculated order parameter
 min    - minimum calculated order parameter
 max    - maximum calculated order parameter
""")
