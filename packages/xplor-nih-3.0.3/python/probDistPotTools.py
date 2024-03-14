"""
 tool for creating a probDistPot
"""
def create_probDistPot(name,
                       targetMap,
                       sel,
                       potType="diff",
                       scale=1):
    """
    creates a <m probDistPot>.ProbDistPot with given name

    targetMap contains the target Map which is an DenstiyGrid object

    potType: switches the energy type.
       "diff" which is the default energy function, is based on the squared
              distance between the target density & current density calculated
       "cross_correlation" is based on a cross correlation
       "correlation" is the correct correlation function
               
    scale: set the scale for the term
    
    """
   
    #from potList import PotList
    #pl = PotList(name)
    
    from probDistPot import ProbDistPot
    probPot=ProbDistPot(name,targetMap,sel)
    probPot.setPotType(potType)
    probPot.setScale(scale)

    #pl.append(probPot)
    return probPot # removed potList to return violation as -1
    
    
def analyze(potList):
    """Perform analysis of ProbDistPot terms and return nicely formatted
    summary.
    """

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'ProbDistPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret+= "%9s  Correlation\n" % " "

    for name in instanceNames:
        p = [x for x in potList if x.instanceName()==name][0]

        print(p.info())

        ret += \
            "%-9s  %6.3f \n" % \
            (name , p.correlation(),)
        pass
    
    return ret

def correlation(term):
    return term.correlation()

from simulationTools import registerTerm, registerExtraStats
registerTerm(analyze,"Probability Distribution Potential","ProbDistPot",
r"""
for each term, report the correlation between calculated and target maps.
""")

registerExtraStats("ProbDistPot","corr",correlation)
