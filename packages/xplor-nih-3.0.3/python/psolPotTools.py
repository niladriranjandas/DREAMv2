"""Tools to aid in setup/analysis of the Solvent Paramagnetic Relaxation
Enhancement potential term.

This module provides functions to simplify the creation, manipulation and
analysis of <m psolPot>.PSolPot potential terms.



"""


def create_PSolPot(name,
                   file=None,
                   restraints="",
                   tauc=0.2,
                   probeR=5.4,
                   probeC=4.0, #units?
                   fixTauc=True,
                   eSpinQuantumNumber=3.5,
                   domainSelection="known and not pseudo",
                   ):
    """Create an <m psolPot>.PSolPot with given name, given the filename of an
    PRE assignment table and/or a string of assignments.

      tauc               - correlation time
      probeR             - radius of probe molecule
      probeC             - probe concentration - units?
      fixTauc            - whether to fix the value of tauc, or to let it float
      eSpinQuantumNumber - electron sping quantum number
      domainSelection    - atoms to use in surface calculation
    """
    from selectTools import convertToAtomSel
    domainSelection = convertToAtomSel(domainSelection)

    #set up radii
    radii=[]
    for atom in domainSelection:
        if     atom.atomName()[0] == 'H': radius = 1.00
        elif   atom.atomName()[0] == 'C': radius = 1.85;
        elif   atom.atomName()[0] == 'N': radius = 1.75;
        elif   atom.atomName()[0] == 'O': radius = 1.60;
        elif   atom.atomName()[0] == 'S': radius = 2.00;
        elif   atom.atomName()[0] == 'P': radius = 2.00;
        elif   atom.atomName()[:2].upper() == 'MN': radius = 1.30
        else:
            radius=1.7
            print("WARNING: Atom " + atom.string(), end=' ')
            print(" is not supported by psolPot. radius=1.7 is used for this.")
            pass
        radii.append(radius)
        pass
    

    import psolPot
    psol = psolPot.PSolPot(name,domainSelection,radii)
    if file != None:
        psol.addRestraints( open(file).read() )
        pass
    psol.addRestraints( restraints )
    
    #psol.setDp(0.8)  #FIX
    psol.setFixTauc(fixTauc)
    psol.setTauC( tauc )
    psol.setSqn(eSpinQuantumNumber)      # electron spin quantum number
    psol.setRho0( probeC )
    
    #psol.setVerbose(True)
    #psol.tessellation().setVerbose(True)
    psol.tessellation().setMoveTol(0.7)
    #psol.tessellation().setSavePrev(True)
    psol.setProbeRadius( probeR )

    return psol
    

def analyze(potList):
    """Perform analysis of PSolPot terms and return nicely formatted summary."""

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'PSol')

    if not potList: return ret

    #set to get rid of duplicates. list for ordering.
    instanceNames = list(set([x.instanceName() for x in potList]))
    instanceNames.sort()

    ret+= "%9s  %6s  %6s  %7s %7s  %4s\n" % \
          (" " , "RMS", "corr.", "Q-fac", "tau_c", "Viols")

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        if term.targetType()=="correlation":
            fitRho0(term)
            pass

        print(term.showViolations())

        print(term.info())

        ret += "%-9s  %6.3f  %6.3f  %7.3f %7.3f  %4d\n" % \
               (name , term.rms(), term.correlation(),
                term.qFactor(),
                term.tauc(), term.violations())
        pass
    
    return ret

def fitRho0(term,weighted=False):
    """For the given <m psolPot>.PsolPot term, set rho0 such that term.rms()
    is optimal. If weighted=True, the computed rho0 instead optimizes a
    weighted rmsd:

      rmsd_w = sum_i w_i (Gamma_i - Gamma_i^obs)^2

    where w_i = 1 / err_i^2. With weighted=True, restraints with err=0 are
    excluded from this calculation.
    
    """

    term.setRho0(1)
    term.calcEnergy() # should be unnecessary
    num=0
    denom=0
    for r in term.restraints():
        w=1.
        if weighted:
            w = r.err()**-2 if r.err()>0. else 0.
            pass
        num += w * r.obs() * r.calcd()
        denom += w * r.calcd()**2
        pass
    rho0 = num / denom
    term.setRho0( rho0 )
    return


import simulationTools
simulationTools.registerTerm(analyze,
                             "Solvent PRE Analysis",
                             "PSol",
r"""
For each term the following are reported:

RMS   - root mean square deviation between calculated and target sPRE values
corr. - correlation between calculated and target sPRE values
Q-fac - Q-factor quality metric
tau_c - fit overall correlation time
Viols - number of violated restraints
""")

def tauc(term):
    if term.potName()=='PSol':
        return term.tauc()
    else:
        # assume we have a list of terms
        from simulationTools import flattenPotList
        prelist = flattenPotList(term)
        if len(prelist)==0:
            return 0.;
        
        from math import sqrt

        sum=0.
        for p in prelist:
            sum += p.tauc()
            pass
        return sum/len(prelist)
    return 

def qFactor(term):
    if term.potName()=='PSol':
        return term.qFactor()
    else:
        # assume we have a list of terms
        from simulationTools import flattenPotList
        prelist = flattenPotList(term)
        if len(prelist)==0:
            return 0.;
        
        from math import sqrt

        sum=0.
        for p in prelist:
            sum += p.qFactor()
            pass
        return sum/len(prelist)
    return 

def correlation(term):
    if term.potName()=='PSol':
        return term.correlation()
    else:
        # assume we have a list of terms
        from simulationTools import flattenPotList
        prelist = flattenPotList(term)
        if len(prelist)==0:
            return 0.;
        
        from math import sqrt

        sum=0.
        for p in prelist:
            sum += p.correlation()
            pass
        return sum/len(prelist)
    return 

simulationTools.registerExtraStats("PSol","correlation",correlation,True)
simulationTools.registerExtraStats("PSol","Q-fac",qFactor,True)
simulationTools.registerExtraStats("PSol","tauc",tauc,True)
