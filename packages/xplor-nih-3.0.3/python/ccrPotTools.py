# -*- coding: utf-8 -*-

"""Tools to aid in setup/analysis of the Cross-Correlation Relaxation
   potential term.

This module provides functions to simplify the creation and
analysis of <m ccrPot>.CCRPot potential terms.
"""

def create_CCRPot(name,file=0,
                  tauc=None,
                  prefactor=None,
                  sim=None,
                  restraints=""):
    """create an <m ccrPot>.CCRPot term with given name. The CCR prefactor
    can be manually specified, or automatically calculated is tauc is
    specified.

    If the file argument is specified, restraints are read from the given
    filename.
    sim is an optional Simulation argument- the current default is used if
    it is not specified.
    """
    from ccrPot import CCRPot
    from simulation import currentSimulation
    if not sim: sim = currentSimulation()

    if file: restraints += open(file).read()

    pot = CCRPot(name,restraints,sim)

    if prefactor==None:
        if tauc==None:
            raise Exception("either prefactor or tauc must be specified")
        r=pot.restraints()[0]
        prefactor = calc_prefac(r.atomX1().atomName(),
                                r.atomX2().atomName(),
                                r.atomY1().atomName(),
                                r.atomY2().atomName(),tauc)
        pass
                                
        
    pot.setPrefactor( prefactor )
        
    # default threshold value
    pot.setThreshold( 0.2 )

    return pot

def calc_prefac(a1,a2,a3,a4,tauC):
    """
    Calculate the CCR prefactor based on atom names and the specified tauc.

    FIX: need to specify units.
    """
    mu     =  1.2566371
    hbar   =  1.054571
    gammaH = 42.576
    gammaN = -4.3156
    gammaC = 10.705
    pi     =  3.141592653 
    rnh      =  1.02
    rcaha    =  1.09
    rcac     =  1.53
    rcbcg    =  1.53
  
    atoms = [a1,a2,a3,a4]
    gamma = []
    for atom in atoms:
      if atom[0] == "N":
        gamma.append(gammaN)
      elif atom[0] == "H":
        gamma.append(gammaH)
      elif atom[0] == "C":
        gamma.append(gammaC)
      else:
        raise Exception("atom gamma not defined: %s"% atom)
        
    pair1 = [a1,a2]
    pair2 = [a3,a4]
    
    dists = []
    
    pairs = [pair1,pair2]
    for pair in pairs:
      if (("N" in pair) and ("HN" in pair)):
        dists.append(rnh)
      elif (("CA" in pair) and ("HA" in pair)):
        dists.append(rcaha)
      elif (("CA" in pair) and ("C" in pair)):
        dists.append(rcac)
      elif (("CB" in pair) and ("CG" in pair)):
        dists.append(rcbcg)
      else:
          raise Exception("distance not defined: %s" %str(pair))
          
    # calculate prefactor
    part1a = ((( mu * hbar * gamma[0] * gamma[1] * pi )/( dists[0] * dists[0] * dists[0] )) * 1e2 ) 
    part1b = ((( mu * hbar * gamma[2] * gamma[3] * pi )/( dists[1] * dists[1] * dists[1] )) * 1e2 )
    part2 = (( 2.0 * tauC ) / 5.0 )
    prefactor = (( part1a * part1b ) * part2 )
    
   # print "prefactor is : ", prefactor, " using Tc : ", tauC
  
    return prefactor

def analyze(potList):
    "perform analysis of CCRPot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'CCRPot')

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
registerTerm(analyze,"CCR terms","CCR",
"""
For each CCR term, print the root mean square deviation between calculated
and observed, the deviation between ensemble members, and the number of
violations.

Please see:

    R. Fenwick, C.D. Schwieters and B. Vögeli, Direct
    Investigation of Slow Correlated Dynamics in Proteins via
    Dipolar Interactions, J. Am. Chem. Soc. 138, 8412-8421
    (2016).
""")
