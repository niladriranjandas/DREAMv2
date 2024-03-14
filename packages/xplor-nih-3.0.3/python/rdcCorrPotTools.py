
def create_RDCCorrPot(name,
                      file=None,
                      restraints="",
                      varTensor=None,
                      simulation=None):
    """
    Construct an RDCCorrPot object from a <m ccrPot>.CCRPot term.
    The arguments are:

       name       - instanceName
       file       - a single restraint table filename, or a sequence of
                    filenames.
       restraints - a string containing restraints
       varTensor  - an optional <m varTensor>.VarTensor object used to define
                    the B-field axis.
       simulation - an optional <m simulation>.Simulation object
       
    Based on
    C. Camilloni and M. Vendruscolo.
    A tensor-free method for the structural and dynamical refinement of proteins
    using residual dipolar couplings.
    J. Phys. Chem. B 9, 653-661 (2015).
    """
    from ccrPot import CCRPot
    from simulation import currentSimulation
    from varTensorTools import create_VarTensor

    if not simulation: simulation = currentSimulation()

    if not varTensor:
        varTensor = create_VarTensor("RDCCorr-default")
        varTensor.setFreedom("fix")
        pass

    pot = CCRPot(name,"",simulation)
    pot.resetPotName("RDCCorrPot")
    pot.setPotType("correlation")

    pot.varTensor = varTensor
    from ccrPot import realCCRPot
    realCCRPot.addXplorRestraints = addXplorRestraints

    if restraints:
        pot.addXplorRestraints(string=restraints)
        pass

    if file:
        files=[]
        if type(file)==type("string"):
            files.append( file)
        else:
            files=file
            pass
        for file in files:
            pot.addXplorRestraints(file)

            pass
        pass

    
    return pot

def addXplorRestraints(self,
                       filename="",
                       string=""):
    """
    Read in Xplor-formatted RDC restraints from the specified filename or
    string.
    """
    restraints = ""
    if filename:
        restraints += open(filename).read()
        pass
    if string:
        restraints += string
        pass
    self.addRestraints( convertFromXplor(restraints,self.varTensor) )
    return

def convertFromXplor(rdcRestraints,
                     varTensor    ):
    """
    Convert XDIPO RDC restraints to CCR restraints using the Z axis from the
    passed varTensor object
    """

    import potUtils
    restraints = potUtils.stripBracketedComments( rdcRestraints )

    import re

    ret=""

    from selectTools import toAtomSelString
    while restraints:
        restraints = restraints.lstrip()
        while re.match("!.*",restraints):
            restraints = re.sub(r'![^\n]*','',restraints,1)
            restraints = restraints.lstrip()
            pass
            
        m=re.match("assi[a-z]*[^(]*",restraints,re.IGNORECASE)
        if m:
            (restraint,
             restraints) = convertFromXplorOne(restraints[len(m.group()):],
                                          toAtomSelString(varTensor.oAtom()),
                                          toAtomSelString(varTensor.zAtom()))
            ret+= 'assign ' + restraint + '\n'
            pass
        else:
            restraints = re.sub(r'\S*','',restraints,1)
        pass
    return ret


def convertFromXplorOne(string,
                        oAtomSel,
                        zAtomSel):
    """
    convert one XPLOR restraint to CCR.
    """
    from parseTools import readFloat, readInt, findNested
    
    isXplorRestraint=False
    string = string.strip()
    i = findNested('(',')',0,string,0)
    sel1,string = string[1:i], string[i+1:]
    string = string.strip()
    i = findNested('(',')',0,string,0)
    sel2,string = string[1:i], string[i+1:]
    string = string.strip()
    if string[0]=='(':
        isXplorRestraint=True
        # read dummy sels
        i = findNested('(',')',0,string,0)
        sel2,string = string[1:i], string[i+1:]
        string = string.strip()
        i = findNested('(',')',0,string,0)
        sel2,string = string[1:i], string[i+1:]
        string = string.strip()

        i = findNested('(',')',0,string,0)
        sel1,string = string[1:i], string[i+1:]
        string = string.strip()
        i = findNested('(',')',0,string,0)
        sel2,string = string[1:i], string[i+1:]
        string = string.strip()
    

        pass

    (obs,string) = readFloat(string)
    (err,string) = readFloat(string)

    try:
        (err2,string) = readFloat(string)
    except ValueError:
        err2 = err
        pass

    if isXplorRestraint:
        #deduce bond length from sel1, sel2
        from atomSel import AtomSel
        atom1=AtomSel(sel1)[0]
        atom2=AtomSel(sel2)[0]
        idName = atom1.atomName()+"_"+atom2.atomName()
        pass
    else:
        pass

    return ("(%s) (%s) (%s) (%s) %f %f" %(sel1,sel2,oAtomSel,zAtomSel,
                                          obs,err),
            string)

def calibrate(term,
              selection="all",
              verbose=False):
    """
    Given a RDCCorrPot term, compute and apply the prefactor
    which optimizes the rmsd fit of observed to calculated RDC. A subset
    restraints can be used in calibrating the prefactor by specifying
    a restricting selection argument.
    """
    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection)

    from cdsMatrix import CDSMatrix_double, inverse
    from cdsVector import CDSVector_double, sum

    rms0 = term.rms()
    term.setPrefactor(1.)

    calcdVals=[]
    targetVals=[]
    from atomSel import intersection
    for r in term.restraints():
        if (selection.containsAtom(r.atomX1()) and
            selection.containsAtom(r.atomX2())    ):
            calcdVals.append(r.calcd())
            targetVals.append(r.obs())
            pass
        pass
    if verbose:
        print("calibrate: using %d (of %d) restraints." % (len(calcdVals),
                                                           term.numRestraints()))
        pass

    calcdVals=CDSVector_double(calcdVals)
    targetVals=CDSVector_double(targetVals)

    n = sum( calcdVals*targetVals)
    d = sum( calcdVals**2)

    if abs(d) < 1e-10:
        print("calibrate: zero denominator. No calibration performed")
        return
    
    p = n / d

    if verbose:
        print("\tprefactor: %.5e" % p)
        #print "\tintercept: %.5e" % intercept
        print("\trmsd:      %.5e --> " %rms0, end=' ')
        pass

    term.setPrefactor(p)
    #term.setIntercept(intercept)
    if verbose:
        print("%.5e" % term.rms())
        pass
    return

    

def Rfactor(pot,
            selection='all'):
    """R-factor (in percent) for an infinite number of
    randomly distributed vectors.
    Eq. 3 from Clore+Garrett, JACS 121, 9008 (1999).

    The selection argument can be used to choose a subset of
    restraints whose atoms lie in selection.
     """
    
    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection,pot.simulation().member())
    import atomSel

    diff2=0
    obs2=0
    for r in pot.restraints():
        if (selection.containsAtom( r.atomX1() ) or
            selection.containsAtom( r.atomX2() ) or
            selection.containsAtom( r.atomY1() ) or
            selection.containsAtom( r.atomY2() )):
            diff2 += (r.calcd() - r.obs())**2
            obs2 +=  r.obs()**2
            pass
        pass
    
    from math import sqrt

    if obs2>0.:
        Rinf = 100 * sqrt(diff2) / sqrt(2*obs2)
    else:
        Rinf=-1
        pass
        
    return Rinf

def correlation(pot):
    """
    Return the RDCCorrPot terms correlation
    """
    return pot.correlation()

    

def analyze(potList):
    "perform analysis of rdcCorrPot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'RDCCorrPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret+= "%9s  %6s  %7s %7s \n" % \
          (" " , "RMS", "corr", "R-factor")

    for name in instanceNames:
        rdc = [x for x in potList if x.instanceName()==name][0]

        print(rdc.showViolations())

        print(rdc.info())

        ret += "%9s  %6.3f  %7.4f %7.3f \n" % \
               (name , rdc.rms(), rdc.correlation(), Rfactor(rdc))
        pass
    
    return ret

import saTensorTools

from simulationTools import registerTerm, registerExtraStats
registerTerm(analyze,"Correlation Dipolar Coupling Analysis","RDCCorr",
r"""
For each term the following are reported::

  RMS      - root mean square deviation between calculated and target sPRE
             values.
  corr     - correlation between calculated and target sPRE values
  R-factor - R-factor quality metric, Eq. 3 from Clore+Garrett, JACS 121,
             9008 (1999).
             
""")

registerExtraStats("RDCCorrPot","correlation",correlation)
registerExtraStats("RDCCorrPot","R-factor",Rfactor)
