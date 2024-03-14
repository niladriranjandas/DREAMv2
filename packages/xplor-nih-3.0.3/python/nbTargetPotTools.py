
"""
Helper tools for nbTargetPot. Used to fit effective surface
accesibility for solvent NOE and PRE data, as described in

   Yu Wang, C.D. Schwieters, N. Tjandra, ``Parameterization of
   Solvent-Protein Interaction and Its Use on NMR Protein
   Structure Determination,'' J. Magn. Res. 221, 76-84 (2012).
"""

def create_NBTargetPot(name,
                       filename=None,
                       restraints="",
                       moveTol=3,
                       invPow=2,
                       cutoffDist=20,
                       slope=1,
                       intercept=0,
                       restraintFormat="residValErr",
                       restraintSegid=None,
                       restraintAtom="N",
                       obsErr=None,
                       selection="not pseudo and not name H*"):
    """
    Create an NBTargetPot term with the instanceName given by the first argument.
    Restraints can be specified by filename or as a string using the restraints
    argument. Default values of <m nbTargetPot>.NBTargetPot's moveTol, invPow,
    and cutoffDist are specified by arguments of those names. The slope and
    intercept arguments are passed to the NBTargetPot object to be applied to
    the calculated values. The restraintFormat argument can be 'residValErr' or
    'xplor'. For the former, lines in the restraint table should be of the
    form
        resid obsVal osbErr
    while for the xplor restraintFormat, the restraints are given using an
    XPLOR atom selection in the format:
        assign (selection) obsVal obsErr
    The restraintSegid is useful when using residValErr restraints so that a
    segid can be specified. restraintAtom specifies the atom in each
    residue to which distances are to be measured - it is only used for the
    residValErr format. If obsErr is specified, it overrides the values
    in the restraint tables.
    The selection argument species which atoms are included in sums as described
    in the documentation for NBTargetPot.
    """
    from nbTargetPot import NBTargetPot
    pot = NBTargetPot(name,selection)
    pot.setMoveTol(moveTol)
    pot.setInvPow(invPow)
    pot.setCutoffDist(cutoffDist)

    pot.setSlope(slope)
    pot.setIntercept(intercept)

    if filename!=None:
        restraints += open(filename).read()
        pass

    table = readRestraints(restraints,restraintFormat,
                           1.,0.,obsErr,restraintAtom,restraintSegid)

    pot.addRestraints(table)
    if restraints and pot.numRestraints()==0:
        print("create_NBTargetPot: Warning: no restraints read.")
        pass

    return pot
    
def readRestraints(restraintString,
                   format='residValErr',
                   slope=1,
                   intercept=0,
                   obsErr=None,
                   atomName="N",
                   segid=None):
    """
    Read restraints from the specified restraintString. A linear transformation
    is applied to the observed target values using the slope and
    intercept arguments. The format argument can be 'residValErr' or
    'xplor'. For the former, lines in the restraint table should be of the
    form
        resid obsVal osbErr
    while for the xplor format, the restraints are given using an
    XPLOR atom selection in the format:
        assign (selection) obsVal obsErr
    The segid is useful when using residValErr restraints so that a
    segid can be specified. If obsErr is specified, it overrides the values
    in the restraint tables.
    """

    if format=="residValErr":
        restraints = readPlainRestraints(restraintString)
        pass
    elif format=="xplor":
        restraints = readXplorRestraints(restraintString)
        pass
    else:
        raise Exception("invalid restraint format: " + format)

    transRestraints=[]
    for (sel,obs,err) in restraints:
        transRestraints.append( (sel, slope*obs+intercept, slope*err) )
        pass

    table=''
    for (sel,obs,err) in transRestraints:
        if obsErr!=None: err=obsErr
        if format=="residValErr":
            if segid:    sel = 'segid "%s" and ' % segid + sel
            if atomName: sel = sel + ' and name "%s" ' % atomName
            pass
        table += 'assign (%s) %f %f\n' % (sel,obs,err)
        pass
    return table

def readPlainRestraints(restraintString):
    import locale # FIX: should atoi/atof really be locale dependent?
    ret = []
    for line in restraintString.split('\n'):
        if not line:
            continue
        fields=line.split()
        try:
            resid = locale.atoi(fields[0])
            val   = locale.atof(fields[1])
            err   = locale.atof(fields[2])
        except ValueError:
            continue
        ret.append( ("resid %d"%resid,val,err) )
        pass

    return ret


def readXplorRestraints(restraintString):
    # first deal with {} comments
    import potUtils
    restraintString = potUtils.stripBracketedComments( restraintString )

    # one assignment statement per line
    #
    import re
    ret=[]
    for line in restraintString.split('\n'):
        line = line.lstrip()
        m=re.match("assi[a-z]*[^ (]",line,re.IGNORECASE)
        if m:
            ret.append( readRestraint(line[len(m.group()):]) )
            pass
        pass
    return ret

def readRestraint(line):
    from parseTools import findNested, readFloat
    
    end = findNested("(",")",0,line,0)
    end += 1
    sel=line[:end]
    line=line[end:]
    obs,line=readFloat(line)
    err,line=readFloat(line)
    line=line.lstrip()
    comment=''
    if line and line[0]=='!':
        comment=line[1:]
        pass


    return (sel,obs,err)

def calibrate(term,
              selection="all",
              cutoff=None,
              verbose=False):
    """
    Given a <m nbTargetPot>.NBTargetPot term, compute and apply the slope and
    intercept so that the computed values best fit the target (experimental)
    values. The optional selection argument can be used to specify a subset
    of restraints to use in the fit. At this time, errors are not used to
    weight the fit.

    In addition to specifying an atom selection to restrict the restraints used
    in performing the fit, one can specify a cutoff: restraints whose target
    values are greater than this cutoff will not be used in calibration.
    """
    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection)

    from cdsMatrix import CDSMatrix_double, inverse
    from cdsVector import CDSVector_double, sum

    rms0 = term.rms()

    term.setSlope(1.)
    term.setIntercept(0.)

    calcdVals=[]
    targetVals=[]
    from atomSel import intersection
    for r in term.restraints():
        if (len(intersection(selection,r.targetSel))>0 and
            (cutoff==None or r.targetVal<=cutoff)         ):
            calcdVals.append(r.calcdVal)
            targetVals.append(r.targetVal)
            pass
        pass
    if verbose:
        print("calibrate: using %d (of %d) restraints." % (len(calcdVals),
                                                            term.numRestraints()))

    calcdVals=CDSVector_double(calcdVals)
    targetVals=CDSVector_double(targetVals)
    
    M=CDSMatrix_double([ [sum(calcdVals**2), sum(calcdVals)],
                         [sum(calcdVals)   , len(calcdVals)]] )
    b=CDSVector_double( [sum(targetVals*calcdVals), sum(targetVals)] )

    slope,intercept=tuple(inverse(M)*b)

    if verbose:
        print("\tslope:     %.5e" % slope)
        print("\tintercept: %.5e" % intercept)
        print("\trmsd:      %.5e --> " %rms0, end=' ')
        pass

    term.setSlope(slope)
    term.setIntercept(intercept)
    if verbose:
        print("%.5e" % term.rms())
        pass
    return

    

    
def correlation(term):
    return term.correlation()


def analyze(potList):
    "perform analysis of NBTargetPot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'NBTargetPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret += "%-9s  %6s  %6s  %6s %7s\n" % \
           ( "", "RMS", "Devia", "Viols", "Correlation")

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        print(term.showViolations())
        print(term.info())

        ret += "%-9s  %6.3f  %6.3f  %6d  %7.3f\n" % \
               (name , term.rms(), term.deviation(), term.violations(),
                term.correlation())
        pass
    
    return ret

from simulationTools import registerTerm, registerExtraStats
registerTerm(analyze,"NBTargetPot terms","NBTarget",
r"""
For each term, the following are reported:
 RMS         - root mean square deivation between calculated and observed
 Devia       - deviation between ensemble members
 Viols       - number of violated terms
 Correlation - correlation between calculated and observed.
""")
registerExtraStats("NBTargetPot","correlation",correlation)
