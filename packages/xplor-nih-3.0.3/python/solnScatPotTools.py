"""tools to aid in setup/analysis of potential terms employing solution
scattering data.

this module provides functions to simplify the analysis of
<m solnScatPot>.SolnScatPot potential terms.  
"""


def analyze(potList):
    """perform analysis of SolnScatPot terms and return nicely formatted
    summary.
    A Chi^2 value is calculated based on the assumption that the weights are
    1/sigma_i^2, where sigma_i is the error in the observed value of I_i.

    The deviation value should be nonzero only for EnsembleSimulation
    calculations with Ne>1. It is calculated as

    Dev^2 = 1/Ne 1/Nk \sum_i \sum_j w_j * (I_ij-I_j)^2

    where Nk is the number of datapoints and all values of I are normalized by
    calcdScale, and multiplied by globCorrect.

    """

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'SolnScatPot')
    
    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret += "%-9s  %6s  %6s  %6s\n" % \
           ( "", "RMS", "Devia", "Chi^2")

    from cdsVector import log, norm
    from math import sqrt
    import ensembleSimulation
    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        print("  Analysis of SolnScatPot term: %s\n" % name)

#        print term.showViolations()
#        print term.showPositions()
        print(term.info())

        Icalc = term.splined()
        Iobs  = [x/term.exptScale() for x in term.expt()]
        globCorrect = term.globCorrect

        calcd = Icalc
#        if term.cmpType()=='log': calcd = log(calcd)
        
        qValues = term.qValuesExpt()

        sim = term.simulation()
        sum=0
        from spline import FloatSpline
        from cdsVector import CDSVector_double
        for index in range(sim.size()):
            contrib = term.I_contrib(index) / term.calcdScale() * globCorrect
            from spline import FloatSpline
            contribSpline=FloatSpline(term.qValues(), contrib)
            splinedContrib= [contribSpline(q) for q in qValues]
            sum += sim.members(index).weight() * \
                   norm(calcd-CDSVector_double(splinedContrib))**2
            pass
        dev = sqrt(sum /len(qValues))

        #term.calcGlobCorrect()
        index = term.normalizeIndex()

        extra=""
        if index==-2:
            extra = "  [normalization calcd: %e  expt: %e]" % \
                    (1/term.calcdScale(), 1/term.exptScale())
        elif index>=0:
            extra = "  [normalized to point %d, q=%e]" %(index,qValues[index])
        print("    Calculated I(q)", extra)
        print("-"*70)
        globHeader=""
#       if term.useGlobs():
        globHeader="gCorrect"
        print("%8s %13s %13s %13s %11s %10s" % ("q",
                                                    "I [calc]","I [obs]",
                                                    "Icalc-Iobs","Obs Error",
                                                    globHeader))
        Icalc = term.splined()
        Iobs  = [x/term.exptScale() for x in term.expt()]

        weights=[w*term.exptScale()**2 for w in term.weights()]
        splinedGlobCorrect=FloatSpline(term.qValues(), globCorrect)
        for i in range(len(qValues)):
            diff = Icalc[i] - Iobs[i]
            error = -1
            if weights[i]>0.: error = 1./sqrt(weights[i])
            print("%8.4f %# -13.4g %# -13.4g %# -13.5g %# -11.5g " % \
                  (qValues[i],Icalc[i],Iobs[i], diff,    error), end=' ')
 #          if term.useGlobs():
            print("%# -10.5g" % splinedGlobCorrect(qValues[i]))
#            else:
#                print
            pass
        print()

        chi2=0
        for i in range(len(Icalc)):
            chi2 += weights[i] * (Icalc[i] - Iobs[i])**2
            pass
        # for denom of chi2, only include points w/ nonzero weight
        N = len( [x for x in term.weights() if x!=0] )
        if N>1:
            chi2 /= N-1
            pass
            
        ret += "%-9s  %6.3f  %6.3f  %6.3f\n" % \
               (name , term.rms(), dev, chi2)
        pass
    
    return ret

def readData(experiment,
             preweighted=False,
             minQ=None,
             maxQ=None):
    """

    read SAXS or SANS data from file, string, or sequence. If a
    sequence, each element should be a sequence of 2 or more values;
    the 1st two are q,I, while the 3'rd is an optional weight
    value. If file, it is read and processed as a string. If a string,
    leading # or ! characters start comment lines which are ignored.
    Noncomment lines must contain at least 2 space-separated values
    corresponding to q and I. An optional third value specifies experimental
    error or weight.
    
    The preweighted and maxQ args specify how string (and file) data is
    interpreted. preweighted=True specifies that the 3rd value is a weight,
    while a False value specifies that it is an experimental error (and
    a weight of 1/e^2 is computed). If maxQ is specified, don't read data with
    larger q values. If minQ is specified, don't read data with
    smaller q values.

    Returns a list of (q,I) tuples and a list of weights for each of these.
    """
    if type(experiment)==type("string"):
        try:
            exptContents=open(experiment).read()
        except:
            exptContents=experiment
            pass
        experiment=[]
        for line in exptContents.split('\n'):
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('!'):
                continue
            
            try:
                (q,I) = list(map(float,line.split()[:2]))
                if maxQ!=None and q>maxQ:
                    continue
                if minQ!=None and q<minQ:
                    continue
                #read weights, if present
                weight=1.0
                if preweighted:
                    try:
                        weight = float(line.split()[2])
                    except IndexError:
                        pass
                    pass
                else:
                    try:
                        error = float(line.split()[2])
                        if error>0:
                            weight=error**-2
                        else:
                            weight = 0
                            pass
                        pass
                    except IndexError:
                        pass
                    pass
                if I<0.:
                    I = 0.
                    if debug:
                        print("warning: I(%f) is <0. Clipping to zero." % q)
                        pass
                    pass
                experiment.append( (q,I,weight) )
                pass
            except:
                pass
            pass
        pass

    weights=[]
    exp0 = experiment
    experiment = []
    for e in exp0:
        experiment.append( e[:2] )
        try:
            weight=e[2]
        except IndexError:
            weight=1
            pass
        weights.append(weight)
        pass
    return experiment, weights



        
def analyze_chi2(term):
    """ compute the Chi^2 value for the <m solnScatPot>.SolnScatPot term
    argument.
    """
    return chi2(term.weights(),
                term.splined(),term.expt(),term.normalizeIndex())

def splineCurve(experiment,weights,numPoints):
    """
    given a curve with elements (q,I) and corresponding weights, return the
    corresponding splined curves sampled uniformly over numPoints.
    """
    from spline import FloatSpline
    splinedI=FloatSpline(experiment)
    qWeights = [(experiment[i][0],weights[i]) for i in range(len(experiment))]
    splinedW=FloatSpline(qWeights)
    qMin = min([q_I[0] for q_I in experiment])
    qMax = max([q_I1[0] for q_I1 in experiment])
    qDelta = float(qMax-qMin)/(numPoints-1)
    qValues = [qMin+i*qDelta for i in range(numPoints)]
    
    experiment = [(qValues[i],splinedI(qValues[i])) for i in range(numPoints)]
    weights= [splinedW(qValues[i]) for i in range(numPoints)]
    print("# datapoints sampled uniformly at %d points" % numPoints, end=' ')
    print("using a Spline representation")
    return (experiment,weights)
    pass

class Interp:
    """
    class for linear interpoation (/extrapolation)
    """
    def __init__(s,listOfTuples):
        """
        specify a list of (x,y) pairs sorted such that the x values are
        monotonically increasing.
        """
        x=listOfTuples[0][0]
        for (xp,yp) in listOfTuples[1:]:
            if xp<=x:
                raise Exception("x-values are not monotonically increasing")
            x=xp
            pass
        s.xy=list(listOfTuples)
        return
    def __call__(s,x):
        """
        return the result of interpolation (/extrapolation) to the point x.
        """
        
        #
        # before list - extrapolate
        #
        if x<=s.xy[0][0]:
            retval = s.xy[0][1] + (x-s.xy[0][0]) * \
                     (s.xy[1][1]-s.xy[0][1]) / (s.xy[1][0]-s.xy[0][0])
            return retval

        #
        # after list - extrapolate
        #
        if x>=s.xy[-1][0]:
            retval = s.xy[-1][1] + (x-s.xy[-1][0]) * \
                     (s.xy[-2][1]-s.xy[-1][1]) / (s.xy[-2][0]-s.xy[-1][0])
            return retval

        #
        # maintain static position counter - this makes this routine most efficient
        # if nearby positions are asked for in succession
        #

        s.pos=0;

        #
        # get x between xlist[pos] and xlist[pos+1]
        # FIX: this is slow for random sampling.
        #
        while s.pos<len(s.xy)-2 and x>s.xy[s.pos+1][0]: s.pos += 1
        while s.pos>0           and x<s.xy[s.pos  ][0]: s.pos -= 1

        retval = s.xy[s.pos][1] + (x-s.xy[s.pos][0]) * \
                 (s.xy[s.pos+1][1]-s.xy[s.pos][1]) / \
                 (s.xy[s.pos+1][0]-s.xy[s.pos][0])

        return retval;
    pass


def interpolateCurve(experiment,weights,numPoints):
    """
    given a curve with elements (q,I) and corresponding weights, return the
    corresponding linearly interpolated curves sampled uniformly over
    numPoints.
    """
    splinedI=Interp(experiment)
    qWeights = [(experiment[i][0],weights[i]) for i in range(len(experiment))]
    splinedW=Interp(qWeights)
    qMin = min([q_I2[0] for q_I2 in experiment])
    qMax = max([q_I3[0] for q_I3 in experiment])
    qDelta = float(qMax-qMin)/(numPoints-1) if numPoints>1 else 0
    qValues = [qMin+i*qDelta for i in range(numPoints)]
    
    experiment = [(qValues[i],splinedI(qValues[i])) for i in range(numPoints)]
    weights= [splinedW(qValues[i]) for i in range(numPoints)]
    print("# datapoints sampled uniformly at %d points" % numPoints, end=' ')
    print("using linear interpolation")
    return (experiment,weights)
    pass


from simulationTools import registerTerm, registerExtraStats
from functools import reduce
registerTerm(analyze,"Solution Scattering Analysis","SolnScat",
             r"""

             For each SolnScat term (using SAXS, WAXS or SANS data), the
             following quantities are reported:

                 RMS    - root mean square deviation between calculated
                          and observed
                 Devia  - deviation between different ensemble memebrs.
                          This is only nonzero for EnsembleSimulation
                          calculations with Ne>1. It is calculated as

                           Dev^2 = 1/Ne 1/Nk \sum_i \sum_j w_j * (I_ij-I_j)^2

                           where Nk is the number of datapoints and
                           all values of I are normalized by
                           calcdScale, and multiplied by globCorrect.


                 Chi^2  - A Chi^2 value is calculated based on the
                          assumption that the weights are 1/sigma_i^2,
                          where sigma_i is the error in the observed value
                          of I_i.

             """)

registerExtraStats("SolnScatPot","Chi^2",analyze_chi2)

def useGlobs(term,globTable=[],globRules=[],
             verbose=0,
             weightFunction=lambda a: 1):
    """set up <m solnScatPot>.SolnScatPot term to use the atom globbing
    approximation.

    globTable contains a list of list of atoms with in user-defined globs.
    Atoms not specified in globTable are glob'ed by the pre-residue definitions
    in the globRules dictionary.

    globRules is a dictionary whose keys are upper case residue names
    each entry containing a list of list of atom names to be globed.

    weightFunction takes an atom as an argument and returns the relative
    weight to give it within the glob.

    Atoms in term.selection() which are not specified by globTable or by
    globRules are placed into single-atom globs.
    """

    from atomSel import AtomSel, intersection

    atomSeen = [0]*term.selection().simulation().numAtoms()
    for glob in globTable:
        for atom in glob:
            atomSeen[ atom.index() ] = 1
            pass
        pass
    
    #1.5) for every residue, apply glob rule
    #FIX: make it understand multiple segments
    from selectTools import getSegsResidues
    segsResidues = getSegsResidues(term.selection(),
                                   term.selection().simulation())
    for segid in list(segsResidues.keys()):
        for (resid,resname) in segsResidues[segid]:
            if not resname in globRules:
                continue
            resAtoms = intersection(term.selection(),
                                    'resid %d and segid "%s"' % 
                                    (resid,segid))
            nameMap={}
            for atom in resAtoms:
                # assume there is a one-to-one atom name -> atom relationship
                nameMap[atom.atomName()] = atom
                pass
            for resGlobs in globRules[resname]:
                glob=[]
                for atomName in resGlobs:
                    if atomName not in nameMap:
                        continue
                    atom = nameMap[atomName]
                        
                    if not atomSeen[ atom.index() ]:
                        glob.append( atom )
                        pass
                    atomSeen[ atom.index() ] = 1
                    pass
                if glob: globTable.append( glob )
                pass
            pass
        pass
            
    #2) make sure all atoms in term.sel are covered
    for atom in term.selection():
        if not atomSeen[ atom.index() ]:
            globTable.append((atom,))
            pass
        pass
    #3) convert list of pairs [(nElec, atom), ...] to
    #                         [(nElec,selIndex), ...]
    globDefs=[]
    singleAtomEntries=[]
    totAtoms=0
    for entry in globTable:
        glob = []
        weightSum=0
        for atom in entry:
            weightSum += weightFunction(atom)
            pass
        for atom in entry:
            glob.append( (weightFunction(atom)/weightSum,atom.index()) )
            pass
        
        if len(glob)==1: singleAtomEntries.append(entry[0])
        totAtoms += len(glob)
        globDefs.append( glob )
        pass
    from simulationWorld import world as simWorld
    if simWorld().logLevel()!='none' or verbose:
        print("number of atoms:", totAtoms)
        print("number of globs:", len(globDefs))
        print("number of single-atom globs:", len(singleAtomEntries))
        pass
    if verbose:
        print('the following atoms are in single-atom globs:')
        for a in singleAtomEntries:
            print(a.string())
            pass
        pass
    #4) enable term.useGlobs
    term.setGlobs( globDefs )
    term.setUseGlobs(1)
    term.calcGlobCorrect()
    return


def SolventAtomParams(volume=-1,radius=-1):
    from math import pi, sqrt
    if radius<0 and volume<0:
        raise Exception("volume and radius are both negative")
    #FIX: rationalize this!
    if radius<0:
        radius = volume**(1./3) / sqrt(pi)
        pass
    if volume<0:
        volume = 4*pi/3 * radius**3
        pass
    return (volume,radius)
        

solventParamSets = {}

#Svergun values - J. Appl. Cryst. 28, 768-773 (1995).
solventParamSets['svergun'] = {}
params = solventParamSets['svergun']
params['H'  ] =  SolventAtomParams( 5.15 , 1.07)
params['C'  ] =  SolventAtomParams(16.44 , 1.58)
params['CH' ] =  SolventAtomParams(21.59 , 1.73)
params['CH2'] =  SolventAtomParams(26.74 , 1.85)
params['CH3'] =  SolventAtomParams(31.89 , 1.97)
params['N'  ] =  SolventAtomParams( 2.49 , 0.84)
params['NH' ] =  SolventAtomParams( 7.64 , 1.22)
params['NH2'] =  SolventAtomParams(12.79 , 1.45)
params['NH3'] =  SolventAtomParams(17.94 , 1.62)
params['O'  ] =  SolventAtomParams( 9.13 , 1.30)
params['OH' ] =  SolventAtomParams(14.28 , 1.50)
params['S'  ] =  SolventAtomParams(19.86 , 1.68)
params['SH' ] =  SolventAtomParams(25.10 , 1.81)
params['Mg' ] =  SolventAtomParams(17.16 , 1.60)
params['P'  ] =  SolventAtomParams( 5.73 , 1.11)
params['Ca' ] =  SolventAtomParams(31.89 , 1.97)
params['Mn' ] =  SolventAtomParams( 9.20 , 1.30)
params['Fe' ] =  SolventAtomParams( 7.99 , 1.24)
params['Cu' ] =  SolventAtomParams( 8.78 , 1.28)
params['Zn' ] =  SolventAtomParams( 9.85 , 1.33)

params['H'  ] =  SolventAtomParams(radius=1.07)
params['C'  ] =  SolventAtomParams(radius=1.58)
params['CH' ] =  SolventAtomParams(radius=1.73)
params['CH2'] =  SolventAtomParams(radius=1.85)
params['CH3'] =  SolventAtomParams(radius=1.97)
params['N'  ] =  SolventAtomParams(radius=0.84)
params['NH' ] =  SolventAtomParams(radius=1.22)
params['NH2'] =  SolventAtomParams(radius=1.45)
params['NH3'] =  SolventAtomParams(radius=1.62)
params['O'  ] =  SolventAtomParams(radius=1.30)
params['OH' ] =  SolventAtomParams(radius=1.50)
params['S'  ] =  SolventAtomParams(radius=1.68)
params['SH' ] =  SolventAtomParams(radius=1.81)
params['Mg' ] =  SolventAtomParams(radius=1.60)
params['P'  ] =  SolventAtomParams(radius=1.11)
params['Ca' ] =  SolventAtomParams(radius=1.97)
params['Mn' ] =  SolventAtomParams(radius=1.30)
params['Fe' ] =  SolventAtomParams(radius=1.24)
params['Cu' ] =  SolventAtomParams(radius=1.28)
params['Zn' ] =  SolventAtomParams(radius=1.33)

params['ZN' ] =   params['Zn']
params['MN' ] =   params['Mn']

params['H2O'] =      SolventAtomParams(19.30) #from Tiede's set below


#Tiede values (private communication)
solventParamSets['tiede'] = {}
params = solventParamSets['tiede']
params['C'] =        SolventAtomParams( 9.00)
params['CH'] =       SolventAtomParams(20.00)
params['CH2'] =      SolventAtomParams(21.00)
params['CH3'] =      SolventAtomParams(33.00)
params['N'] =        SolventAtomParams(17.00)
params['NH'] =       SolventAtomParams(17.00)
params['NH2'] =      SolventAtomParams(25.00)
params['NH3'] =      SolventAtomParams(30.00) #FIX: need real NH3 value
params['O'] =        SolventAtomParams(22.00)
params['OH'] =       SolventAtomParams(25.00)
params['S'] =        SolventAtomParams(25.00)    
params['SH'] =       SolventAtomParams(34.00)
params['P'] =        SolventAtomParams(24.00)
params['Fe(2)'] =    SolventAtomParams( 8.30)
params['Cu(2)'] =    SolventAtomParams( 9.20)
params['Ca'] =       SolventAtomParams(31.00)
params['Mg(2)'] =    SolventAtomParams(21.70)
params['Mn(2)'] =    SolventAtomParams( 7.20)
params['Zn(2)'] =    SolventAtomParams(11.49)
params['Fe(3)'] =    SolventAtomParams( 8.30)
params['H2O'] =      SolventAtomParams(19.30)  	
params['Se'] =       SolventAtomParams(28.70)
params['Ni'] =       SolventAtomParams(18.14)
params['Br'] =       SolventAtomParams(26.52)
params['I'] =        SolventAtomParams(32.51)
params['Pt'] =       SolventAtomParams(21.31)
params['Re'] =       SolventAtomParams( 9.20)
params['Cl'] =       SolventAtomParams(22.45)
params['ZN' ] =   params['Zn(2)']
params['MN' ] =   params['Mn(2)']

#Xiaobing modified values (added Nov 9, 2007)
#scaling factor using 0.97 for RNA
solventParamSets['xiaobing'] = {}
params = solventParamSets['xiaobing']
params['C'] =        SolventAtomParams(16.44)
params['CH'] =       SolventAtomParams(21.59)
params['CH2'] =      SolventAtomParams(26.74)
params['CH3'] =      SolventAtomParams(31.89)
params['N'] =        SolventAtomParams(13.00)
params['NH'] =       SolventAtomParams(18.20)
params['NH2'] =      SolventAtomParams(23.40)
params['NH3'] =      SolventAtomParams(30.00) 
params['O'] =        SolventAtomParams( 9.13)
params['OH'] =       SolventAtomParams(14.28)
params['S'] =        SolventAtomParams(19.86)
params['SH'] =       SolventAtomParams(25.10)
params['P'] =        SolventAtomParams( 5.73)
params['Fe(2)'] =    SolventAtomParams( 8.30)
params['Cu(2)'] =    SolventAtomParams( 9.20)
params['Ca'] =       SolventAtomParams(31.89)
params['Mg(2)'] =    SolventAtomParams(17.16)
params['Mn(2)'] =    SolventAtomParams( 7.20)
params['Zn(2)'] =    SolventAtomParams(11.49)
params['Fe(3)'] =    SolventAtomParams( 8.30)
params['H2O'] =      SolventAtomParams(19.30)
params['Se'] =       SolventAtomParams(28.70)
params['Ni'] =       SolventAtomParams(18.14)
params['Br'] =       SolventAtomParams(26.52)
params['I'] =        SolventAtomParams(32.51)
params['Pt'] =       SolventAtomParams(21.31)
params['Re'] =       SolventAtomParams(28.70)
params['Cl'] =       SolventAtomParams(22.45)

#for backwardscompatibility
solventVolumeSets = solventParamSets

#def patterson(pot,
#                    rMax=50,
#                    rDelta=0.1):
#    """given a <m solnScatPot>.SolnScatPot, calculate the spherically-averaged
#    Patterson function from 0 to rMax in increments of rDelta.
#    This function is
#
#    <P(r)> = 1/pi \int_0^\infty dq sin(qr) q/r I(q)
#    """
#
#    from math import pi, sin, ceil
#
#
#    qValues = pot.qValues()
#    Ivalues = pot.calcd()
#
#    qDelta = qValues[1]
#
#    qMax = qValues[-1]
#
#    qCut = 2*pi/rDelta
#
#    iMax = ceil(qCut / qDelta)
#
#    qcMax = qValues
#
#    Imax = Ivalues[-1]
#
#    import sys
#    sys.stderr.write("iMax: %d ; qmax= %f\n" %(iMax,qMax))
#
#    j=0
#    Pr=[]
#    while 1:
#        r = j * rDelta
#
#        if r>rMax: break
#        
#        yValues = []
#        i=0
#        if j>0:
#            while i<iMax:
#                if i<len(qValues):
#                    q = qValues[i]
#                    I = Ivalues[i]
#                else:
#                    #FIX: for nonuniform q
#                    q = i*qDelta
#                    #asymptotic expression
#                    I = Imax * qMax**4 / q**4
#                    
##                sys.stderr.write("%f %f\n" % (q,I))
#                yValues.append( sin(q*r) * q/r * I )
#                i += 1
#                pass
#            pass
#        else:
#            for i in range(len(qValues)):
#                yValues.append( qValues[i]**2 * Ivalues[i] )
#                pass
#            pass
#
#        sum = 0
#        for i in range(len(qValues)-1):
#            dq = qValues[i+1]-qValues[i]
#            y0 = yValues[i]
#            y1 = yValues[i+1]
#            
#            sum += 0.5 * dq * (y0+y1)
#            pass
#
#        sum *= 4 * pi**3
#
#        Pr.append((r,sum))
#
#        j += 1
#        pass
#    
#
#    return Pr
#
#
#def radialPatterson(pot,
#                    rMax=50,
#                    rDelta=0.1):
#    """given a <m solnScatPot>.SolnScatPot, calculate the radial
#    Patterson function from 0 to rMax in increments of rDelta. This function is
#
#    U(r) = r^2 <P(r)> = 1/pi \int_0^\infty dq sin(qr) q I(q)
#    """
#    Pr = map(lambda (q,P): P, patterson(pot,rMax,rDelta))
#
#    j=0
#    Ur=[]
#    while 1:
#        r = j * rDelta
#
#        if r>rMax: break
#
#        P = Pr[j]
#
#        Ur.append((r,r**2 * P))
#        
#        j += 1
#        pass
#
#    return Ur

def Pr(numBins,
       rMax,
       selection="known",
       weights=None):
    """return the unnormalized pairwise distance distribution function.

    numBins specifies the granulaity of the distribution, and rMax specifies
    the upper distance cutoff.

    The calculation is performed using all atoms in selection
      [all atoms whose positions have been defined, by default].

    The optional weights array specifies a weight for each atom.
    """

    if type(selection) == type("string"):
        import atomSel
        selection = atomSel.AtomSel(selection)
        pass

    if not weights:
        from cdsVector import CDSVector_double
        weights = CDSVector_double(selection.simulation().numAtoms(),1.)
        pass

    from solnScatPot import pairDistribution
    ret = pairDistribution(selection,weights,numBins,rMax)

    # do proper averaging for ensemble simulations
    sim = selection.simulation()
    if ( sim.type() == "EnsembleSimulation"):
        from ensembleSimulation import EnsembleSimulation_currentSimulation
        esim = EnsembleSimulation_currentSimulation()

        yVals = list([x_y[1]*esim.member().weight() for x_y in ret])
        
        from ensembleSharedObj import SharedObj
        sharedObj = SharedObj()

        #collect all Pr arrays
        from cdsVector import CDSVector_double
        Pr = CDSVector_double(len(yVals),0.)
        for i in range( esim.size() ):
            if esim.member().memberIndex()==i:
                sharedObj.set(yVals)
                pass
            Pr += CDSVector_double( sharedObj.barrierGet() )
            pass

        for i in range(len(ret)):
            ret[i] = (ret[i][0], Pr[i] )
            pass
        pass

    return ret

# list of glob definitions indexed by residue name
globRules={}

# protein: from A. Grishaev, J. Wu, J. Trewhella, and A. Bax
#  JACS 127, 16621 (2005)
# FIX: should add, N-, C-termini
globRules['GLY']=[('C','O','N','HN'),
                  ('CA','HA1','HA2')]
globRules['ALA']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2','HB3')]
globRules['VAL']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB','CG1','HG11','HG12','HG13',
                   'CG2','HG21','HG22','HG23')]
globRules['LEU']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2'),
                  ('CG','HG','CD1','HD11','HD12','HD13'
                   ,'CD2','HD21','HD22','HD23')]
globRules['ILE']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB','CG2','HG21','HG22','HG23'),
                  ('CG1','HG11','HG12','CD1','HD11','HD12','HD13')]
globRules['PHE']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2'),
                  ('CG','CD1','HD1','CD2','HD2','CE1','HE1','CE2','HE2',
                   'CZ','HZ')]
globRules['TYR']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2'),
                  ('CG','CD1','HD1','CD2','HD2','CE1','HE1','CE2','HE2',
                   'CZ','OH','HH')]
globRules['TRP']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2'),
                  ('CG','CD1','HD1','CD2','NE1','HE1','CE2','CE3','HE3',
                   'CZ2','HZ2','CZ3','HZ3','CH2','HH2')]
globRules['ASP']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2'),
                  ('CG','OD1','OD2')]
globRules['GLU']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2','CG','HG1','HG2'),
                  ('CD','OE1','OE2')]
globRules['SER']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2','OG','HG')]
globRules['THR']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB','OG1','HG1','CG2','HG21','HG22','HG23')]
globRules['ASN']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2'),
                  ('CG','OD1','ND2','HD21','HD22')]
globRules['GLN']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2','CG','HG1','HG2'),
                  ('CD','OE1','NE2','HE21','HE22')]
globRules['LYS']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2','CG','HG1','HG2'),
                  ('CD','HD1','HD2','CE','HE1','HE2','NZ','HZ1','HZ2','HZ3')]
globRules['ARG']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2','CG','HG1','HG2'),
                  ('CD','HD1','HD2','NE','HE','CZ','NH1','NH2','HH11','HH12',
                   'HH21','HH22')]
globRules['HIS']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2'),
                  ('CG','ND1','HD1','CD2','HD2','CE1','HE1','NE2','HE2')]
globRules['MET']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2','CG','HG1','HG2'),
                  ('SD','CE','HE1','HE2','HE3')]
globRules['CYS']=[('C','O','N','HN'),
                  ('CA','HA','CB','HB1','HB2','SG','HG')]
globRules['PRO']=[('C','O','N'),
                  ('CA','HA','CB','HB1','HB2','CG','HG1','HG2',
                   'CD','HD1','HD2')]


# nucleic acid
globRules["CYT"] = [('P','O1P','O2P'),
                  ("O5'","C5'"),
                  ("O3'","C3'"),
                  ("O4'","C4'"),
                  ("C1'","C2'"),
                  ("N1","C2","O2"),
                  ("N3","C4","N4"),
                  ("C5","C6")   ]
globRules["GUA"] = [('P','O1P','O2P'),
                    ("O5'","C5'"),
                    ("O3'","C3'"),
                    ("O4'","C4'"),
                    ("C1'","C2'"),
                    ("N1","C2","N2"),
                    ("N3","C4"),
                    ("C5","C6","O6"),
                    ("N7","C8","N9")]
globRules["ADE"] = [('P','O1P','O2P'),
                    ("O5'","C5'"),
                    ("O3'","C3'"),
                    ("O4'","C4'"),
                    ("C1'","C2'"),
                    ("N1","C2"),
                    ("N3","C4"),
                    ("C5","C6","N6"),
                    ("N7","C8","N9")]
globRules["THY"] = [('P','O1P','O2P'),
                    ("O5'","C5'"),
                    ("O3'","C3'"),
                    ("O4'","C4'"),
                    ("C1'","C2'"),
                    ("N1","C2","O2"),
                    ("N3","C4","O4"),
                    ("C5","CM","C6")   ]

def normalize(weights,I,expt,normalizeIndex):
    """
    """
    if normalizeIndex==-3:
        # use normalization which minimizes chi^2:
        nom=0
        denom=0
#        for j in range(len(I)):
#            nom   += weights[j] * I[j]**2
#            denom += weights[j] * I[j] * expt[j]
#            pass
        # instead use numpy-type notation for speed
        from cdsVector import sum
        nom = sum(weights * I**2)
        denom = sum(weights * I * expt)
        N = nom / denom
        return 1./N * I
    elif normalizeIndex==-2:
        # normalize to average value
        from cdsVector import sum
        nom = sum(weights * I)
        denom = sum(weights)
        N = nom / denom
        return 1./N * I
    elif normalizeIndex>=0:
        return I * expt[normalizeIndex] / I[normalizeIndex]
    else:
        #no normalization
        return I
    
def chi2(weights,I,expt,normalizeIndex=-3):
    " return chi^2 using the specified normalization method"
    I = normalize(weights,I,expt,normalizeIndex)
    ret=0
    #for i in range(len(I)):
    #    ret += weights[i] * (I[i] - expt[i])**2
    #    pass
    #instead use numpy-type notation for speed
    from cdsVector import sum
    ret = sum( weights * (I-expt)**2 )
    ret /= (len(I)-1)
    return ret

def fitParams(scat,
              r0Start=None  , r0End=None  , r0Num=None,
              V0Start=None  , V0End=None  , V0Num=None,
              rhobStart=None, rhobEnd=None, rhobNum=None,
              bgStart=None  , bgEnd=None  , bgNum=None,
              writer=None,
              verbose=False):
    """
    given a SolnScatPot instance, optimize the boundary layer
    solvent density, the excluded volume effective radius, and excluded
    volume parameters. These parameters are the same as those optimized in
    the crysol program (J. Appl. Cryst. 28, 768-773 (1995)]

    scat is a SolnScatPot instance

       r0Start=0.9,  r0End=1.1,    r0Num=21,
       V0Start=0.92, V0End=1.08,   V0Num=21,
       rhobStart=0,  rhobEnd=0.02, rhobNum=20,
       bgStart=-0.5, bgEnd=0.5,    bgNum=0,

    specify the ranges for the three fitting parameters r0, V0 and rhob,
    and a constant background present in the experimental spectrum. The r0
    and V0 parameters specify a multiplicative factors relative to rm and
    Vm, respectively. By default, bg=0 and no background fitting is performed.
    bgStart and bgEnd specify values relative to the average value of the
    last 10 experimental datapoints. Default values of all parameters are
    specified above. 

    If writer is specified, it should be a function to take informational
    output.
    
    """

    if r0Start==None:   r0Start=0.9
    if r0End==None:     r0End=1.1
    if r0Num==None:     r0Num=21
    if V0Start==None:   V0Start=0.92
    if V0End==None:     V0End=1.08
    if V0Num==None:     V0Num=21
    if rhobStart==None: rhobStart=0
    if rhobEnd==None:   rhobEnd=0.02
    if rhobNum==None:   rhobNum=21
    if bgStart==None:   bgStart=-0.5
    if bgEnd==None:     bgEnd=0.5
    if bgNum==None:     bgNum=0

    if writer==None:
        import sys
        writer=sys.stdout.write
        pass

    from solnScatPot import volContrib, calcIfromF
    import cdsVector

    rho0=scat.rho0()
    rhob=scat.rhob()
    scat.setCalcBoundary(False)

    #vacuum contribution
    scat.setRho0(0)
    scat.setRhob(0)
    scat.calcEnergy()
    F = scat.getF()
    Fv = scat.getF()
    
    #excluded solvent contribution
    scat.setFormScale(0)
    scat.setRho0(rho0)
    scat.setVolumeScale(1)
    scat.setRadiusScale(1)
    scat.calcEnergy()
    Fs = scat.getF()
    
    #boundary solvent contribution
    scat.setFormScale(0)
    scat.setRho0(0)
    scat.setRhob(1)
    scat.setCalcBoundary(True)
#    scat.setSolventRadius(1.44)
#    scat.setBoundaryThickness(3)
    scat.calcBoundaryPoints()
    scat.calcEnergy()
    bVol = scat.boundaryVol()
    Fb = scat.getF()
    
    

    p1=1
    rm=scat.aveRadius()
    r0=rm
    
    p3=1
    
    #r0Start=0.9*rm; r0Delta=0.01*rm; r0Num=21
    r0Start = r0Start * rm
    r0End = r0End * rm
    r0Delta=(r0End - r0Start) / (r0Num-1) if r0Num>1 else 0
    r0Values = [r0Start+x*r0Delta for x in range(max(1,r0Num))]
    
    #p1Start=0.9; p1Delta=0.01; p1Num=21
    V0Delta = (V0End-V0Start) / (V0Num-1) if V0Num>1 else 0
    V0Values = [V0Start+x*V0Delta for x in range(max(1,V0Num))]
    
    #p3Start=0; p3Delta=0.001; p3Num=20
    rhobDelta = ( rhobEnd - rhobStart) / (rhobNum-1) if rhobNum>1 else 0
    rhobValues = [rhobStart+x*rhobDelta for x in range(max(1,rhobNum))]

    #background correction
    bgValues=[0]
    if bgNum>0:
        scat.setBG(0)
        try:
            from cdsVector import sum
            bgVal = sum(scat.expt()[-10:]) / 10
            bgStart*=bgVal; bgEnd*=bgVal
            bgDelta = ( bgEnd - bgStart) / (bgNum-1) if bgNum>1 else 0
            bgValues = [bgStart+x*bgDelta for x in range(bgNum)]
        except IndexError:
            pass
        pass
    
    
    minRMSD = None

    bg=bgMin=0
    expt=scat.expt()
    from cdsMatrix import CDSMatrix_complex
    for rhob in rhobValues:
        Fvb = Fv + rhob * Fb
        for r0 in r0Values:
            Fr0 = volContrib(r0,rm,scat.qValues(),Fs)
            for V0 in V0Values:
                #F = Fvb + V0 * Fr0
                F = CDSMatrix_complex(Fr0)
                F.scale( V0)
                F += Fvb
                I = calcIfromF(F,scat)

                for bg in bgValues:
                    lexpt= expt - bg
                    
                    # this next line sums up intensity contributions from all
                    # angles and does ensemble averaging using the
                    # weights scat.ensWeights
                    curRMSD = chi2(scat.weights(),I,lexpt,
                                   scat.normalizeIndex())
                    if  minRMSD==None or curRMSD<minRMSD:
                        minRMSD = curRMSD
                        minI = I
                        r0Min = r0
                        V0Min = V0
                        rhobMin = rhob
                        bgMin=bg
                        pass
                    if verbose:
                        print("cur param vals: ", r0,V0,rhob,bg,"cost:",curRMSD)
                        pass
                    pass
                pass
            pass
        pass
    
    
    writer("solvent fit results:\n")
    writer(" rm: %.4f\n" % rm)
    writer(" shell volume: %.2f\n" % bVol)
    from cdsVector import sum
    writer(" boundary layer contrib: %.2f\n" %
           (rhobMin*sum(calcIfromF(Fb,scat))))
    writer(" background: %f\n" % bgMin)
    writer(" min: r0/rm: %f V/V0: %f rhob: %f chi^2: %f\n" %
           (r0Min/rm,V0Min,rhobMin,minRMSD))

    scat.setBG(bgMin)
    scat.setRho0(rho0)
    scat.setFormScale(1)
    scat.setVolumeScale(V0Min)
    scat.setRhob(rhobMin)
    scat.setRadiusScale(r0Min/rm)

def fitSolventBuffer(scat,
                     qValues=None,
                     sampleData=None,
                     bufferData=None,
                     Iweights=None,
                     minQ=None,
                     maxQ=None,
                     alpha0=1,
                     sigmaAlpha=0.01,
                     sigmaRho0=0.01,
                     sigmaRhob=0.1,
                     sigmaBG=0.01,
                     rhob0=1e-3,
                     bg0=0,
                     N0=None,
                     optimizeRho0=False,
                     optimizeRhob=True,
                     writer=None,
                     computeContributions=False,
                     verbose=False):
    """
    given a SolnScatPot instance, and input sample and buffer filenames,
    optimize the buffer subtraction, background contribution and boundary
    layer contribution

    scat         -a SolnScatPot instance.
    sampleData   -a sequence containing the sample scattering curve, without
                  buffer subtraction. The length must be the same as
                  scat.expt().
    bufferData   -a sequence containing the buffer scattering curve. The
                  length must be the same as scat.expt().

    This is an adaptation of the AXES algorithm:  A. Grishaev, L.A. Guo,
    T. Irving, and A. Bax, ``Improved fitting of solution X-ray
    scattering data to macromolecular structures and structural
    ensembles by explicit water modeling'',  J. Am. Chem. Soc. 132,
    15484-15486 (2010).


    The experimental curve is computed as:

       I_expt(q) = I_samp(q) - alpha I_buff(q) + bg

    The calculated curve is given by:

       I(q) = N * < |A(q) + rhob Ab(q)|^2 >

    where N is overall normalization, A(q) is the scattering amplitude from
    solute molecule and excluded solvent,  Ab(q) is scattering from
    the solvent boundary layer, and rhob is the coefficient of the boundary
    layer scattering.

    The fit range can be reduced from the entire range by specifying
    minQ or maxQ.
    
    Initial values of the following parameters may be specified: 

        alpha0     - initial value of coefficient for buffer subtraction
        sigmaAlpha - uncertainty in alpha
        sigmaRho0  - uncertainty in rho0 - only used when fitting rho0.
        bg0        - initial value for isotropic background
        N0         - initial value of normalization for the calculated curve.
                     By default, this is taken by an overall chi^2 fit to
                     the initial experimental curve.

    The density of excluded solvent can be optimized by setting
    optimizeRho0=True.

    The computation of the optimal bound-solvent density can be disabled by
    setting optimizeRhob=False.

    The return value is a dictionary with optimal values in elements
    indexed by 'alpha', 'bg', 'N', 'rho0', and 'rhob', in addition to
    'chi2' and 'cost', the total function value where the parameters
    take their optimal values. If the computeContributions argument is
    set to True, the dictionary will additionally include 'Ivac',
    'Isol', 'Ib'and 'Inob', corresponding to intensities from the
    vacuum, excluded solvent boundary, and boundaryless contributions.

    If writer is specified, it should be a function to take informational
    output.
    
    """

    if writer==None:
        import sys
        writer=sys.stdout.write
        pass

    if not qValues:    qValues    = scat.qValuesExpt()
    if not sampleData: sampleData = scat.sampleData
    if not bufferData: bufferData = scat.bufferData
    if not Iweights:   Iweights   = scat.weights()

    if not hasattr(sampleData[0],"__len__"):
        qValues=[qValues]
        sampleData=[sampleData]
        bufferData=[bufferData]
        Iweights=  [Iweights]
        pass
    numDatasets = len(qValues)


    numSample=reduce(lambda x,y: x+y,[len(l) for l in sampleData])
    if numSample!=len(scat.expt()):
        raise Exception("sampleData has incorrect number of points: %d."
                        " Should be %d." % (numSample,len(scat.expt())))
    numBuffer=reduce(lambda x,y: x+y,[len(l) for l in bufferData])
    if numBuffer!=len(scat.expt()):
        raise Exception("bufferData has incorrect number of points: %d."
                        " Should be %d." % (numBuffer,len(scat.expt())))
    numWeights=reduce(lambda x,y: x+y,[len(l) for l in Iweights])
    if numWeights!=len(scat.expt()):
        raise Exception("Iweights has incorrect number of points: %d."
                        " Should be %d." % (numWeights,len(scat.weights())))

    startIndex=None
    if minQ!=None:
        for i,q in enumerate(scat.qValuesExpt()):
            if q>minQ:
                break
            startIndex=i
            pass
        pass
    endIndex=None
    if maxQ!=None:
        for i,q in enumerate(scat.qValuesExpt()):
            endIndex=i
            if q>maxQ:
                break
            pass
        pass
    if verbose and (minQ!=None or maxQ!=None):
        writer("fitting q=%.3f..%.3f, " %
               (minQ if minQ else 0.,
                maxQ if maxQ else scat.qValuesExpt()[-1]))
        writer("indices %d..%d\n" % (startIndex if startIndex else 0,
                                     endIndex if endIndex else endIndex))
        pass
               
    from solnScatPot import volContrib, calcIfromF, calc_dIdrhob
    import cdsVector

    #contribution without boundary layer
    rhob = scat.rhob()
    scat.setRhob(0)
    scat.calcEnergy()
    scat.setCalcBoundary(False)
    F = scat.getF()

    #vacuum contribution
    rho00=scat.rho0()
    rho0 = rho00
    scat.setRho0(0)
    scat.calcEnergy()
    Fv = scat.getF()

    #excluded solvent contribution
    scat.setFormScale(0)
    scat.setRho0(-1)
#    scat.setVolumeScale(1)
#    scat.setRadiusScale(1)
    scat.calcEnergy()
    Fs = scat.getF()
        

    #

    #boundary solvent contribution
    scat.setFormScale(0)
    scat.setRho0(0)
    scat.setRhob(1)
    scat.setCalcBoundary(True)
    scat.calcBoundaryPoints()
    scat.calcEnergy()
    bVol = scat.boundaryVol()
    writer( 'bad voxel fraction: %.3f\n'% scat.solnScat().badVoxelFraction)
    Fb = scat.getF()

    from cdsVector import CDSVector_double
#    global Isamp, Ibuf
    Isamp    = [CDSVector_double(data) for data in sampleData]
    Ibuf     = [CDSVector_double(data) for data in bufferData]
    weights0 = [CDSVector_double(data) for data in Iweights]
    weights  = [CDSVector_double(w) for w in weights0]

    def sortByQ(arrayOfArrays):
        """
        Given a sequence of sequences corresponding to datasets with
        different q-values, return a CDSVector_double with the data, sorted by
        the corresponding q-value.
        """
        singleArray=[]
        for i in range(len(qValues)):
            for j,q in enumerate(qValues[i]):
                singleArray.append( (q,arrayOfArrays[i][j]) )
                pass
            pass
        singleArray.sort( key=lambda x: x[0] ) # sort by q
        ret = CDSVector_double( [val for q,val in singleArray] )
        return ret

    exptScale=[1]*len(qValues)
    alpha0=[alpha0]*len(qValues)
    bg0=[bg0]*len(qValues)
    sigmaAlpha=[sigmaAlpha]*len(qValues)
    sigmaBG=[sigmaBG]*len(qValues)

    def computeExpts(exptScale,alpha,bg):
        """compute experimental scattering curves, taking into account
        background and buffer subtraction.
        """
        ret=[]
        for i in range(len(Isamp)):
            ret.append(  exptScale[i]*(Isamp[i] -
                                       alpha[i] * Ibuf[i] +
                                       bg[i]) )
            pass
        return ret
        

    if N0==None:
        
        expts=computeExpts(exptScale,alpha0,bg0)

        from cdsMatrix import CDSMatrix_complex
        Ffull = CDSMatrix_complex(Fb)
        Ffull.scale(rhob)
        Ffull += F

        # this next line sums up intensity contributions from all
        # angles and does ensemble averaging using the
        # weights scat.ensWeights
        I = calcIfromF(Ffull,scat)

        Inormalized = normalize(sortByQ(weights),I,sortByQ(expts),
                                normalizeIndex=-3)
        N0 = Inormalized[0] / I[0]
        pass

    
    def Jaxes(x):
        """ J = chi^2 + \sum (alpha-alpha0)^2 / sigma_alpha^2
                     [ + (rho0-rho00)^2 / sigma_rho0^2 ]
    
        where the final term is enabled only for optimizeRho0=True
        """
        numDatasets = len(qValues)
        end=3*numDatasets
        exptScale= x[:numDatasets]
        alpha    = x[numDatasets:2*numDatasets]
        bg       = x[2*numDatasets:3*numDatasets]
        M = x[end] ; end += 1
        if optimizeRhob:
            rhob = x[end] ; end += 1
            pass
        if optimizeRho0:
            rho0 = x[end] 
            pass

        expts = computeExpts(exptScale,alpha,bg)
        for i,S in enumerate(exptScale):
            weights[i] = weights0[i] * S**-2
            
    
        from cdsMatrix import CDSMatrix_complex
        Ffull = CDSMatrix_complex(Fb)
        Ffull.scale(rhob)
        if optimizeRho0:
            Ffull += Fv - rho0 * Fs
        else:
            Ffull += F
            pass
    
        # this next line sums up intensity contributions from all
        # angles and does ensemble averaging using the
        # weights scat.ensWeights
        from solnScatPot import calcIfromF
        I = calcIfromF(Ffull,scat)
        I.scale(N0*M)
    
        expt=sortByQ(expts)
        diff = (I-expt)[startIndex:endIndex]
            
        from cdsVector import sum
        curChi2 =1. / (len(diff)-1) * sum(sortByQ(weights)[startIndex:endIndex] *
                                          diff**2)
        
        ret = curChi2
        for i,alphai in enumerate(alpha):
            ret += (alphai-alpha0[i])**2/sigmaAlpha[i]**2
            pass
        for i,bgi in enumerate(bg):
            ret += bgi**2/sigmaBG[i]**2
            pass
        if optimizeRhob:
            ret += (rhob-rhob0)**2 / sigmaRhob**2
            pass
        if optimizeRho0:
            ret += (rho0-rho00)**2 / sigmaRho0**2
            pass
        return ret

    def dJaxes(x):
        " return d_J/d_alpha, d_J/d_bg, d_J/d_N,d_J/d_rhob"
    
        end=3*numDatasets
        exptScale= x[:numDatasets]
        alpha    = x[numDatasets:2*numDatasets]
        bg       = x[2*numDatasets:3*numDatasets]
        M = x[end] ; end += 1
        if optimizeRhob:
            rhob = x[end]
            end += 1
            pass
        if optimizeRho0:
            rho0 = x[end]
            pass

    
        expts = computeExpts(exptScale,alpha,bg)
    
        from cdsMatrix import CDSMatrix_complex
        Ffull = CDSMatrix_complex(Fb)
        Ffull.scale(rhob)
        if optimizeRho0:
            Ffull += Fv - rho0 * Fs 
        else:
            Ffull += F
            pass
    
        # this next line sums up intensity contributions from all
        # angles and does ensemble averaging using the
        # weights scat.ensWeights
        from solnScatPot import calcIfromF
        I = calcIfromF(Ffull,scat)
        I.scale(N0*M)

        from spline import FloatSpline
        Ispline=FloatSpline(scat.qValuesExpt(),I)
        splined=[]
        diffs=[]
        for i in range(numDatasets):
            splined.append( CDSVector_double([Ispline(q) for q in qValues[i]]) )
            diffs.append(splined[i]-expts[i])
            pass
    
        expt=sortByQ(expts)
        for i,S in enumerate(exptScale):
            weights[i] = weights0[i] * S**-2

        diff = (I-expt)[startIndex:endIndex]
        numq = len(I)

        from cdsVector import sum
        #work from here
        d_exptScale=[0] #exptScale[0]=1 always
        for i in range(1,numDatasets):
            d_exptScale.append( -2./(numq-1)/exptScale[i] *
                                sum(weights[i]*diffs[i]*expts[i]) -
                                2./exptScale[i]/(numq-1)*
                                sum(weights[i]*diffs[i]**2) )
            pass
        d_alpha=[]
        for i in range(numDatasets):
            d_alpha.append( 2./(numq-1)*exptScale[i] *
                                sum(weights[i]*diffs[i]*Ibuf[i]) +
                            2 * (alpha[i]-alpha0[i])/sigmaAlpha[i]**2 )
            pass
        d_bg=[]
        for i in range(numDatasets):
            d_bg.append(  -2./(numq-1)*exptScale[i] *   
                                sum(weights[i]*diffs[i]) +
                          2 * bg[i]/sigmaBG[i]**2         )
            pass

        lweights = sortByQ(weights)[startIndex:endIndex]
        lIbuf    = Ibuf[startIndex:endIndex]
        lI       = I[startIndex:endIndex]
    
    
        d_M = 2./M / (numq-1) * sum(lweights * (diff) * lI)
        
        ret = d_exptScale + d_alpha  + d_bg + [d_M]
        from solnScatPot import calc_dIdrhob
        if optimizeRhob:
            dIdrhob = calc_dIdrhob(Fb,F,rhob,scat)[startIndex:endIndex]
            d_rhob = 2. / (len(diff)-1)*N0*M*sum(lweights * (diff) * dIdrhob)
            d_rhob += 2 * (rhob-rhob0)/sigmaRhob**2
            ret.append(d_rhob)
            pass
    
        if optimizeRho0:
            dIdrho0 = calc_dIdrhob(Fs,Fv+rhob*Fb,-rho0,scat)[startIndex:endIndex]
            d_rho0 = -2. / (len(diff)-1)*N0*M*sum(lweights * (diff) * dIdrho0)
            d_rho0 += 2 * (rho0-rho00)/sigmaRho0**2
            ret.append(d_rho0)
            pass
    
        return ret

#    #testing
#    x =exptScale + alpha0 + bg0 + [1.0]
#    if optimizeRhob:
#        x.append(float(rhob))
#        pass
#    if optimizeRho0:
#        x.append(float(rho0))
#        pass
##    
##    x=[alpha0,bg0,1.0,rhob0]
#    writer("gradient test:\n")
#    writer("  close form: " + str(dJaxes(x)) + "\n")
#    J0=Jaxes(x)
#    eps0=1e-6
#    for i in range(len(x)):
#        xp=list(x)
#        eps = abs(eps0 * x[i]) if abs(x[i])>1e-9 else eps0
#        #writer ("eps: " + str(eps) + " " + str(x[i]) + "\n")
#        
#        xp[i] += eps
#        dJ_fd = (Jaxes(xp) - J0) / eps
#        writer ("element: " + str((i,dJ_fd, eps)) + "\n")
#        pass
#    import sys; sys.exit(0)

    x0 =exptScale + alpha0 + bg0 + [1.0]
    if optimizeRhob:
        x0.append(float(rhob0))
        pass
    if optimizeRho0:
        x0.append(float(rho00))
        pass
    from minimize import bfgs
    x=x0
    iters=0
    numMinIters=5
    for i in range(numMinIters):
        x,Jx,liters = bfgs(x,Jaxes,dJaxes,verbose=verbose,
                          stepsize0=1e-8,
                          costTol=1e-16,
                          maxIters=1000,
                          writer=writer)
        iters += liters
        pass

    if iters<2:
        writer("WARNING! small number (%d) of bfgs iterations taken in" %iters +
               " fitSolventBuffer.")
        pass




    end=3*numDatasets
    exptScale= x[:numDatasets]
    alpha    = x[numDatasets:2*numDatasets]
    bg       = x[2*numDatasets:3*numDatasets]
    M = x[end] ; end += 1
    if optimizeRhob:
        rhob = x[end]
        end += 1
        pass
    if optimizeRho0:
        rho0 = x[end]
        pass
    N=M*N0
    
    
    
    writer(" fit results:\n")
    writer(" exptScale: " + ", ".join(["%.5f" % v for v in exptScale])+"\n")
    writer(" alpha:     " + ", ".join(["%.5f" % v for v in alpha])+"\n")
    writer(" bg:        " + ", ".join(["%.5f" % v for v in bg])+"\n")
    writer(" M    : %.4f\n" % M)
    writer(" N    : %.4e\n" % N)
    writer(" rhob : %.3e\n" % rhob)
    writer(" rho0 : %.4f\n" % rho0)
    writer(" iterations taken : %d\n" % iters)

    writer(" shell volume: %.2f\n" % bVol)
    from cdsVector import sum
    writer(" boundary layer contrib: %.2f\n" %
           (rhob*sum(calcIfromF(Fb,scat))))

    expt = sortByQ( computeExpts(exptScale,alpha, bg) )
    I = sortByQ(expts)
    for i,S in enumerate(exptScale):
        weights[i] = weights0[i] * S**-2

    scat.setExpt(expt)
    scat.setWeights( sortByQ(weights) )

    if optimizeRho0:
        F = Fv - rho0 * Fs
        pass        

    Ffull = F + rhob * Fb
    I = calcIfromF(Ffull,scat)
    I.scale(N)
    diff = (I-expt)[startIndex:endIndex]

    from cdsVector import sum
    curChi2 =1. / (len(diff)-1) * sum(sortByQ(weights)[startIndex:endIndex] *
                                      diff**2)

    writer(" chi^2: %f\n" % curChi2)

    scat.setRho0(rho0)
    scat.setFormScale(1)
    scat.setRhob(rhob)
    ret={ 'S'     : exptScale,
          'alpha' : alpha,
          'bg'    : bg,
          'N'     : N,
          'rhob'  : rhob, 
          'rho0'  : rho0,
          'chi2'  : curChi2,
          'cost'  : Jx}
    if computeContributions:
        ret['Ivac'] = N * calcIfromF(Fv,scat)
        ret['Ib']   = N * rhob**2 * calcIfromF(Fb,scat)
        ret['Inob'] = N * calcIfromF(F,scat)
        ret['Isol'] = N * rho0**2 * calcIfromF(Fs,scat)
        pass
    
    return ret


