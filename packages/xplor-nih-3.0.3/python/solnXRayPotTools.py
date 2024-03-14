"""tools to aid in setup/analysis of potential terms employing solution
X-Ray scattering data.

this module provides functions to simplify the creation and
analysis of <m solnScatPot>.SolnScatPot potential terms. 

See C.D. Schwieters and G.M. Clore, ``Using Small Angle Solution
          Scattering Data in Xplor-NIH Structure Calculations ,''
          Progr. NMR Spectroscopy, accepted (2014).
"""

import solnScatPotTools

from atomSel import AtomSel
bondsByAtomIndex={}   # contains all bonds for every atom- indexed by Simulation

def create_SolnXRayPot(instanceName,
                       aSelection="not PSEUDO and not name H*",
                       experiment=0,
                       radiusScale=1,
                       volumeScale=None,
                       radiusType='volume',
                       rho0=0.334,
                       normalizeIndex=-3,
                       preweighted=True,
                       numPoints=None,
                       sampleData=None,
                       bufferData=None,
                       useInternalSpline=False,
                       minQ=None,
                       maxQ=None,
                       formFactors="original",
                       useNM=False,
                       verbose=None,
                       simulation=0):
    """
    create a <m solnScatPot>.SolnScatPot term appropriate for refining
    against solution X-Ray scattering data.

    experiment can be specified as the name of a file whose first two
    columns are the experimental q, I(q). An optional third column specifies
    a weight to use for I(q). Lines starting with # or ! are
    skipped. 

    If preweighted=False, the third entry is standard deviation instead of
    a weight, and the weight is computed as 1/sigma^2 if sigma!=0, else 0.
    If normalizeIndex>=0, the weight is also normalized [if this behavior is
    not desired, set normalizeIndex after this function has been called.]

    experiment can also be specified as a sequence of (q,I(q)) tuples.

    experimental value of q less than minQ will be ignored.
    experimental value of q greater than maxQ will be ignored.

    normalizeIndex specifies which grid point to use to normalize data.
    normalizeIndex=-2 specifies average normalization .
    normalizeIndex=-3 specifies normalization which minimizes the Chi^2 value.
    Other negative values mean no normalization.

    numPoints specifies the number of datapoints to sample. The experimental
    I(q) and errors/weights are linearly interpolated and sampled uniformly
    over the entire range. No smoothing is performed.

    If the q values are evenly spaced or numPoints is specified,
    calcType is set to 'uniform'.

    <m solnScatPot>.SolnScatPot.numAngles is initialized to 500

    <m solnScatPot>.SolnScatPot.cmpType is initialized to 'plain'.

    radiusScale is used to scale the effective atomic radii- perhaps due to
    the presence of a solvent layer.

    If specified, volumeScale is a separate scale factor for the effective
    volume. If it is not specified, this scale factor will be consistent with
    radiusScale

    rho0 is the solvent electron density. The default value is that of water.

    this function initializes effective atomic form factors including
    an excluded volume term.

    For atom j
      feff_j(q) = f_j(q) - g_j(q)

    where f_j(q) is the atomic scattering factor of group j, evaluated using a
    5-Gaussian sum.

    g_j(q) = rho0 * V_j * exp(-s**2 * PI * V_j**(2./3)) *
             volumeScale * exp(-s**2 * PI * (4PI/3)^{1/3} * (r0^2 - rm^2))

    where s=q/(2*PI), rho0 is the solvent electron density, V_j is the
    solvent-displaced volume. The second line of the expression is an overall
    expansion factor, in which 

    rm is the radius corresponding to the average atomic group volume
    (radiusType=volume) or the average atomic group radius (radiusType=radius),
    and r0 is an effective radius.

      r0 = radiusScale * rm

    and radiusScale is specified as a number approximately equal to 1.

    Specify useNM=True to indicate that q is in units of inverse nm rather than
    the default inverse angstroms.

    """

    from simulation import syncAllSimulations
    syncAllSimulations()

    if formFactors not in ("original","JBW"):
        raise Exception("formFactors %s is not valid" % formFactors)

    errors=[]
    if sampleData:
        experiment=sampleData
        preweighted=False       #FIX: should be able to override these
        useInternalSpline=True  #with argument values
        Isamp, weights = solnScatPotTools.readData(sampleData,False,
                                                   minQ,maxQ)
        from math import sqrt
        for weight in weights:
            errors.append( 1./sqrt(weight) )
        pass

    if bufferData:
        Ibuf, weights = solnScatPotTools.readData(bufferData,False,
                                                  minQ,maxQ)
        from math import sqrt
        for i,weight in enumerate(weights):
            errors[i] +=  1./sqrt(weight) 
        pass

    if (sampleData and not bufferData) or (not sampleData and bufferData):
        raise Exception("either neither or both sampleData and bufferData "
                        "arguments must be specified.")

    if not experiment:
        raise Exception("experiment or sampleData argument must be specified")

    experiment, weights = solnScatPotTools.readData(experiment,preweighted,
                                                    minQ,maxQ)

    if len(weights)==0:
        raise Exception("no datapoints read")

    if errors:
        weights=[]
        for error in errors:
            weights.append( error**-2 )
            pass
        pass

    if numPoints!=None and not useInternalSpline:
        from solnScatPotTools import interpolateCurve
        (experiment,weights) = interpolateCurve(experiment,weights,
                                                numPoints)
        pass

    qValuesExpt=[]
    IValues=[]
    for q,I in experiment:
        if useNM: q /= 10.
        qValuesExpt.append( q )
        IValues.append( I )
        pass

    qValues=[]
    if numPoints!=None and useInternalSpline:
        qmin=qValuesExpt[0]
        qmax=qValuesExpt[-1]
        delta_q = (qmax-qmin)/(numPoints-1)
        for i in range(numPoints):
            qValues.append( qValuesExpt[0] + i * delta_q )
            pass
        pass
    else:
        qValues=qValuesExpt
        pass

    from atomSel import AtomSel
    if type(aSelection)==type("string"):
        aSelection = AtomSel(aSelection)
        pass

    if not simulation:
        from simulation import currentSimulation
        simulation=currentSimulation()
        pass

    aSelection = AtomSel("(%s) and not PSEUDO" % aSelection.string(),
                         simulation)
    if not aSelection.allValid():
        raise Exception("selection contains atoms with uninitialized " +
                        "coordinates")
    

    
    #FIX: change to simply world
    from simulationWorld import SimulationWorld_world
    debug=verbose if verbose!=None else True \
           if SimulationWorld_world().logLevel()!='none' else False

    use5GaussianFormFactors()

    (fValues,
     fValues_indexed,
     atomGroups,
     groupMap       ) = genFormFactors(qValues,aSelection)


    if volumeScale==None: volumeScale = radiusScale**3

    
    if formFactors=="JBW":
        print("#using JBW form factors")
        fValues, volumeScale = convertFormFactors_jbw(aSelection, fValues,
                                                      qValues, groupMap,rho0,
                                                      radiusType, radiusScale,
                                                      solventVolume, volumeScale)     
    

        

    from solnScatPot import SolnScatPot as SolnXRayPot
    pot = SolnXRayPot(instanceName,aSelection,qValues,
                      groupMap,
                      fValues,fValues_indexed,
                      solventVolume,
                      rho0,
                      radiusScale,
                      volumeScale,
                      [],
                      qValuesExpt,
                      IValues,simulation)

    if radiusType!=pot.radiusType():
        pot.setRadiusType(radiusType)
        pass

    pot.setCmpType('plain')

    #determine whether we have uniformly spaced q values - so we can use
    #  a faster algorithm
    qtol=1e-5
    uniformQ=0
    numQ=len(qValues)
    uniformQ=1
    for qi in range(1,numQ):
        if abs((qValues[qi]-qValues[qi-1]) -
               (qValues[1]-qValues[0])     )>qtol:
            uniformQ=0
            break
        pass

    if uniformQ:
        pot.setCalcType('uniform')
    else:
        print('WARNING: q values do not take uniform values.', end=' ')
        print('This will be slow.')
        pass

    pot.setWeights(weights)

    pot.setNormalizeIndex( normalizeIndex )
    pot.setNumAngles( 500 )

    if debug:
        print("create_SolnXRayPot: read I(q) at %d values of q" % \
              len(experiment))
        if uniformQ: print("\tusing uniform q algorithm")
        print("\taverage group radius=%f" % pot.aveRadius())

        keys=list(atomGroups.keys())
        keys.sort()

        totVolume = 0
        totElectrons = 0
        totEffScat0  =0
        totCount = 0
        print("%-8s %7s %9s %9s %9s" % ('group','count',
                                        'electrons','eff. f(0)','volume'))
        for key in keys:
            electrons = calc_f(key,0.)
            effScat0 = fValues[0][key]
            totCount += atomGroups[key]
            totElectrons += atomGroups[key] * electrons
            totEffScat0 += atomGroups[key] * effScat0
            totVolume    += atomGroups[key] * solventVolume[key][0]
            print("%8s %7d %9.2f %9.2f %9.2f" % (key,atomGroups[key],
                                                 electrons,fValues[0][key],
                                                 solventVolume[key][0]))
            pass
        print("%-8s %7d %9.2f %9.2f %9.2f" % ('total',totCount,
                                              totElectrons,totEffScat0,
                                              totVolume))
        
        pass


    if sampleData:
        from cdsVector import CDSVector_double
        pot.sampleData = CDSVector_double([I for q,I in Isamp])
        pot.bufferData = CDSVector_double([I for q,I in Ibuf])


    if SimulationWorld_world().logLevel()=='debug':
        pot.setVerbose(1)
        pass

    return pot

create_solnXRayPot = create_SolnXRayPot # for backwards compatibility

def genFormFactors(qValues,aSelection):

    from math import pi, sqrt

    atomGroups={}# "H2O" : 0 } #initialize with water
    sim=aSelection.simulation()
    groupMap=[""]*sim.numAtoms()
    bbAI=[[] for i in range(sim.numAtoms())]
    for i in range(sim.numBonds()):
        bPair = sim.bondPairByID(i)
        bbAI[bPair[0]].append( bPair[1])
        bbAI[bPair[1]].append( bPair[0])
        pass
    global bondsByAtomIndex
    bondsByAtomIndex[sim.name()] = bbAI
    
    for atom in aSelection:
        name = getCombinedAtoms2(atom,bbAI,atomGroups)
        groupMap[atom.index()] = name
        pass

    fValues=[]
    fValues_indexed=[]
    for q in qValues:
        fValues.append( calcFormFactors(q,atomGroups) )
        fValues_indexed.append( calcFormFactors_indexed(q) )
        pass

    return fValues, fValues_indexed, atomGroups, groupMap



def calcFormFactors(q,atomGroups):
    """
    calculate atomic form factors for combined atoms in atomGroups.
    """
    ret={}
    for group in list(atomGroups.keys()):
        ret[group] = calc_f(group,q)
        pass
    return ret

def calcFormFactors_indexed(q):
    """
    calculate atomic form factors for special combined atoms identified by
    an atom selection and references by atom index.
    """
    ret={}
    # FIX: selection string: should reference pot's
    # DNA's O1P/O2P are special: they share 1 H and have -1 charge
    for atom in AtomSel("name O1P and chemical XO2"):
        ret[atom.index()] = calc_f('O-',q)
        pass
        
    return ret


####### contributed by Jinbu Wang prior to 2015/11/12
####### added parameters and functions for tunning SAXS paramenters
 
raw_corr_facts = [
  (  0.000, 0.9882 ),(  0.020, 0.9908 ),(  0.040, 1.0087 ),(  0.060, 1.0226 ),(  0.080, 1.0242 ),
  (  0.100, 1.0161 ),(  0.120, 1.0087 ),(  0.140, 1.0092 ),(  0.160, 1.0074 ),(  0.180, 1.0115 ),
  (  0.200, 1.0256 ),(  0.220, 1.0600 ),(  0.240, 1.1042 ),(  0.260, 1.1333 ),(  0.280, 1.1394 ),
  (  0.300, 1.1171 ),(  0.320, 1.0862 ),(  0.340, 1.0438 ),(  0.360, 1.0128 ),(  0.380, 0.9815 ),
  (  0.400, 0.9566 ),(  0.420, 0.9156 ),(  0.440, 0.8926 ),(  0.460, 0.8681 ),(  0.480, 0.8552 ),
  (  0.500, 0.8409 ),(  0.520, 0.8252 ),(  0.540, 0.8135 ),(  0.560, 0.8065 ),(  0.580, 0.8209 ),
  (  0.600, 0.8630 ),(  0.620, 0.9233 ),(  0.640, 0.9888 ),(  0.660, 1.0182 ),(  0.680, 1.0088 ),
  (  0.700, 0.9889 ),(  0.720, 0.9643 ),(  0.740, 0.9655 ),(  0.760, 0.9793 ),(  0.780, 1.0092 ),
  (  0.800, 1.0410 ),(  0.820, 1.0848 ),(  0.840, 1.1034 ),(  0.860, 1.1120 ),(  0.880, 1.1280 ),
  (  0.900, 1.1297 ),(  0.920, 1.1406 ),(  0.940, 1.1524 ),(  0.960, 1.1571 ),(  0.980, 1.1548 ),
  (  1.000, 1.1426 ),(  1.020, 1.1243 ),(  1.040, 1.0919 ),(  1.060, 1.0877 ),(  1.080, 1.0906 ),
  (  1.100, 1.0962 ),(  1.120, 1.1032 ),(  1.140, 1.1099 ),(  1.160, 1.1168 ),(  1.180, 1.1184 ),
  (  1.200, 1.1359 ),(  1.220, 1.1322 ),(  1.240, 1.1400 ),(  1.260, 1.1410 ),(  1.280, 1.1480 ),
  (  1.300, 1.1637 ),(  1.320, 1.1837 ),(  1.340, 1.2090 ),(  1.360, 1.2233 ),(  1.380, 1.2403 ),
  (  1.400, 1.2533 ),(  1.420, 1.2506 ),(  1.440, 1.2448 ),(  1.460, 1.2534 ),(  1.480, 1.2665 ),
  (  1.500, 1.2850 ),(  1.520, 1.3130 ),(  1.540, 1.3078 ),(  1.560, 1.2999 ),(  1.580, 1.2892 ),
  (  1.600, 1.2841 ),(  1.620, 1.2872 ),(  1.640, 1.2992 ),(  1.660, 1.3217 ),(  1.680, 1.3346 ),
  (  1.700, 1.3528 ),(  1.720, 1.3497 ),(  1.740, 1.3392 ),(  1.760, 1.3151 ),(  1.780, 1.2857 ),
  (  1.800, 1.2662 ),(  1.820, 1.2672 ),(  1.840, 1.2645 ),(  1.860, 1.2495 ),(  1.880, 1.2466 ),
  (  1.900, 1.2305 ),(  1.920, 1.2058 ),(  1.940, 1.1776 ),(  1.960, 1.1412 ),(  1.980, 1.1172 ),
  (  2.000, 1.0763 ),(  2.020, 1.0666 ),(  2.040, 1.0768 ),(  2.060, 1.0707 ),(  2.080, 1.0914 ),
  (  2.100, 1.1116 ),(  2.120, 1.1386 ),(  2.140, 1.1491 ),(  2.160, 1.1564 ),(  2.180, 1.1809 ),
  (  2.200, 1.1803 ),(  2.220, 1.1836 ),(  2.240, 1.1935 ),(  2.260, 1.1840 ),(  2.280, 1.1921 ),
  (  2.300, 1.2033 ),(  2.320, 1.2357 ),(  2.340, 1.2414 ),(  2.360, 1.2490 ),(  2.380, 1.2299 ),
  (  2.400, 1.2232 ),(  2.420, 1.2277 ),(  2.440, 1.2195 ),(  2.460, 1.1966 ),(  2.480, 1.2053 ),
  (  2.500, 1.1878 ),(  2.520, 1.1751 ),(  2.540, 1.1849 ),(  2.560, 1.1929 ),(  2.580, 1.1883 ),
  (  2.600, 1.2266 ),(  2.620, 1.2503 ),(  2.640, 1.2219 ),(  2.660, 1.2339 ),(  2.680, 1.2391 ),
  (  2.700, 1.2175 )]

corr_facts = [
  (  0.000, 1.0000 ), 
  (  0.200, 1.0000 ),(  0.220, 1.0600 ),(  0.240, 1.1042 ),(  0.260, 1.1333 ),
  (  0.380, 0.9815 ),
  (  0.680, 1.0088 ),(  0.780, 1.0092 ),
  (  0.920, 1.1406 ),(  0.940, 1.1524 ),(  0.960, 1.1571 ),(  0.980, 1.1548 ),
  (  1.300, 1.1637 ),(  1.320, 1.1837 ),(  1.340, 1.2090 ),(  1.360, 1.2233 ),
  (  3.000, 1.2375 )]

def getAveRadius_jbw(sel, groupMap,radius_type,radius_scale,solvent_volume):
    Vm = 0.
    rm = 0.
    #print solvent_volume
    for atom in sel:
        groupName = groupMap[ atom.index() ]
        Vm += solvent_volume[groupName][0]
        rm += solvent_volume[groupName][1]
     ##// calc Vm (ave group vol) ; rm is corresponding radius
    Vm /= sel.size()
    rm /= sel.size()
    print("# num_atom= %4d rm(rm)= %8.4f" % (sel.size(), rm ))
    if ( radius_type=="volume"):
        rm = pow(3*Vm/(4*pi),1./3)
    else:
        rm = rm/sel.size()
        print("# rm(rm)= %8.4f" % rm) 
    #return rm 
    r0 = rm * radius_scale    
    print("#*** r0=%8.4f  rm= %8.4f  radius_scale=%8.4f Vm=%8.4f" % (r0, rm, radius_scale, Vm))
    return rm, r0

def calc_g_jbw(rho0, volume_scale, V,q,r0,rm):
    q /= 2*pi
    g = rho0 * V * exp(-pi * q**2 * pow(V,2./3)) *		\
		volume_scale *  exp(-pi * q**2 * pow(4*pi/3,2./3) *
                                    (r0**2 - rm**2))
    return g

def convertFormFactors_jbw(sel, fValues, qValues, groupMap,rho0,
                           radiusType, radius_scale, solvent_volume,volum_scale):
    print("# ---------------------")
    #print "#***  num_f: ", len(fValues)
    #print "# fvalues: ", fValues
    rm,r0 = getAveRadius_jbw(sel, groupMap, radiusType,
                             radius_scale,solvent_volume) 
    from solnScatPotTools import Interp
    corr_interp = Interp(corr_facts)
    
    num_q = len(fValues)
    print("# rh0= %8.4f" % rho0)
    for i in range(num_q):
        q = qValues[i]
        fv = fValues[i]
        q_corr = corr_interp(q)
        #print "# q_corr= %8.4f" % q_corr
        #print "# q=%8.4f  " % q, fv 
        for group_name in list(fv.keys()):
            fvg = fv[group_name]
            gg = calc_g_jbw(rho0,volum_scale, solvent_volume[group_name][0],
                            q, r0, rm)
            fv[group_name] -= gg 
            fv[group_name] = fv[group_name] *  q_corr
            if i == 0:
                print("# %4s: %8.4f  %8.4f  %8.4f v=%8.4f r=%8.4f" % (group_name, fvg, gg, fv[group_name], solvent_volume[group_name][0], solvent_volume[group_name][1]))
        #break 
    print("# ---------------------")
    return fValues, 0.0

def updateQValues(pot,newQ):
    """
    update the <m solnScatPot>.SolnScatPot term with new values of scattering
    vector amplitude newQ.

    """

    if len(newQ) != len(pot.qValues()):
        raise exception("newQ is the wrong size")
    

    aSelection = pot.selection()
    simulation = aSelection.simulation()
    
    atomGroups={}
    groupMap=[""]*simulation.numAtoms()
    bbAI = bondsByAtomIndex[simulation.name()]
    for atom in aSelection:
        name = getCombinedAtoms2(atom,bbAI,atomGroups)
        groupMap[atom.index()] = name
        pass


    numQ = len(newQ)

    fValues=[]
    fValues_indexed=[]
    for qi in range(numQ):
        q = newQ[qi]
        fValues.append( calcFormFactors(q,atomGroups) )
        fValues_indexed.append( calcFormFactors_indexed(q) )
        pass

    pot.setQValues(newQ)
    pot.setFormFactors( fValues )
    pot.setIFormFactors( fValues_indexed )

    pot.convertFormFactors(fValues ,fValues_indexed,  groupMap)

    return

from math import exp, pi

def use4GaussianFormFactors():

    global calc_f
    def calc_f(name,q):
    
        (a,b,c) = atomicFormFactor[ name[0] ]

        #    print 'calc_f: a:', name, a
        #    print 'calc_f: b:', name, b

        f = c
        q /= 2*pi
        for i in range(4):
            f += a[i] * exp(-b[i] * q**2)
            pass
        return f

    global atomicFormFactor
    atomicFormFactor = {}

    #values used with XPLOR - from International Tables
    # for 4-gaussian sum + constant
    #   a[n]                               b[n]                             c
    atomicFormFactor['C'] = \
     ((2.31000,1.02000,1.58860,.865000),(20.8439,10.2075,.568700,51.6512),.215600)
    atomicFormFactor['N'] = \
     ((12.2126,3.13220,2.01250,1.16630),(.005700,9.89330,28.9975,.582600),-11.529)
    atomicFormFactor['O'] = \
     ((3.04850,2.28680,1.54630,.867000),(13.2771,5.70110,.323900,32.9089),.250800)
    atomicFormFactor['S'] = \
     ((6.90530,5.20340,1.43790,1.58630),(1.46790,22.2151,.253600,56.1720),.866900)
    atomicFormFactor['P'] = \
     ((6.43450,4.17910,1.78000,1.49080),(1.90670,27.1570,0.52600,68.1645),1.11490)

    return

def use5GaussianFormFactors():

    global calc_f
    def calc_f(name,q):
    
        (a,b) = atomicFormFactor[ name ]

        #    print 'calc_f: a:', name, a
        #    print 'calc_f: b:', name, b

        f = 0
        q /= 4*pi
        for i in range(5):
            f += a[i] * exp(-b[i] * q**2)
            pass
        return f

    global atomicFormFactor
    atomicFormFactor = {}
    z = atomicFormFactor #alias
    # values from D. Tiede (private communication)
    #  uses a 5-Gaussian sum
    #
    z[  "C"] = ((  1.81015,   2.03806,   1.25283,   0.57653,   0.32247,) , (  30.18233, 13.07207,  0.69408,   0.16703, 69.95614 ))
    z[ "CH"] = ((  1.28885,   3.24706,   1.31939,   0.49963,   0.64505,) , (  10.65674, 29.77388,  0.66138,   0.14458, 66.27304 ))
    z["CH2"] = ((  0.91282,   4.22234,   0.36957,   1.431159,  1.06372,) , (   8.58306, 30.48037,  0.09802,   0.61042, 66.31128 ))
    z["CH3"] = ((  1.57514,   3.46444,   0.1925,    2.74389,   1.02131,) , (   0.53977, 49.90266,  0.0002,   33.43777,  7.80685 ))
    z[  "N"] = ((  2.2414,    1.33805,   0.47731,   2.40679,   0.53634,) , (   8.54383,  0.47709,  0.09489,  19.70556, 46.95773 ))
    z[ "NH"] = ((  0.32295,   2.13372,   1.48215,   0.92588,   3.13508,) , (   0.04482,  8.03069,  0.43969,  51.76493, 22.9704  ))
    z["NH2"] = ((  0.2621,    1.9553,    1.53407,   1.26507,   3.98324,) , (   0.02,     7.5474,   0.4225,   56.70841, 25.02174 ))
#FIX: need real NH3 value - this obtained by fitting 5 Gaussians to
#     glob expression for q=0..3 A^{-1}
    z["NH3"] = (
        [1.9600690174078528, 1.2613566493167097, 1.2200364734935332, 1.0088784078784709, 4.5496594519],
        [11.887194879202859, 6.9505330081410479, 4.1413965977016662, 3.4675129926353132, 43.924152848529658])

    z[  "O"] = ((  0.68831,   1.10353,   0.70943,   2.58052,   2.91796,) , (  35.4895,   0.40474,  0.11101,   6.02703, 14.45603 ))
    z[ "OH"] = ((  1.1298,    1.02317,   0.79278,   2.9258,    3.12788,) , (  44.60773,  0.42288,  0.1228,    6.24215, 18.21448 ))
    z[  "S"] = ((  6.94238,   0.83089,   1.42609,   5.12232,   1.67781,) , (   1.45773,  0.0001,   0.23462,  21.94413, 54.89569 ))
    z[ "SH"] = ((  2.55447,   1.82612,   4.87293,   4.23984,   3.50432,) , (   1.000,    0.0744,   1.64766,  21.56838, 50.4351  ))
    z[  "P"] = ((  5.95119,   1.81595,   1.56221,   4.15127,   1.41842,) , (   1.95905,  0.81481,  0.05359,  27.05272, 67.60129 ))
    z["Fe2"] = (( 11.0424,    7.37400,   4.13460,   0.43990,   1.00970,) , (   4.65380,  0.30530, 12.0546,   31.2809 ,  0.00001 ))
    z["Cu2"] = (( 11.8168,    7.11181,   5.78135,   1.14523,   1.14431,) , (   3.37484,  0.244078, 7.98760,  19.8970 ,  0.00001 ))
    z["Ca2"] = (( 15.6348,    7.95180,   8.43720,   0.85370, -14.875,  ) , (  -0.00740,  0.608900,10.3116,   25.9905 ,  0.00001 ))
    z["Mg2"] = ((  3.49880,   3.83780,   1.32840,   0.84970,   0.48530,) , (   2.16760,  4.75420,  0.18500,  10.1411 ,  0.00001 ))
    z["Mn2"] = (( 10.8061,    7.36200,   3.52680,   0.21840,   1.08740,) , (   5.27960,  0.34350, 14.3430,   41.3235 ,  0.00001 ))
    z["Zn2"] = (( 11.9719,    7.38620,   6.46680,   1.39400,   0.78070,) , (   2.99460,  0.20310,  7.08260,  18.0995 ,  0.00001 ))
    z["Fe3"] = (( 11.1764,    7.38630,   3.39480,   0.07240,   0.97070,) , (   4.61470,  0.30050, 11.6729,   38.5566 ,  0.00001 ))
    z["H2O"] = ((  3.68835,   1.32318,   0.48264,   1.5622,    2.94289,) , (  20.80753,  0.36266,  0.07079,  48.264  ,  6.15175 ))
    z[ "Se"] = (( 17.0006,    5.81960,   3.97310,   4.35430,   2.84090,) , (   2.40980,  0.27260, 15.2372,   43.8163 ,  0.00001 ))
    z[ "Ni"] = (( 12.8376,    7.29200,   4.44380,   2.38000,   1.03410,) , (   3.87850,  0.25650, 12.1763,   66.3421 ,  0.00001 ))
    z[ "Br"] = (( 17.1789,    5.2358,    5.6377,    3.98510,   2.95570,) , (   2.17230, 16.5796,   0.2609,   41.4328 ,  0.00001 ))
    z[  "I"] = (( 20.1472,   18.9949,    7.5138,    2.27350,   4.07120,) , (   4.34700,  0.38140, 27.7660,   66.8776 ,  0.00001 ))
    z[ "Pt"] = (( 13.7128,    4.04642,  14.6598,   30.12358,  15.43967,) , (   0.0001,  48.54183,  2.90269,   0.81053, 11.38207 ))
    z[ "Re"] = (( 14.79924,  31.38331,   5.78461,   2.76129,  20.26665,) , (   7.10155,  1.40708, 23.36672,  79.77528,  0.09357 ))
    z[ "Cl"] = ((  2.55403,   1.95231,   4.5842,    6.13071,   1.77844,) , (   1.1686,   0.0765,   1.16857,  18.27448, 46.11295 ))
    z[  "H"] = ((  0.2117,    0.49227,   0.02521,   0.00722,   0.26364,) , (   6.04991, 17.64621,  1.36048, 119.67414, 42.7858  ))
    z["Cl1"] = (( 18.2915,    7.2084,    6.5337,    2.33860, -16.378,  ) , (   0.00660,  1.17170, 19.5424,   60.4486 ,  0.00001 ))
    z["Br1"] = (( 17.1718,    6.3338,    5.5754,    3.72720,   3.1776, ) , (   2.20590, 19.3345,   0.2871,   58.1535 ,  0.00001 ))
    z[ "I1"] = (( 20.2332,   18.9970,    7.8069,    2.88680,   4.07140,) , (   4.35790,  0.38150, 29.5259,   84.9304 ,  0.00001 ))
    z[ "O-"] = ((  4.1916,    1.63969,   1.52673, -20.307,    21.9412, ) , (  12.8573,   4.17236, 47.0179,   -0.01404,  0.00001 ))

    z["ZN"] = z["Zn2"]
    z["MN"] = z["Mn2"]
    z["MG"] = z["Mg2"]
    return

def sinc(x):
    from math import sin
    if x<1e-7:
        return 1.-1./6*x**2
    else:
        return sin(x) / x
    pass


def f_glob(q,atoms,dist):
    """ determine globbed value of f (using calc_f) given a list of atom names
    and a dictionary containing atom-atom distances. Can only be used with very
    simple glob defs.
    """
    ret=0.
    from math import sqrt
    for i in range(len(atoms)):
        ret += calc_f(atoms[i],q)**2
        for j in range(i+1,len(atoms)):
            d = dist[atoms[i] + '-' + atoms[j]]
            ret +=  2 * calc_f(atoms[i],q) * calc_f(atoms[j],q) * sinc(q*d)
            pass
        pass
    return sqrt(ret)

metalAtoms=['ZN','MN','MG']

def getCombinedAtoms2(atom,bondsByAtomIndex,atomGroups=0):
    """ given an atom object, determine the combination of heavy and light
    atoms in this group.

    atomGroups is used to hold per-group counts.
    """
    ret=''
    name = atom.atomName()
    if name[:2] in metalAtoms:
        ret= name[:2]
    else:
        ret=name[0]
        pass
    sim=atom.simulation()
    cnt=0
    for index in bondsByAtomIndex[atom.index()]:
        if sim.atomByID(index).atomName()[0]=="H": cnt += 1
        pass
    if cnt>0:
        ret += "H"
        if cnt>1:
            ret += "%d" % cnt
            pass
        pass
    if atomGroups==0:
        return ret

    try:
        atomGroups[ret] += 1
    except KeyError:
        atomGroups[ret] = 1
        pass
    return ret

import solnScatPotTools
#FIX: solventVolume is more appropriately named solventParams
solventVolume = solnScatPotTools.solventVolumeSets['svergun']

from solnScatPotTools import globRules

# for backwards compatibility
from solnScatPotTools import solventVolumeSets


def numElectrons(atom):
    """
    return the number of electrons of a given atom.
    """
    name = getCombinedAtoms2(atom,bondsByAtomIndex[atom.simulation().name()])
    return calc_f(name,0.)


def useGlobs(term,globTable=[],globRules=globRules,verbose=0):
    """set up <m solnScatPot>.SolnScatPot term to use the atom globbing
    approximation for solution X-ray scattering.

    globTable contains a list of list of atoms with in user-defined globs.
    Atoms not specified in globTable are glob'ed by the pre-residue definitions
    in the globRules dictionary.

    globRules is a dictionary whose keys are upper case residue names
    each entry containing a list of list of atoms to be globed.

    Atoms in term.selection() which are not specified by globTable or by
    globRules are placed into single-atom globs.
    
    """
    import solnScatPotTools

    solnScatPotTools.useGlobs(term,globTable,globRules,verbose,numElectrons)

    return


    
