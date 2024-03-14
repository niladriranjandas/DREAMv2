
import solnScatPotTools



def create_SANSPot(instanceName,
                   aSelection="all",
                   experiment=0,
                   normalizeIndex=-3,
                   preweighted=True,
                   numPoints=None,
                   minQ=None,
                   maxQ=None,
                   weightByExpt=False,
                   globProtons=0,
                   radiusScale=1,
                   volumeScale=None,
                   cmpType='log',
                   fractionD2O=1,
                   exchangeFraction=0.9,
                   fractionDeuterated=0,
                   verbose=False,
                   simulation=0,
                   altDeuteratedSels=[]):
    """
    create a <m solnScatPot>.SolnScatPot term appropriate for refining
    against solution neutron scattering data.

    For efficient refinement against scattering curves, the experimental SANS
    data should be extrapolated to q=0 and downsampled at constant spacing in
    q. 

    aSelection specifies the atoms which the SANS experiment measures.

    normalizeIndex specifies which grid point to use to normalize data.
    A value of -2 invokes normalization based on the average (calculated 
    and experimental) values. A value of -3 causes normalization to be chosen
    such that the chi^2 is minimized. Other values mean no normalization. 

    experiment is a string or filename which contains the lines of the form
    q I(q) w(q)
    the # symbol introduces a comment. If preweighted=False, the third
    entry is standard deviation instead of a weight, and the weight is computed
    as 1/sigma^2 if sigma!=0, else 0. If normalizeIndex>=0, sigma is also
    normalized.

    numPoints specifies the number of datapoints to sample. The experimental
    I(q) and errors/weights are linearly interpolated and sampled uniformly
    over the entire range. No smoothing is performed.

    experimental value of q less than minQ will be ignored.
    experimental value of q greater than maxQ will be ignored.

    If the q values are evenly spaced or numPoints is specified,
    calcType is set to 'uniform'.

    If weightByExpt is true, weights are multiplied by (I(0)/I(q))^2, where
    I(q) is the scattering intensity at scattering amplitude q.

    globProtons specifies whether or not to glob protons on to their bonded
    heavy atoms

    radiusScale specifies a radius correction factor

    If specified, volumeScale is a separate scale factor for the effective
    volume. If it is not specified, this scale factor will be consistent with
    radiusScale.

    cmpType should be 'plain' or 'log' and corresponds to cmpType specified
    for <m solnScatPot>.SolnScatPot.

    fractionD2O specifies the fraction of D2O in the solvent.

    exchangeFraction specifies the percentage of exchangable protons which
    actually exchange.

    fractionDeuterated specifies the deuteration fraction of non-exchangable
    protons. The altDeuteratedSels argument allows one to specify different
    deuteration fractions for specified <m atomSel>.AtomSels. This argument
    should be a sequence of tuples (atom selection, value), where the value
    is the deuteration fraction of non-exchangable protons.
    """

    debug=verbose
    #FIX: change to simply world
    from simulationWorld import SimulationWorld_world
    if SimulationWorld_world().logLevel()!='none': debug=True


    experiment, weights = solnScatPotTools.readData(experiment,preweighted,
                                                    minQ,maxQ)
    
    if weightByExpt:
        for cnt in range(len(experiment)):
            weights[cnt] /= experiment[cnt][1]**2
            pass
        pass
    
    if len(weights)==0:
        raise Exception("no datapoints read")

    if numPoints!=None:
        from solnScatPotTools import interpolateCurve
        (experiment,weights) = interpolateCurve(experiment,weights,numPoints)
        pass
    
    if weightByExpt:
        index=0
        if normalizeIndex>=0: index = normalizeIndex
        I0 = experiment[index][1]
        weights = [w*I0**2 for w in weights]
        pass

# now handled by the potential
#    if not preweighted and normalizeIndex>=0:
#        I0 = experiment[normalizeIndex][1]
#        weights = map(lambda w: w*I0**2, weights)
#        pass

    from atomSel import AtomSel
    if type(aSelection)==type("string"):
        aSelection = AtomSel(aSelection)
        pass

    aSelection = AtomSel("(%s) and not PSEUDO" % aSelection.string())
    
    if not simulation:
        from simulation import currentSimulation
        simulation=currentSimulation()
        pass

    #FIX: solvent (d2o) electron density
    global rho0
    #rho0 = 0.00562# water 10**(-12) cm/A^3
    #rho0 = 0.06404 # d2o

    rho0 = -(1. - fractionD2O) * 0.00562 + fractionD2O * 0.06404

    numQ = len(experiment)
    qValues=[]
    IValues=[]
    for qi in range(numQ):
        (q,I) = experiment[qi]
        qValues.append( q )
        IValues.append( I )
        pass

    (fValues,
     fValues_indexed,
     atomGroups,
     groupMap,
     excludedVol)  = genFormFactors(qValues,aSelection,exchangeFraction,
                                    fractionDeuterated,fractionD2O,
                                    altDeuteratedSels)
                                  
    if volumeScale==None:
        volumeScale = radiusScale**3

    from solnScatPot import SolnScatPot
    pot = SolnScatPot(instanceName,aSelection,qValues,
                      groupMap,
                      fValues,fValues_indexed,
                      solventVolume,
                      rho0,
                      radiusScale,
                      volumeScale,
                      [],qValues,IValues,simulation)

    pot.setCmpType(cmpType)

    #determine whether we have uniformly spaced q values - so we can use
    #  a faster algorithm
    qtol=1e-5
    uniformQ=0
    if len(qValues)>1:
        uniformQ=1
        for qi in range(1,numQ):
            if abs((qValues[qi]-qValues[qi-1]) -
                   (qValues[1]-qValues[0])    )>qtol:
                uniformQ=0
                break
            pass
        pass

    if uniformQ: pot.setCalcType('uniform')

    pot.setWeights(weights)

    if pot.cmpType()=='log' and normalizeIndex<-2:
        raise Exception('log cmpType and normalizeIndex<-2 are not compatible')

    pot.setNormalizeIndex( normalizeIndex )
    pot.setNumAngles( 500 )

    if debug:
        print("create_SANSPot: read I(q) at %d values of q" % \
              len(experiment))
        if uniformQ: print("\tusing uniform q algorithm")

        keys=list(atomGroups.keys())
        keys.sort()

        excludedVol=0
        for atom in aSelection: 
            name = atom.atomName()
            if name in metalAtoms:
                name=name
            else:
                name=name[0]
                pass
            xName = exchangableName(atom)
            if xName: name = xName
            V = solventVolume[name][0]
            excludedVol += V
            pass

        print('\tuncorrected excluded Vol: %.2f Angstrom^3'%excludedVol)
        print('\taverage atomic radius:    %.3f Angstrom  '%pot.aveRadius())

        totCount = 0
        print("%-8s %7s" % ('group','count'))
        for key in keys:
            totCount += atomGroups[key]
            print("%8s %7d" % (key,atomGroups[key]))
            pass
        print("%-8s %7d" % ('total',totCount))
        
        pass

    
    if globProtons:
        globTable = groupProtons(aSelection)
        useGlobs(pot,globTable)
        pass

    if SimulationWorld_world().logLevel()=='debug':
        pot.setVerbose(1)
        pass

    return pot


def genFormFactors(qValues,aSelection,
                   exchangeFraction,fractionDeuterated,fractionD2O,
                   altDeuteratedSels):
    simulation=aSelection.simulation()
    # this use of exchangeFraction is not consistent with the documentation
    #for Crysol
    fracD = exchangeFraction * fractionD2O
    scattLength['exchH1'] =  fracD*scattLength['D']+(1 - fracD)*scattLength['H']

    # backbone HN atoms: frac is reduced by factor 0.9
    fracD = 0.9 * exchangeFraction * fractionD2O
    scattLength['exchH2'] =  fracD*scattLength['D']+(1 - fracD)*scattLength['H']

#    print 'd-lengths', scattLength['exchH1'],    scattLength['exchH2']


    solventVolume['exchH1'] = solventVolume['H'] 
    solventVolume['exchH2'] = solventVolume['H']

    from math import pi, sqrt
    groupMap=[""]*simulation.numAtoms()
    atomGroups={}
    excludedVol=0
    for atom in aSelection:
        #NOTE: scattering length lookup is currently based on the
        #first character of the atomName. This is fragile, and should be fixed.
        name = atom.atomName()
        if name in metalAtoms:
            name=name
        else:
            name=name[0]
            pass
        xName = exchangableName(atom)
        if xName: name = xName
        groupMap[atom.index()] = name
        if name in atomGroups:
            atomGroups[name] += 1
        else:
            atomGroups[name] = 1
        pass

    lScattLength = dict(scattLength)
    f=fractionDeuterated
    lScattLength['H'] = (1-f) * lScattLength['H'] + f * scattLength['D']

    #deal with alternatively deuterated selections
    cnt=0
    from selectTools import convertToAtomSel
    from atomSel import intersection
    for (sel,val) in altDeuteratedSels:
        hSel=convertToAtomSel(sel)
        if val<0 or val>1:
            raise Exception("bad value for altDeuteratedSel with associated" +
                            " selection: "+hSel.string())
        hSel = intersection("name H*",hSel)
        if len(hSel)==0:
            raise Exception("no protons in  altDeuteratedSel selection: " +
                            hSel.string())
        cnt += 1
        name = 'D%d' % cnt
        lScattLength[name] = ((1-val) * lScattLength['H'] + 
                              val * scattLength['D'])
        solventVolume[name] = solventVolume['H'] 
        for atom in hSel:
            groupMap[atom.index()] = name
            if name in atomGroups:
                atomGroups[name] += 1
            else:
                atomGroups[name] = 1
                pass
            pass
        pass
        
    fValues=[]
    fValues_indexed=[]
    gValues=[]
    for qi in range(len(qValues)):
        fValues.append( {} )
        fValues_indexed.append( {} )
        gValues.append( {} )

        for group in list(atomGroups.keys()):
            fValues[qi][group] = lScattLength[group]
            pass
        pass
    return fValues,fValues_indexed,atomGroups,groupMap,excludedVol


metalAtoms=['ZN','MN','MG']

def groupProtons(sel="known"):
    """ return a list of heavy atom/proton groupings.
    """
    ret=[]
    from atomSel import AtomSel
    sel = AtomSel("(%s) and not hydro" % sel.string())
    for atom in sel:
        hsel = AtomSel("bondedto id %d and hydro" % (atom.index()+1))
        if len(hsel):
            entry = [atom]
            for hatom in hsel:
                entry.append( hatom )
                pass
            ret.append(entry)
            pass
        pass
    return ret

from math import pi, exp

from solnScatPotTools import solventVolumeSets
solventVolume = solventVolumeSets['svergun']

exchangableProtons={}
exchangableProtons['GLY']=['HN']
exchangableProtons['ALA']=['HN']
exchangableProtons['VAL']=['HN']
exchangableProtons['LEU']=['HN']
exchangableProtons['ILE']=['HN']
exchangableProtons['PHE']=['HN']
exchangableProtons['TYR']=['HN','HH']
exchangableProtons['TRP']=['HN','HE1']
exchangableProtons['ASP']=['HN']
exchangableProtons['GLU']=['HN']
exchangableProtons['SER']=['HN','HG']
exchangableProtons['THR']=['HN','HG1']
exchangableProtons['ASN']=['HN','HD21','HD22']
exchangableProtons['GLN']=['HN','HE21','HE22']
exchangableProtons['LYS']=['HN','HZ1','HZ2','HZ3']
exchangableProtons['ARG']=['HN','HE','HH21','HH22','HH11','HH12']
exchangableProtons['HIS']=['HN','HD1']
exchangableProtons['MET']=['HN']
exchangableProtons['CYS']=['HN','HG']
exchangableProtons['PRO']=[]
exchangableProtons['MG']=[]

def exchangableName(atom):
    """return True if the atom is a proton exchangable with solvent
    """
    if atom.atomName()=='HN': return 'exchH2'
    if atom.residueName() not in exchangableProtons:
        import sys
        sys.stderr.write('residue name %s not yet supported\n' %
                         atom.residueName());
        return False
    if atom.atomName() in exchangableProtons[atom.residueName()]:
        return 'exchH1'
    return False


scattLength={}  #units are 10^{-12} cm
#values from B. Jacrot, Rep. Prog. Phys. 39, 911-953 (1976)
scattLength['H'] = -.3742
scattLength['D'] =  .6671         
scattLength['C'] =  .6651         #12C
scattLength['N'] =  .940          #14N
scattLength['O'] =  .5804         #16O
scattLength['P'] =  .517          #31P
scattLength['S'] =  .2847         #32S
# values from Cantor and Schimmel
scattLength['Mn'] = -.360
scattLength['Fe'] =  .951  
scattLength['Pt'] =  .95
# from Bauspiess, W., Bonse, U., Rauch H.: Nucl. Instr. Meth. 157 (1978) 495.
# found from http://www.ati.ac.at/~neutropt/scattering/table.html
scattLength['Mg'] = .5375

scattLength['MG'] = scattLength['Mg']


from solnScatPotTools import globRules

def absScattLength(atom):
    """return the abs values of the atom's scattering length
    """
    #FIX: bad way of doing scattLength lookup
    name = atom.atomName()
    if name in metalAtoms:
        name=name
    else:
        name=name[0]
        pass
    xName = exchangableName(atom)
    if xName: name = xName
    return abs(scattLength[name])

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

    solnScatPotTools.useGlobs(term,globTable,globRules,verbose,absScattLength)

    return

