"""Tools to aid in setup/analysis of dipolar coupling potential term.

This module provides functions to simplify the creation, manipulation and
analysis of <m rdcPot>.RDCPot potential terms.

Note on dipolar coupling normalization and sign convention:

  Protein residual dipolars are usually normalized relative to
  measurements of N-HN - this is achieved by calling the function
  scale_toNH.  This normalization involves the gamma_A * gamma_B /
  rAB^3 prefactor (see the docs for <m rdcPot>).  15N has a negative
  gyromagnetic ratio, but it was deemed simpler to treat as positive,
  so that if there is an experiment which does not include 15N, the
  dipolar coupling sign must be flipped.  If you would rather not use
  this convention, call the function correctGyromagneticSigns() at the
  beginning of your script, and then the correct signs will be used
  for the Da prefactors.  """

from cosRatioPot import Cos2RatioPot, CosRatioPot
from atomSel import AtomSel
from avePot import AvePot
from potList import PotList

def create_twoDaRatiosPot(name,
                          rdc1,rdc2,rdc3,rdc4):
    """NOTE: needs to be updated and transferred to <m varTensorTools>.

       Create a potential term which restrains the ratio of
         rdc1.Da/rdc2.Da to be that of rdc3.Da/rdc4.Da
       For EnsembleSimulations, three terms are returned in a PotList. These
       are named 'inter', 'ratio1', and 'ratio2'. The first restrains the
       ensemble average of the two ratios to be the same. The 'ratio' terms
       restrain the ratio of each ensemble member to be the same. The terms
       can be weighted separately using the scale member."""
    from simulation import Simulation_currentSimulation
    from ensembleSimulation import EnsembleSimulation_currentSimulation
    if Simulation_currentSimulation().type() == "EnsembleSimulation":
        ret = PotList(name)
        sim = EnsembleSimulation_currentSimulation().subSim()
        ret.add( AvePot(Cos2RatioPot,
                        "inter",
                        AtomSel("id %d" %(1+rdc1.oAtom().index() ),sim ),
                        AtomSel("id %d" %(1+rdc1.xAtom().index() ),sim ),
                        AtomSel("id %d" %(1+rdc1.p1Atom().index()),sim ),
                        AtomSel("id %d" %(1+rdc2.oAtom().index() ),sim ),
                        AtomSel("id %d" %(1+rdc2.xAtom().index() ),sim ),
                        AtomSel("id %d" %(1+rdc2.p1Atom().index()),sim ),
                        AtomSel("id %d" %(1+rdc3.oAtom().index() ),sim ),
                        AtomSel("id %d" %(1+rdc3.xAtom().index() ),sim ),
                        AtomSel("id %d" %(1+rdc3.p1Atom().index()),sim ),
                        AtomSel("id %d" %(1+rdc4.oAtom().index() ),sim ),
                        AtomSel("id %d" %(1+rdc4.xAtom().index() ),sim ),
                        AtomSel("id %d" %(1+rdc4.p1Atom().index()),sim )) )
        pot2 = CosRatioPot("ratio1",
                           AtomSel("id %d" %(1+rdc1.oAtom().index() ) ),
                           AtomSel("id %d" %(1+rdc1.xAtom().index() ) ),
                           AtomSel("id %d" %(1+rdc1.p1Atom().index()) ),
                           AtomSel("id %d" %(1+rdc2.oAtom().index() ) ),
                           AtomSel("id %d" %(1+rdc2.xAtom().index() ) ),
                           AtomSel("id %d" %(1+rdc2.p1Atom().index()) ))
        pot2.setTargetType( "average" )
        ret.add( pot2 )
        pot3 = CosRatioPot("ratio2",
                           AtomSel("id %d" %(1+rdc3.oAtom().index() ) ),
                           AtomSel("id %d" %(1+rdc3.xAtom().index() ) ),
                           AtomSel("id %d" %(1+rdc3.p1Atom().index()) ),
                           AtomSel("id %d" %(1+rdc4.oAtom().index() ) ),
                           AtomSel("id %d" %(1+rdc4.xAtom().index() ) ),
                           AtomSel("id %d" %(1+rdc4.p1Atom().index()) ))
        pot3.setTargetType( "average" )
        ret.add( pot3 )
        return ret
    else:
        return Cos2RatioPot(name,
                            AtomSel("id %d" %(1+rdc1.oAtom().index())),
                            AtomSel("id %d" %(1+rdc1.xAtom().index())),
                            AtomSel("id %d" %(1+rdc1.p1Atom().index())),
                            AtomSel("id %d" %(1+rdc2.oAtom().index())),
                            AtomSel("id %d" %(1+rdc2.xAtom().index())),
                            AtomSel("id %d" %(1+rdc2.p1Atom().index())),
                            AtomSel("id %d" %(1+rdc3.oAtom().index())),
                            AtomSel("id %d" %(1+rdc3.xAtom().index())),
                            AtomSel("id %d" %(1+rdc3.p1Atom().index())),
                            AtomSel("id %d" %(1+rdc4.oAtom().index())),
                            AtomSel("id %d" %(1+rdc4.xAtom().index())),
                            AtomSel("id %d" %(1+rdc4.p1Atom().index())))
    pass
                        
def create_DaRatioPot(name,
                      rdc1,rdc2,ratio):
    """NOTE: needs to be updated and transfered to <m varTensorTools>.

       create a potential term which restrains the ratio of
         rdc1.Da/rdc2.Da to be ratio
         """
    from simulation import Simulation_currentSimulation
    from ensembleSimulation import EnsembleSimulation_currentSimulation
    if Simulation_currentSimulation().type() == "EnsembleSimulation":
        sim = EnsembleSimulation_currentSimulation().subSim()
        pot =  CosRatioPot(name,
                           AtomSel("id %d" %(1+rdc1.oAtom().index())),
                           AtomSel("id %d" %(1+rdc1.xAtom().index())),
                           AtomSel("id %d" %(1+rdc1.p1Atom().index())),
                           AtomSel("id %d" %(1+rdc2.oAtom().index())),
                           AtomSel("id %d" %(1+rdc2.xAtom().index())),
                           AtomSel("id %d" %(1+rdc2.p1Atom().index())),sim)
        pot.setRatio( ratio )
        ret = AvePot(pot)
    else:
        ret = CosRatioPot(name,
                          AtomSel("id %d" %(1+rdc1.oAtom().index())),
                          AtomSel("id %d" %(1+rdc1.xAtom().index())),
                          AtomSel("id %d" %(1+rdc1.p1Atom().index())),
                          AtomSel("id %d" %(1+rdc2.oAtom().index())),
                          AtomSel("id %d" %(1+rdc2.xAtom().index())),
                          AtomSel("id %d" %(1+rdc2.p1Atom().index())))
        ret.setRatio( ratio )
        pass
    return ret

def create_RDCPot(name,file=0,oTensor=0,esim=0,restraints="",
                  subSel=""):
    """Create an <m rdcPot>.RDCPot1 with given name, given orientational
    <m varTensor>.VarTensor, the filename of an rdc assignment table and/or
    a string of assignments, and an ensemble simulation.

    The file argument can optionally be a sequence of filenames.

    The subSel argument is passed to as the selectionFiler argument to
    <m rdcPot>.RDCPot.addRestraints so that only the subset of restraints
    which match the selection will be loaded.
    """
##    # Alternative docstring. Misses explanation of esim argument.
##    """Return an <m rdcPot>.RDCPot1 instance.
##
##    name is a user-specified identifier for the instance.  file is a string with
##    the filename of an RDC assignment table (or a sequence containing such
##    strings if more than one table is available).  The RDC assignments can also
##    be input directly as a string using the restraints argument.  oTensor is the
##    associated <m varTensor>.VarTensor instance.
##
##    The subSel argument is passed to as the selectionFiler argument to
##    <m rdcPot>.RDCPot.addRestraints so that only the subset of restraints
##    which match the selection will be loaded.
##
##    """
    from rdcPot import RDCPot1
    from simulation import currentSimulation
    from varTensorTools import create_VarTensor
    
    if not esim: esim = currentSimulation()

    needOTensor=False
    if not oTensor:
        needOTensor=True

        #creation of VarTensor will fail if the name already exists:
        # generate a unique name
        otName = "RDC-default-%s"%name
        from simulationTools import getRegisteredTopoTerms
        otTerms=getRegisteredTopoTerms("VarTensor",esim)
        cnt=0
        while otName in [term.instanceName() for term in otTerms]:
            otName = "RDC-default-%s-%d"%(name,cnt)
            cnt += 1
            pass
        oTensor = create_VarTensor(otName)
        pass

    rdc = RDCPot1(name,oTensor,"",esim)

    if needOTensor: rdc.oTensor_ref = oTensor

    if restraints:
        rdc.addRestraints(restraints,selectionFilter=subSel)
        pass

    if file:
        files=[]
        if type(file)==type("string"):
            files.append( file)
        else:
            files=file
            pass
        for file in files:
            rdc.addRestraints( open(file).read(),
                               selectionFilter=subSel)
            pass
        pass
    

    #default averaging type
    rdc.setAveType("sum")

    # default threshold value
    rdc.setThreshold( 0.5 )

    from varTensorTools import registerExptToVarTensor
    registerExptToVarTensor(oTensor,rdc)
    return rdc

def scaleByErrs(rdc):
    """For the specified RDCPot term, set the per-restraint energy
    scale factor to be 1/plusErr^2.
    """
    for r in rdc.restraints():
        r.setScale( 1./r.plusErr()**2 )
        pass
    pass

def correlation(rdc):
    """Return correlation of calculated to observed RDC values.
    """
    
    from math import sqrt
    obs = [r.obs() for r in rdc.restraints()]
    calcd = [r.calcd() for r in rdc.restraints()]
    aveCalcd = sum(calcd) / len(calcd)
    aveObs   = sum(obs) / len(obs)

    aveCalcd2 = sum([(c-aveCalcd)**2 for c in calcd])
    aveObs2  = sum( [(o-aveObs)**2 for o in obs] )
    cross=0.
    for i in range(len(obs)):
        cross += (calcd[i]-aveCalcd) * (obs[i]-aveObs)
        pass
    
    ret = cross / sqrt(aveCalcd2*aveObs2) \
          if abs(aveCalcd2*aveObs2)>1e-30 else -1.
    return ret

def Rfactor(rdc):
    """R-factor (in percent).

    Eq. 5 from Clore+Garrett, JACS 121, 9008 (1999).
    Da is scaled using the current scaling convention.
    """
    from math import sqrt
    return rdc.rms()* 100 / \
           sqrt(2*(rdc.oTensor.aveDa()*rdc.gyroA())**2 *
                (4 + 3 * rdc.oTensor.aveRh()**2)/5)

def Rfactor_infinite(rdc,
                     selection='all'):
    """R-factor (in percent) for an infinite number of
    randomly distributed vectors.
    Eq. 3 from Clore+Garrett, JACS 121, 9008 (1999).

    The selection argument can be used to choose a subset of
    restraints whose atoms lie in selection
    """
    from math import sqrt
    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection,rdc.simulation())

    
    rNum=0
    mean2=0
    dev2=0
    for restraint in rdc.restraints():
        skip=False
        for selPair in restraint.selPairs():
            for aAtom in selPair.a:
                if not selection.containsAtom(aAtom):
                    skip=True
                    continue
                for bAtom in selPair.b:
                    if not selection.containsAtom(bAtom):
                        skip=True
                        continue
                    pass
                pass
            pass
        if not skip:
            mean2 += restraint.obs()**2
            dev2 += (restraint.diff())**2
            rNum += 1
            pass
        pass
    if rNum==0: return 0
    rms = sqrt(dev2/rNum)
    ret = 0
    if rNum>0:
        mean2 /= rNum
        if mean2>0:
            ret = 100*rms/ sqrt(2*mean2)
        pass
    return ret

def composite_Rfactor_infinite(rdcs,
                               selection='all',
                               normalize=True):
    """R-factor (in percent) for an infinite number of
    randomly distributed vectors.
    Eq. 3 from Clore+Garrett, JACS 121, 9008 (1999).

    This version takes a list of <m rdcPot>.RDCPot terms and produces a
    composite R-factor by simply performing sums over all (mixed) RDCs. For
    this to work the RDCs should be unscaled.

    The selection argument can be used to choose a subset of
    restraints whose atoms lie in selection.

    The normalize argument specifies whether or not to rescale RDCs by
    dividing by the Da prefactor to that they are on the same scale
    prior to calculating the R-factor.
    """
    #to do: get subsets working
    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection,rdcs[0].simulation())

    nom2 =0
    denom2 = 0
    for rdc in rdcs:
        for restraint in rdc.restraints():
            skip=False
            for selPair in restraint.selPairs():
                for aAtom in selPair.a:
                    if not selection.containsAtom(aAtom):
                        skip=True
                        continue
                    for bAtom in selPair.b:
                        if not selection.containsAtom(bAtom):
                            skip=True
                            continue
                        pass
                    pass
                pass
            if not skip:
                norm = rdc.gyroA() if normalize else 1.
#                nom2   += (restraint.diff() / rdc.gyroA())**2
#                denom2 += (restraint.obs() / rdc.gyroA())**2
                nom2   += ( restraint.diff() / norm )**2
                denom2 += ( restraint.obs() / norm )  **2
                pass
            pass
        pass
    from math import sqrt
    return 100 / sqrt(2) * sqrt( nom2 / denom2 )    


def chi2(pot,
         selection='all'):
    """Compute the Chi^2 for the restraints associated with the
    specified <m rdcPot>.RDCPot term.

    The selection argument can be used to choose a subset of
    restraints whose atoms lie in selection
    """
    from math import sqrt
    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection,pot.simulation())

    
    rNum=0
    sum2=0
    for restraint in pot.restraints():
        skip=False
        for selPair in restraint.selPairs():
            for aAtom in selPair.a:
                if not selection.containsAtom(aAtom):
                    skip=True
                    continue
                for bAtom in selPair.b:
                    if not selection.containsAtom(bAtom):
                        skip=True
                        continue
                    pass
                pass
            pass
        if not skip:
            err = 0.5 * (restraint.plusErr() + restraint.minusErr())
            if err>0.:
                sum2 += (restraint.diff() / err)**2
                rNum += 1
                pass
            pass
        pass
    if rNum==0: return 0

    chi2 = sum2 /rNum
    return chi2

def composite_chi2(rdcs,
                   selection='all'):
    """Compute the Chi^2 for the restraints associated with the
    specified squence of <m rdcPot>.RDCPot terms.

    The selection argument can be used to choose a subset of
    restraints whose atoms lie in selection
    """
    from math import sqrt
    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection,rdcs[0].simulation())

    
    rNum=0
    sum2=0
    for pot in rdcs:
        for restraint in pot.restraints():
            skip=False
            for selPair in restraint.selPairs():
                for aAtom in selPair.a:
                    if not selection.containsAtom(aAtom):
                        skip=True
                        continue
                    for bAtom in selPair.b:
                        if not selection.containsAtom(bAtom):
                            skip=True
                            continue
                        pass
                    pass
                pass
            if not skip:
                err = 0.5 * (restraint.plusErr() + restraint.minusErr())
                if err>0.:
                    sum2 += (restraint.diff() / err)**2
                    rNum += 1
                    pass
                pass
            pass
        pass
    if rNum==0: return 0

    chi2 = sum2 /rNum
    return chi2
def deviation_percent(rdc):
    """Deviation as a percent of total range of dipolar coupling values.
    Da is scaled using the current scaling convention.
    """
    Da = rdc.oTensor.aveDa()*rdc.gyroA() # using NH convention
    return abs(100*rdc.deviation() /
               ( 3*Da + 1.5*Da* rdc.oTensor.aveRh()))


def orderParameter(aSel,bSel,simulation=0):
    """Calculate the order parameter S^2 given two atom selections defining
    a bond vector. The simulation parameter specifies an EnsembleSimulation
    and defaults to the current simulation. It is assumed that the atomSels
    correspond to single atoms: averaging is not performed"""
    from simulation import Simulation_currentSimulation
    from ensembleSimulation import EnsembleSimulation_currentSimulation
    from vec3 import unitVec, dot
    if not simulation:
        simulation = Simulation_currentSimulation()
        pass
    if simulation.type() != "EnsembleSimulation":
        return 0.
    esim = EnsembleSimulation_currentSimulation()

    bondVecs = [] 
    for index in range( esim.size() ):
        memberSim = esim.members(index)
        aAtom = AtomSel( aSel.string(), memberSim )[0]
        bAtom = AtomSel( bSel.string(), memberSim )[0]
        bondVecs.append( aAtom.pos()-bAtom.pos() )
        pass
    S2 = 0.
    for i in range( esim.size() ):
        S2 += esim.weight(i)**2              # diagonal term
        for j in range(i+1,esim.size() ):
            cosTheta = dot(unitVec(bondVecs[i]),
                           unitVec(bondVecs[j]))
            S2 += esim.weight(i) * esim.weight(j) * \
                  (3*cosTheta**2 - 1)        #two times off-diagonal term
            pass
        pass
    return S2



#def calcTensorSize(rdc,rdcType="NHN"):
#    """
#    NOTE: this needs to be updated and moved to the varTensorTools module.
#
#    given a known structure, and experimental RDCs, and tensor orientation,
#    determine the dipolar coupling Da, and rhombicity. These quantities are
#    returned in a structure.
#    """
#    # dipolar coupling restraints for protein amide NH.
#    import cdsMatrix
#    from cdsVector import CDSVector_double
#    from cdsMatrix import RMat, transpose
#    from vec3 import unitVec
#
#    m = RMat( len(rdc.restraints()) , 5)
#    b = CDSVector_double( len(rdc.restraints()) )
#
#    row=0
#    for restraint in rdc.restraints():
#        for aAtom in restraint.aSelection():
#            for bAtom in restraint.bSelection():
#                uVec = unitVec( bAtom.pos() - aAtom.pos() )
#                m.set(row,0,uVec[1]**2 - uVec[0]**2)
#                m.set(row,1,uVec[2]**2 - uVec[0]**2)
#                m.set(row,2,2 * uVec[0] * uVec[1]  )
#                m.set(row,3,2 * uVec[0] * uVec[2]  )
#                m.set(row,4,2 * uVec[1] * uVec[2]  )
#
#                b[row] = restraint.obs() / DaMax(rdcType)
#
#                if rdc.aveType()=="average":
#                    b[row] *= restraint.aveSize();
#                    pass
#
#                row+=1
#                pass
#            pass
#        pass
#
#    U = []
#    U.append( unitVec( rdc.xAtom().pos() - rdc.oAtom().pos() ))
#    U.append( unitVec( rdc.yAtom().pos() - rdc.oAtom().pos() ))
#    U.append( unitVec( rdc.zAtom().pos() - rdc.oAtom().pos() ))
#    U = transpose( RMat(3,3).fromList(U) )
#
#    A = RMat( len(rdc.restraints()), 2 )
#    row=0
#    for restraint in rdc.restraints():
#        for aAtom in restraint.aSelection():
#            for bAtom in restraint.bSelection():
#                A.set(row,0,
#                      m[row][0]*(U[1][1]**2   -U[1][0]**2   ) +
#                      m[row][1]*(U[1][2]**2   -U[2][0]**2   ) +
#                      m[row][2]*(U[0][1]*U[1][1]-U[0][0]*U[1][0]) +
#                      m[row][3]*(U[0][1]*U[2][1]-U[0][0]*U[2][0]) +
#                      m[row][4]*(U[1][1]*U[2][1]-U[1][0]*U[2][0]))
#                A.set(row,1,
#                      m[row][0]*(U[1][2]**2-U[1][0]**2) +
#                      m[row][1]*(U[2][2]**2-U[2][0]**2) + 
#                      m[row][2]*(U[0][2]*U[1][2]-U[0][0]*U[1][0]) +
#                      m[row][3]*(U[0][2]*U[2][2]-U[0][0]*U[2][0]) +
#                      m[row][4]*(U[1][2]*U[2][2]-U[1][0]*U[2][0]))
#                row+=1
#                pass
#            pass
#        pass
#
#    svdResults = cdsMatrix.svd(A)
#
#    v = transpose( svdResults.vT )
#    uT = transpose( svdResults.u )
#    diag = RMat(2,2,0)
#    print "sigma:",svdResults.sigma[:]
#    for i in range(2):
#        diag.set(i,i, 1.0/svdResults.sigma[i])
#        pass
#
#    # x is [ lambda_2, lambda_3 ]
#    x = v * diag * uT * b
#    values = []
#    values.append( -x[0] - x[1] )
#    values.append( x[0] )
#    values.append( x[1] )
#
#    diag = RMat(3,3,0)
#    diag.set(0,0,values[0])
#    diag.set(1,1,values[1])
#    diag.set(2,2,values[2])
#    tensor = transpose(U) * diag * U
#
#    xInd = 0; yInd=1; zInd=2;
#
#    if abs(values[zInd]) < abs(values[yInd]):
#        zInd,yInd = yInd,zInd
#        pass
#    if abs(values[zInd]) < abs(values[xInd]):
#        zInd,xInd = xInd,zInd
#        pass
#    if abs(values[yInd]) < abs(values[xInd]):
#        yInd,xInd = xInd,yInd
#        pass
#
#    Azz = 2.0/3.0 * values[zInd]
#    Ayy = 2.0/3.0 * values[yInd]
#    Axx = 2.0/3.0 * values[xInd]
#
#    Aa = 3 * Azz / 2.
#    Ar = (Axx-Ayy)
#
#    rdc.setDa( 0.5*DaMax(rdcType) * Aa )
#    rdc.setRhombicity( Ar/Aa )
#
#    return tensor

def getRDCType(rdc):
    """Return the atoms involved in the RDC experiment from the atom selections.  

    """
    import re
    ret = "NHN"
    r = rdc.rawRestraints()
    if not r:
        print("getRDCType: WARNING: no restraints: assuming type:", ret)
        pass
    aName = r[0].selPairs()[0].a[0].atomName()
    bName = r[0].selPairs()[0].b[0].atomName()
    cName = aName+bName
    cTable={}
    cTable['CAC']  = cTable['CCA']  = 'CACO'
    cTable['CACB'] = cTable['CBCA'] = 'CACB'
    cTable['CAHA'] = cTable['HACA'] = 'CAHA'
    cTable['CN']   = cTable['NC']   = 'NCO'
    cTable['NHN']  = cTable['HNN']  = 'NHN'
    cTable['HNC']  = cTable['CHN']  = 'HNCO'

## Edited 11/18/04 by Alex Maltsev ##########################
    cTable['HNCA']  = cTable['CAHN']  = 'CAHN'
#### End correction #############################################

    #nucleic acids
    if re.search(r"C[1-9]'?H[1-9]'?",cName) or \
       re.search(r"H[1-9]'?C[1-9]'?",cName):
        return 'CH'
    if re.search(r"N[1-9]'?H[1-9]'?",cName) or \
       re.search(r"H[1-9]'?N[1-9]'?",cName):
        return 'NH'
    if re.search(r"PH[1-9]'?",cName) or \
       re.search(r"H[1-9]'?P",cName):
        return 'HP'
    if re.search(r"H[1-9]'?H[1-9]'?",cName):
        return 'HH'
    if re.search(r"C5CM",cName): #special treatment of THY methyl groups
        return 'methyl'
    if re.search(r"C5CA",cName): #special treatment of THY methyl groups
        return 'methyl'

    #probably can find better CH params
    if not cName in list(cTable.keys()):
        if ( (aName.startswith("H") and bName.startswith("C")) or
             (aName.startswith("C") and bName.startswith("H")) ):
            cName='CAHA'
            pass
        pass
    
    
    try:
        ret = cTable[cName]
    except:
        raise Exception("getRDCType: could not determine RDC type from: %s"%
                        cName)
    return ret

def makeTable(rdc):
    """Return the assignment table (a string) corresponding to the 
    restraints associated with the specified <m rdcPot>.RDCPot. 

    """
    ret=""
    t = rdc.oTensor
    for restraint in rdc.rawRestraints():
        ret += "assign "
        ret += "( resid %4d and name OO)\n" % t.oAtom().residueNum()
        ret += "\t( resid %4d and name Z )\n"  % t.zAtom().residueNum()
        ret += "\t( resid %4d and name X )\n"  % t.xAtom().residueNum()
        ret += "\t( resid %4d and name Y )\n"  % t.yAtom().residueNum()
        ret += "\t( %s )\n"  % restraint.selPairs()[0].a.string()
        ret += "\t( %s )  "  % restraint.selPairs()[0].b.string()
        ret += "%7.4f %7.4f %7.4f\n\n" % (restraint.obs(),
                                          restraint.minusErr(),
                                          restraint.plusErr())
        for selPair in restraint.selPairs()[1:]:
            ret += "OR\n"
            ret += "\t( %s )\n"  % selPair.a.string()
            ret += "\t( %s )  "  % selPair.b.string()
            pass            
        pass
    return ret
        
        
    

# The following version of makeTable should work, but after coding it I found it
# not very useful.  I'm leaving the original, simpler version to avoid
# confusion. Guillermo A. Bermejo.

##def makeTable(rdc, scale=1.0):
##    """Return the assignment table (a string) associated to an <m rdcPot>.RDCPot.
##
##    rdc is the <m rdcPot>.RDCPot in question.
##
##    The optional argument scale is a number with which to multiply each RDC
##    value in the returned table.  Alternatively, scale can be set to a string
##    with one of the supported nuclei types (e.g., 'NHN'; a list with all the
##    types is given by Da_prefactor.keys()).  If a nuclei pair type is specified,
##    it is used to normalize all RDCs in the returned table.  For example, if
##    scale='NHN', each P-Q RDC will be multiplied by the factor:
##
##    gammaN * gammaH * rNHN**(-3) / gammaP * gammaQ * rPQ**(-3)
##
##    where gammaX is the gyromagnetic ratio of nucleus X, and rXY is the X-Y
##    internuclear distance.
##    
##    """
##    ret=""
##    t = rdc.oTensor
##    
##    if type(scale) is str:
##        scale = Da_prefactor[scale]/Da_prefactor[getRDCType(rdc)]
##        
##    for restraint in rdc.rawRestraints():
##        ret += "assign "
##        ret += "( resid %4d and name OO)\n" % t.oAtom().residueNum()
##        ret += "\t( resid %4d and name Z )\n"  % t.zAtom().residueNum()
##        ret += "\t( resid %4d and name X )\n"  % t.xAtom().residueNum()
##        ret += "\t( resid %4d and name Y )\n"  % t.yAtom().residueNum()
##        ret += "\t( %s )\n"  % restraint.selPairs()[0].a.string()
##        ret += "\t( %s )  "  % restraint.selPairs()[0].b.string()
##        ret += "%7.4f %7.4f %7.4f\n\n" % (restraint.obs()*scale,
##                                          restraint.minusErr(),
##                                          restraint.plusErr())
##        for selPair in restraint.selPairs()[1:]:
##            ret += "OR\n"
##            ret += "\t( %s )\n"  % selPair.a.string()
##            ret += "\t( %s )  "  % selPair.b.string()
##            pass            
##        pass
##    return ret
    


#
# function to relate NH and non-NH expts.
#
def scale_toNH(rdc,type=0):
    """Scale the given expression for dipolar coupling such that
    given Da = scale * Da_NHN

    For thymine, methyl RDCs are detected, and the
    appropriate -3.17 D_CC value [Ottinger and Bax JACS 121, 4690 (1999)] is
    used. If this behavior is not desired please manually specify the type
    argument.

    Also, for proton 'HH' experiments, enables the rdc distance dependence.

    Valid types are:
         'NHN', 'NCO', 'HNCO', 'CACO', 'CACB', 'CAHA', 'HH', 'HP'
    """
    global Da_prefactor
    
    if not type:
        type=getRDCType(rdc)

    scale=0.
    try:
        scale = Da_prefactor[type] / Da_prefactor["NHN"]
    except KeyError:
        print("valid types:", list(Da_prefactor.keys()))
        raise Exception("scale_toNH: no data for coupling type " + type)
    rdc.setGyroA( scale )
    if type=="HH" or type=="HP":
        rdc.setUseDistance(1)
    return

def scale_toCH(rdc,type=0):
    """ scales the given expression for dipolar coupling such that
    given Da = scale * Da_CH

    For thymine, methyl RDCs are detected, and the
    appropriate -3.17 D_CC value [Ottinger and Bax JACS 121, 4690 (1999)] is
    used. If this behavior is not desired please manually specify the type
    argument.

    Also, for proton 'HH' experiments, enables the rdc distance dependence.

    Valid types are:
         'NHN', 'NCO', 'HNCO', 'CACO', 'CACB', 'CAHA', 'HH', 'HP'
    """
    global Da_prefactor

    if not type:
        type=getRDCType(rdc)
    
    scale=0.
    try:
        scale = Da_prefactor[type] / Da_prefactor["CH"]
#        print "type:", type, " scale:", scale
    except KeyError:
        print("valid types:", list(Da_prefactor.keys()))
        raise Exception("scale_toCH: no data for coupling type " + type)
    rdc.setGyroA( scale )
    if type=="HH" or type=="HP":
        rdc.setUseDistance(1)
    return

def DaMax(type,dist=0):
    """Return scale factor to scale other experiments to the NNH norm
     convention.
    """
    global Da_prefactor
    scale=1.0
    try:
        scale *= Da_prefactor[type] / Da_prefactor["NHN"]
    except KeyError:
        raise Exception("DaMax: no data for coupling type " + type)
    if (type=="HH" or type=="HP"):
        if not dist:
            raise Exception("DaMax: distance must be specified.")
        scale /= dist **3
        pass
    return scale

def DaScale(term,dist=0):
    """return scale factor to scale other experiments to the current
    convention, set by calling either scale_toNH or scale_toCH.
    """
    scale = term.gyroA()

    if term.useDistance():
        if not dist:
            raise Exception("DaMax: distance must be specified.")
        scale /= dist **3
        pass
    return scale





#
# constants
#

# gyromagnetic ratios    

gamma_H = 2.675198e4
gamma_N = 2.7116e3    #note that this should be negative: we take this
                      #into account in the sign of Da.
gamma_C = 6.7283e3
gamma_P = 1.0841e4

# canonical bond distances

r_nhn  = 1.042
#r_caha = 1.119
r_caha = 1.105 #CDS 2013/09/16
               #obtained from scales relative to NHN in JACS 122, 10143 (2000)
r_nco  = 1.329
#r_cohn = 2.069
r_cohn = 2.107 #CDS 2013/09/16
               #obtained from scales relative to NHN in JACS 122, 10143 (2000)
r_caco = 1.525
r_cacb = 1.525

## Edited 11/18/04 by Alex Maltsev ##########################
r_cahn = 2.560  # !!!! Value based on an average distance from very
                # few samples; distance from HN(i) to CA(i-1)
#### End correction #############################################

#nucleic acids
# chosen to match factor used in eginput/dna_refi/refine_full.inp
r_nh = r_caha / (gamma_C/gamma_N*0.48125)**(1./3)

def recalculatePrefactors():
    """Recalculate gamma_A * gamma_B / rAB^3 prefactors.

    Call this function after changing gyromagnetic ratios or bond distance
    values. 
    """
    global Da_prefactor
    Da_prefactor = {}
    Da_prefactor["NHN" ] = gamma_N * gamma_H / r_nhn**3
    Da_prefactor["NCO"]  = gamma_N * gamma_C / r_nco**3
    Da_prefactor["HNCO"] = gamma_H * gamma_C / r_cohn**3
    Da_prefactor["CACO"] = gamma_C * gamma_C / r_caco**3
    Da_prefactor["CACB"] = gamma_C * gamma_C / r_cacb**3
    Da_prefactor["CAHA"] = gamma_C * gamma_H / r_caha**3
    Da_prefactor["HH"]   = gamma_H * gamma_H 
    Da_prefactor["HP"]   = gamma_H * gamma_P

    ## Edited 11/18/04 by Alex Maltsev ##########################
    Da_prefactor["CAHN"] = gamma_C * gamma_H / r_cahn**3
    #### End correction #############################################

    #nucleic acids
    Da_prefactor["CH"] = gamma_C * gamma_H / r_caha**3 #FIX: distances!
    Da_prefactor["NH"] = gamma_N * gamma_H / r_nh**3

    #Ottinger and Bax JACS 121, 4690 (1999).
    Da_prefactor["methyl"] = -3.17 * Da_prefactor["CACB"] 
#    Da_prefactor["CC"] = Da_prefactor["CH"]

    return

#initial prefactor setup
recalculatePrefactors()

def correctGyromagneticSigns():
    """Use the correct (negative) sign for 15N's gyromagnetic ratio.

    Call this function before any other in this module.
    """
    global gamma_N
    gamma_N *= -1
    recalculatePrefactors()
    return

headerHelpString=r"""
RMS
  root mean square difference between calculated and observed.
R-fac   
  Eq. 5 from Clore+Garrett, JACS 121, 9008 (1999).  
  Da is scaled using the current scaling convention.       
Devia
  Deviation of calculated values in an EnsembleSimulation calculation.
Da
  Alignment tensor Da.
Rh
  rhombicity of the alignment tensor
R-inf
  R-factor (in percent) for an infinite number of
  randomly distributed vectors.
  Eq. 3 from Clore+Garrett, JACS 121, 9008 (1999)
Viols
  number of RDCs with :math:`|calcd - obs|` values larger than threshold().
Chi^2
  A normalized :math:`\chi^2 =  1/N \sum_i (calcd_i-obs_i)^2 / err_i^2`. It
  is not meaningful if errors are not realistic.
corr
  Pearson's correlation coefficient between calculated and observed RDCs.
"""


def RDC_analyze(potList):
    """Perform analysis of RDCPot terms and return nicely formatted
    summary.

    For each RDC term, the following is printed in the PDB header:

    """

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'RDCPot1')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret+= "%9s  %6s  %6s  %5s  %7s  %5s  %6s %5s %6s %6s\n" % \
     (" " , "RMS", "R-fac", "Devia", "Da", "Rh", "R-inf", "Viols", "Chi^2",
      "corr")

    rdcs = []
    for name in instanceNames:
        rdc = [x for x in potList if x.instanceName()==name][0]
        rdcs.append(rdc)

        print(rdc.showViolations())

        print(rdc.info())

        varTensor = rdc.oTensor
        ret += \
            "%-9s  %6.3f  %6.3f  %5.3f  %7.3f  %5.3f  %6.3f %5d %6.2f %6.3f\n" % \
            (name , rdc.rms(), Rfactor(rdc), deviation_percent(rdc),
             varTensor.Da(), varTensor.Rh(),
             Rfactor_infinite(rdc), rdc.violations(), chi2(rdc), correlation(rdc))
        pass
    
    return ret
RDC_analyze.__doc__ += headerHelpString


from simulationTools import registerTerm, registerExtraStats
registerTerm(RDC_analyze,"Dipolar Coupling Analysis","RDC",
             headerHelpString)

registerExtraStats("RDCPot1","R-fac",Rfactor)
registerExtraStats("RDCPot1","R-infinite",Rfactor_infinite)
registerExtraStats("RDCPot1","Chi^2",chi2)
registerExtraStats("RDCPot1","Correlation",correlation)

def Da(potTerm): return potTerm.oTensor.aveDa()
def Rh(potTerm): return potTerm.oTensor.aveRh()
registerExtraStats("RDCPot1","Da",Da)
registerExtraStats("RDCPot1","Rhombicity",Rh)

def readNEF(potList,block,
            oTensor=None,
            subSel="all",
            minErr=0.1):
    """
    Given a RDCPot object, and a saveframe from a NEF record, read in
    restraints. Each saveframe corresponds to a single <m varTensor>.VarTensor 
    object which is returned. The passed potList argument should be a
    <m potList>.PotList object. Restraints with with the same value of
    scale, weight and distance_dependent are grouped into <m rdcPot>.RDCPot
    objects which are appended to the potList object.
    minErr gives the minimum error value allowed.

    oTensor is an optional <m varTensor>.VarTensor object to use as an alignment
    tensor. If not specified one is created.

    The subSel argument is passed as the selectionFilter argument to
    the RDCPot's addRestraints method.
    """

    rl = block.nef_rdc_restraint_list

    tname = rl.sf_framecode[0][len(rl.sf_category[0])+1:]
    vtens=oTensor
    if not vtens:
        from varTensorTools import create_VarTensor
        vtens = create_VarTensor(tname)
        pass
    try:
        Da = float( rl.tensor_magnitude[0] )
        vtens.setDa(Da)
    except AttributeError:
        pass

    try:
        Rh = float( rl.tensor_rhombicity[0] )
        vtens.setRh(Rh)
    except AttributeError:
        pass

        

    classes = {} #separate RDCPot for each scale/weight/distance dependent

    #TODO: deal with through-space RDCs (distance_dependent true/false)
    #TODO: weight will be interpreted as energy scale

    restraints = block.nef_rdc_restraint
    ids = [int(id) for id in restraints.restraint_id]
    idIndices={}
    for i,id in enumerate(ids):
        if id in idIndices:
            idIndices[id].append(i)
        else:
            idIndices[id] = [i]
            pass
        pass

    from nefTools import fromNefAtomname
    restraintString=""
    if hasattr(restraints,"xplor_nih_use_sign"):
        useSign=None
    else:
        useSign=True
        pass
    for id in sorted(idIndices.keys()):
        first=True
        fromSels=[]
        toSels=[]
        for i in idIndices[id]:

            if useSign==None:
                useSign = restraints.xplor_nih_use_sign[i]
            elif hasattr(restraints,"xplor_nih_use_sign"):
                if useSign != restraints.xplor_nih_use_sign[i]:
                    raise Exception("all entries within a table must have the "+
                                    "same value of xplor_nih_use_sign")
                pass

            segid1=restraints.chain_code_1[i]
            resid1=int(restraints.sequence_code_1[i])
            name1=fromNefAtomname(restraints.atom_name_1[i])
            sel1 = 'resid %d and name %s' % (resid1,name1)
            if segid1!='.': sel1 += ' and segid "%s"' % segid1

            segid2=restraints.chain_code_2[i]
            resid2=int(restraints.sequence_code_2[i])
            name2=fromNefAtomname(restraints.atom_name_2[i])
            sel2 = 'resid %d and name %s' % (resid2,name2)
            if segid2!='.': sel2 += ' and segid "%s"' % segid2

            try:
                scale = restraints.scale[i]
            except AttributeError:
                scale = "1"
                pass
            
            try:
                weight = restraints.weight[i]
            except AttributeError:
                weight = "1"
                pass
            try:
                useDistance = restraints.distance_dependent[i]
            except AttributeError:
                useDistance = "false"
                pass
                
            if not first:
                continue #???
            try:
                val  = float(restraints.target_value[i])
            except ValueError:
                print("Error adding restraint with ID",id)
                continue
            
            try:
                err  = restraints.target_value_uncertainty[i]
            except AttributeError:
                err="."
                pass
            
            try:
                lower = restraints.lower_limit[i]
                upper = restraints.upper_limit[i]
            except:
                lower="."
                upper="."
                pass
            
            if err!='.':
                lower = float(err)
                upper = lower
            elif lower!='.':
                lower = val - float(lower)
                upper = float(upper) - val
            else:
                lower = minErr
                upper = minErr
                pass
            pass

        key = scale+"_"+weight+"_"+useDistance

        if not key in list(classes.keys()):
            name = "%s-%s%s" % (tname,name1,name2)
            while name in [term.instanceName() for term in list(classes.values())]:
                name = name + "1"
                pass
            rdc = create_RDCPot(name,oTensor=vtens)
            rdc.setUseDistance( useDistance.lower()=='true' )
            rdc.setScale( float(weight) )
            rdc.setGyroA( float(scale) )
            classes[key] = rdc
            potList.append(rdc)
        else:
            rdc = classes[key]
            pass


        restraint = """assign (BOGUS)  (BOGUS)  (BOGUS)  (BOGUS) 
                         (%s) (%s) %f %f %f""" % (sel1,sel2,val,lower,upper)
        restraint +=  " ! NEF ID: %d\n" % id
        rdc.addRestraints(restraint,
                          selectionFilter=subSel)
        rdc.setUseSign( useSign )
        pass
    return vtens

def writeNEF_VarTensor(vten,name):
    """
    Return a formatted NEF record for the restraints associated with the
    <m varTensor>.VarTensor object vten as a string with the specified saveframe
    name. 
    """
    from cif import Cif, CifDatablock, CifCategory
    block = CifDatablock()

    cat = CifCategory()
    cat["sf_category"] = "nef_rdc_restraint_list"
    cat["sf_framecode"] = name
    cat["tensor_magnitude"] = str(vten.Da())
    cat["tensor_rhombicity"] = str(vten.Rh())

    block["nef_rdc_restraint_list"] = cat

    cat = CifCategory()
    for key in (
      "restraint_id",
#      "restraint_combination_id",
      "chain_code_1",
      "sequence_code_1",
      "residue_name_1",
      "atom_name_1",
      "chain_code_2",
      "sequence_code_2",
      "residue_name_2",
      "atom_name_2",
      "weight",
      "target_value",
      "xplor_nih_use_sign",
      "lower_limit",
      "upper_limit",
      "scale",
      "distance_dependent",
      ):
        cat.addKey(key)
        pass
    
    from nefTools import toNefAtomname
    id=0
    for rdc in vten.expts():
        if rdc.potName()!="RDCPot1":
            continue
        for r in rdc.restraints():
            id += 1
            obs=r.obs()
            lower = obs - r.minusErr()
            upper = obs + r.plusErr()

            scale = str(rdc.gyroA())
            useDistance = "true" if rdc.useDistance() else "false"
    
            for selPair in r.selPairs():
                for a in selPair.a:
                    for b in selPair.b:
                        cat.addValue("restraint_id",str(id))
                        cat.addValue("chain_code_1",a.segmentName())
                        cat.addValue("sequence_code_1",str(a.residueNum()))
                        cat.addValue("residue_name_1",a.residueName())
                        cat.addValue("atom_name_1",toNefAtomname(a.atomName()))
                        cat.addValue("chain_code_2",b.segmentName())
                        cat.addValue("sequence_code_2",str(b.residueNum()))
                        cat.addValue("residue_name_2",b.residueName())
                        cat.addValue("atom_name_2",toNefAtomname(b.atomName()))
                        cat.addValue("target_value",str(obs))
                        cat.addValue("lower_limit",str(lower))
                        cat.addValue("upper_limit",str(upper))
                        cat.addValue("xplor_nih_use_sign","True" if r.useSign() else "False")
                        cat.addValue("weight",str(rdc.scale()) )
                        cat.addValue("scale",scale)
                        cat.addValue("distance_dependent",useDistance)
                        pass
                    pass
                pass
            pass
        pass
    block["nef_rdc_restraint"] = cat

    cif=Cif()
    block.setIsSaveframe(True)
    cif[name] = block
    cif.setUseTrailingStop(True)
    cif.setUseTrailingSave(True)
    return cif.asString()


def writeNEF(rdc,name):
    """
    Return a formatted NEF record for the restraints in the given RDCPot
    object as a string with the specified saveframe name.
    """

    if rdc.potName()=="VarTensor":
        return writeNEF_VarTensor(rdc,name)
        pass

    #FIX: just call writeNEF_VarTensor(rdc.oTensor)


    from cif import Cif, CifDatablock, CifCategory
    block = CifDatablock()

    cat = CifCategory()
    cat["sf_category"] = "nef_rdc_restraint_list"
    cat["sf_framecode"] = name
    cat["tensor_magnitude"] = '.'#str(11.0000)
    cat["tensor_rhombicity"] = '.'#str(0.0670)

    block["nef_rdc_restraint_list"] = cat

    cat = CifCategory()
    for key in (
      "restraint_id",
#      "restraint_combination_id",
      "chain_code_1",
      "sequence_code_1",
      "residue_name_1",
      "atom_name_1",
      "chain_code_2",
      "sequence_code_2",
      "residue_name_2",
      "atom_name_2",
      "weight",
      "target_value",
      "xplor_nih_use_sign",
      "lower_limit",
      "upper_limit",
      "scale",
      ):
        cat.addKey(key)
        pass
    
    from nefTools import toNefAtomname
    id=0
    for r in rdc.restraints():
        id += 1
        obs=r.obs()
        lower = obs - r.minusErr()
        upper = obs + r.plusErr()

        for selPair in r.selPairs():
            for a in selPair.a:
                for b in selPair.b:
                    cat.addValue("restraint_id",str(id))
                    cat.addValue("chain_code_1",a.segmentName())
                    cat.addValue("sequence_code_1",str(a.residueNum()))
                    cat.addValue("residue_name_1",a.residueName())
                    cat.addValue("atom_name_1",toNefAtomname(a.atomName()))
                    cat.addValue("chain_code_2",b.segmentName())
                    cat.addValue("sequence_code_2",str(b.residueNum()))
                    cat.addValue("residue_name_2",b.residueName())
                    cat.addValue("atom_name_2",toNefAtomname(b.atomName()))
                    cat.addValue("target_value",str(obs))
                    cat.addValue("lower_limit",str(lower))
                    cat.addValue("upper_limit",str(upper))
                    cat.addValue("xplor_nih_use_sign","True" if r.useSign() else "False")
                    cat.addValue("weight","1.0")
                    cat.addValue("scale","1.0")
                    pass
                pass
            pass
        pass
    block["nef_rdc_restraint"] = cat

    cif=Cif()
    block.setIsSaveframe(True)
    cif[name] = block
    cif.setUseTrailingStop(True)
    cif.setUseTrailingSave(True)
    return cif.asString()

def spaceSeparatedToRestraint(inString,
                              atom1Name="N",
                              atom2Name="HN",
                              residCol=1,
                              rdcCol=2,
                              errCol=None,
                              defaultErr=0.1,
                              segid=None,
                              ):
    """
    Convert string restraint table (inString) consisting of columns
    for resid and RDC values to an Xplor-NIH readable restraint
    table. Column numbers start with 1. This function assumes that the RDCs are
    intraresidue with atom names spacified by atom1Name and atom2Name.
    """

    lines = [line for line in inString.split('\n') if not line.startswith('#')]

    ret=""
    for line in lines:
        cols = line.split()
        if not cols:
            continue
        try:
            try:
                resid = int(cols[residCol-1])
                rdc = float(cols[rdcCol-1])
                err = float(cols[errCol-1]) if errCol!=None else defaultErr
            except IndexError:
                print("Warning: could not read line: %s" % line)
                continue
            pass
        except:
            print('error processing line:', line)
            raise
            
        atom1Sel = "name %s and resid %d" % (atom1Name,resid)
        atom2Sel = "name %s and resid %d" % (atom2Name,resid)
        if segid != None:
            atomSel += ' and segid "%s"' % segid
            pass
        ret += "assign () () () ()\n (%s) (%s) %f %f ! %s\n" % (atom1Sel,atom2Sel,
                                                                rdc,err,line)
        pass
    return ret
