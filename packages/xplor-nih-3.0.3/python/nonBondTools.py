"""
 Tools for nonbonded interaction analysis.
"""

from cdsMatrix import CDSMatrix_double
chemTypeDist = CDSMatrix_double(0,0,0.)
from cdsVector import CDSVector_int
lookupVector = CDSVector_int(0)
exclList=[]

from xplorSimulation import getXplorSimulation

def initializeRadii():
    import os
    xSim = getXplorSimulation()
    from simulationTools import mktemp
    tmpFilename=mktemp(suffix=".xplor")
    outputState=xSim.disableOutput()
    xSim.fastCommand("write params output=%s end" % tmpFilename)
    xSim.enableOutput(outputState)
    readXplorRadii(tmpFilename)
    os.unlink(tmpFilename)
    return
def readXplorRadii(filename):
    import re
    global chemTypeDist
    chemTypeRadius=[]
    chemTypeLookup = {}
    for line in open(filename).readlines():
        match = re.search(r"^\s*nonb.*\s+([a-z0-9_+]+)\s+([0-9.]+)\s+([0-9.]+)",
                          line,re.IGNORECASE)
        if match:
            chemType = match.group(1)
            sigma = float( match.group(3) )
            i = len(chemTypeRadius)
            chemTypeRadius.append( sigma )
            chemTypeLookup[chemType] = i
            chemTypeDist.resize(i+1,i+1)
            for j in range(len(chemTypeRadius)):
                chemTypeDist[i,j] = 0.5**(5./6) * (chemTypeRadius[i]+
                                                   chemTypeRadius[j])
                chemTypeDist[j,i] = chemTypeDist[i,j]
                pass
            pass

        match = re.search(
            r"^\s*nbfi.*\s+([a-z0-9_+]+)\s+([a-z0-9_+]+)\s+([0-9.]+)\s+([0-9.]+)",
                          line,re.IGNORECASE)
        if match:
            i = chemTypeLookup[ match.group(1) ]
            j = chemTypeLookup[ match.group(2) ]
            A = float( match.group(3) )
            B = float( match.group(4) )
            sigma=0
            if B>0: sigma = (A/B)**(1./6)
            chemTypeDist[i,j] = 0.5**(5./6) * sigma
            chemTypeDist[j,i] = chemTypeDist[i,j]
            pass
        pass

    xSim = getXplorSimulation()
    lookupVector.resize( xSim.numAtoms() )
    for i in range(xSim.numAtoms()):
        lookupVector[i] = chemTypeLookup[ xSim.atomByID(i).chemType() ]
        pass

    global exclList
    exclList=[]
    bondList=[]
    for i in range(xSim.numAtoms()):
        exclList.append([])
        bondList.append([])
        pass

    outputState=xSim.disableOutput()
    xSim.fastCommand("set print none end")
    nbondNBXMod = int(xSim.fastCommand("param nbond ? end end","NBXMOD")[0])
    xSim.fastCommand("set print $prev_print_file end")
    xSim.enableOutput(outputState)
    if abs(nbondNBXMod)>1:
        for i in range(xSim.numBonds()):
            pair = xSim.bondPairByID(i)
            exclList[pair[0]].append( pair[1] )
            exclList[pair[1]].append( pair[0] )
            bondList[pair[0]].append( pair[1] )
            bondList[pair[1]].append( pair[0] )
            pass
        pass
    #add in atoms with 1-3 relationship
    from nbond import addExclusions
    if abs(nbondNBXMod)>2:
        exclList = addExclusions(exclList,bondList)
        pass
    #add in atoms with 1-4 relationship
    if abs(nbondNBXMod)==4:
        exclList = addExclusions(exclList,bondList)
        pass
    #add in explicitly excluded interactions
    if nbondNBXMod>0:
        #write out psf
        outputState=xSim.disableOutput()
        from simulationTools import mktemp
        tmpFilename=mktemp(suffix=".xplor")
        xSim.fastCommand("write psf output=%s end" % tmpFilename)
        xSim.enableOutput(outputState)
        # read in nbondexcl list
        psfFile=open(tmpFilename)
        psfFile.readline();        psfFile.readline()
        numTitle=int(psfFile.readline().split()[0])
        for i in range(numTitle):
            psfFile.readline()
            pass
        psf=psfFile.read()
        (before,after)=psf.split('!NNB')
        numNNB=int(before.split()[-1])
        combTable = [int(e) for e in after.split('!NGRP')[0].split()][:-2]
        (nnbTable,idxTable)=(combTable[:numNNB],combTable[numNNB:])

        #update exclList
        idx=0
        for i in range(xSim.numAtoms()):
            if idx!=idxTable[i]:
                for j in range(idx,idxTable[i]):
                    exclList[i].append(nnbTable[j]-1)
                    pass
                pass
            idx=idxTable[i]
            pass
        pass
                    
        
        
    return

defaultRadiusScale=0.8

def vdwViolations(threshold,
                  selection="not PSEUDO",
                  selection2=None,
                  radiusScale=None):
    """determine nbonded violations:
        atoms pairs in selection for which threshold plus the distance
        between the atom positions is less than the sums of the appropriate,
        scaled vdw radii.

      The scale factor is taken from the XPLOR REPEl parameter, if it is used.
      Otherwise the variable defaultRadiusScale is used (default value is 0.8).

      The XPLOR NBXMod parameter is consulted and the appropriate 1-2,3,4
      interactions are excluded. Also, explicitly excluded interactions are
      excluded, if NBXMod>0.

      If selection2 is specified, only distances between atoms in the two
      selections are considered.

      Radius scale is normally set from the value of the XPLOR REPEL
      parameter. An alternative value can be specified using the radiusScale
      argument.
    """
    xSim = getXplorSimulation()
    if len(lookupVector) != xSim.numAtoms():
        initializeRadii()
        pass

    from selectTools import convertToAtomSel
    atoms1 = convertToAtomSel(selection)
    if selection2!=None:
        atoms2=convertToAtomSel(selection2)
    else:
        atoms2=atoms1
        pass

    if radiusScale==None:
        outputState=xSim.disableOutput()
        xSim.fastCommand("set print none end")
        try:
            radiusScale = float(xSim.fastCommand("param nbond ? end end",
                                                 "repel")[0])
        except:
            radiusScale = defaultRadiusScale
            pass
        xSim.fastCommand("set print $prev_print_file end")
        xSim.enableOutput(outputState)
        pass
    
    distanceMatrix = CDSMatrix_double(chemTypeDist)
    distanceMatrix.scale( radiusScale )

    violations=[]
    from vec3 import norm
    from nbond import findCloseAtoms2
    for (i,j) in findCloseAtoms2(atoms1,atoms2,lookupVector,
                                 distanceMatrix, threshold):
        if not j in exclList[i]:
            atomi = atoms1[i]
            atomj = atoms2[j]
            violations.append(
                Violation(atomi, atomj,
                          norm(atomi.pos()-atomj.pos()),
                          distanceMatrix[lookupVector[atomi.index()],
                                         lookupVector[atomj.index()]]) )
            pass
        pass

    return violations

#def vdwViolationsPrint(threshold):
#    """print out vdw violation info, and return the number of violations.
#    """
#    violations = vdwViolations(threshold)
#
#    print
#    print "  Nonbonded Violations"
#    print "        threshold: %f" % threshold
#    print
#    print "%26s %26s %6s %6s" % ("atom1     ", "atom2     ",
#                                 "dist", "vdwDist")
#    print "_"*70
#
#    for v in violations:
#        print "%26s %26s %6.2f %6.2f" % (v.atomi.string(), v.atomj.string(),
#                                         v.dist, v.dist0)
#        pass
#        
#    return len(violations)

class Violation:
    def __init__(s,atomi,atomj,dist,dist0):
        s.atomi = atomi
        s.atomj = atomj
        s.dist  = dist
        s.dist0 = dist0
        return
    def name(s):
        return "( %s ) ( %s )" % (s.atomi.string(),s.atomj.string())
    def diff(s):
        return s.dist - s.dist0
    
        
def getNBParamsFromXplor(simulation):
    """
    Read nonbonded parameters from underlying XplorSimulation, via the PSF
    and loaded parameters.
    """

    from xplorSimulation import getXplorSimulation
    from simulationTools import mktemp
    xSim = getXplorSimulation(simulation)


    outputState=xSim.disableOutput()
    psfFilename=mktemp('xplor-psf')
    #xSim.fastCommand("set print %s end" % psfFilename)
    xSim.fastCommand("write psf output=%s end"%psfFilename )
    #xSim.fastCommand("set print $prev_print_file end" )

    from collections import namedtuple
    Ret = namedtuple('getNBParamsFromXplorRet',
                     ['nbExclusions','nbParams', 'nbfixParams',
                      'startIndices','resExclude'])


    lines=open(psfFilename).readlines()
    processGroups=False
    processExclusions=False
    processExclIndices=False
    exclusions=[]
    exclIndices=[]
    groupsString=""
    for line in lines:
        line = line.rstrip()
        if line.endswith(' !NATOM'):
            psfAtoms = int(line.split()[0])
        elif line.endswith(' !NNB'):
            numExclusions = int(line.split()[0])
            processExclusions=True
        elif processExclusions:
            lineExclusions=[int(i) for i in line.split()]
            exclusions += lineExclusions
            if len(exclusions)==numExclusions:
                processExclusions=False
                processExclIndices=True
                pass
            pass
        elif processExclIndices:
            lineIndices=[int(i) for i in line.split()]
            exclIndices += lineIndices
            if len(exclIndices)==psfAtoms:
                processExclIndices=False
                pass
        elif processGroups:
            groupsString+= line
        elif line.endswith(' !NGRP'):
            processGroups=True
            pass
        pass

    # probably should directly generate a CDSList< CDSList<int> >
    # nbExclusions is a list of size numAtoms
    # each entry is a list of exclusions for that atom
    nbExclusions=[[]]
    for i in range(psfAtoms-1):
        nbExclusions.append( [] )
        for j in range(exclIndices[i],exclIndices[i+1]):
            nbExclusions[i+1].append( exclusions[j]-1 )
            pass
        pass
    nbExclusions.append([])
    if exclIndices:
        for j in range(exclIndices[-1],numExclusions):
            nbExclusions[-1].append( exclusions[j]-1 )
            pass
        pass
    groups = groupsString.split()

    startIndices=[]
    for cnt,startIndex in enumerate(groups):
        if cnt%3!=0:
            continue
        startIndices.append(int(startIndex))
        pass

##=======initialize vdw params
    import os
    parFilename=mktemp('param-output')
    xSim.fastCommand("write param output=%s end" % parFilename)

    xSim.enableOutput(outputState)
    paramOutput = [line.strip() for line in open(parFilename).readlines()]
    os.unlink(parFilename)

    nonbondLines=[line for line in paramOutput if line.startswith('NONBon')]

    nbParams={}
    for line in nonbondLines:
        contents=line.split()
        chemtype=contents[1]
        nbParams[chemtype] = [float(v) for v in contents[2:6]]
        pass
    
    nbfixLines=[line for line in paramOutput if line.startswith('NBFIx')]
    nbfixParams={}
    for line in nbfixLines:
        contents=line.split()
        chemtypes="*".join(contents[1:3])
        nbfixParams[chemtypes] = [float(v) for v in contents[3:]]
        pass

    #go through nbExclusions
    # for each, id residue name and atom name
    # for each exclusion add (atomi, atomj) tuple, if it doesn't already exist.
    # also record resid, if it changes can just ignore the exclusions
    from collections import defaultdict
    resExclude=defaultdict(list)
    resNameNumExclMap={}
    from atom import Atom

    for indexi,excl in enumerate( nbExclusions ):
        if not excl: continue
        atomi = Atom(xSim,indexi)
        namei = atomi.atomName()
        resname = atomi.residueName()
        if (resname in list(resNameNumExclMap.keys()) and
            resNameNumExclMap[resname] != atomi.residueNum()):
            continue
        else:
            resNameNumExclMap[resname] = atomi.residueNum()
            pass
        for indexj in excl:
            atomj = Atom(xSim,indexj)
            namej = atomj.atomName()
            resExclude[resname].append( (namei,namej) )
            pass
        pass

    # special case for group definitions for SymSimulations
    if simulation.type()=="SymSimulation":
        from symSimulation import fromSimulation
        symSim = fromSimulation(simulation)
        xIndices = startIndices
        startIndices = []
        totCopies=1
        # for nested SymSimulations: assume that the selections
        # each reference the same set of XplorSimulation atoms
        while True:
            sel = symSim.subSel() #the simulation should be or mirror xSim
            totCopies *= symSim.numCopies()
            if sel.simulation().type() != "SymSimulation":
                break
            symSim = fromSimulation(sel.simulation())
            pass
        for n in range(totCopies):
            startIndex = sel[0].index()
            selSize = len(sel)
            for index in xIndices:
                if index>=startIndex and index<sel[-1].index():
                    startIndices.append( index + n * selSize - startIndex )
                    pass
                pass
            pass
        pass

    return Ret(nbExclusions,nbParams, nbfixParams, startIndices, resExclude)

  #  for i in startIndices:
  #      print i
  
