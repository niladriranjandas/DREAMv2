"""
Functions to help in creating and manipulating AnglePots
"""

def findAngles(centerAtom):
    """
    Find all bond angles with the specified central atom.
    """
    from selectTools import convertToAtom
    centerAtom = convertToAtom(centerAtom)
    
    centerID = centerAtom.index()+1

    from xplorSimulation import getXplorSimulation
    from simulationTools import mktemp
    sim = centerAtom.simulation()
    xSim = getXplorSimulation(sim)


    outputState=xSim.disableOutput()
    psfFilename=mktemp('xplor-psf')
    #xSim.fastCommand("set print %s end" % psfFilename)
    xSim.fastCommand("write psf output=%s end"%psfFilename )
    #xSim.fastCommand("set print $prev_print_file end" )


    lines=open(psfFilename).readlines()
    processAngles=False
    angleIDlist=[]
    numSeen=0
    for line in lines:
        #        line = line.rstrip()
        if "!" in line:
            if line.split("!")[1].startswith('NATOM'):
                psfAtoms = int(line.split()[0])
            elif line.split("!")[1].startswith('NTHETA'):
                numAngles = int(line.split()[0])
                processAngles=True
                pass
            pass
        elif processAngles:
            ids=[int(i) for i in line.split()]
            while ids:
                angleIDs = ids[:3]
                ids=ids[3:]
                numSeen+=1
                if angleIDs[1]==centerID:
                    angleIDlist.append( angleIDs )
                    pass
                pass
            if numSeen==numAngles:
                break
            pass
        pass

    

    print("numAtoms:", psfAtoms)
    print("centerID:", centerID, centerAtom.atomName())
    
    print("numAngles:", numAngles)

    print("matching angles:", len(angleIDlist))
    #    print angleIDlist

    #get parameters
    import os
    parFilename=mktemp('param-output')
    xSim.fastCommand("write param output=%s end" % parFilename)

    xSim.enableOutput(outputState)
    paramOutput = [line.strip() for line in open(parFilename).readlines()]
    os.unlink(parFilename)
    
    #FIX: multiline angle statements will cause problems
    angleLines=[line for line in paramOutput if line.startswith('ANGLe')]
    from potUtils import stripBracketedComments
    angleLines=[stripBracketedComments(line) for line in angleLines]

    #type-based angle lookup
    angleParams=[]    
    for line in angleLines:
        if "(" in line:
            print("findAngles: Warning: found an atom-based angle parameter:")
            print("  ",line)
            continue
        a,b,c,forcec,theta0 = line.split()[1:]
        angleParams.append( (a,b,c,forcec,theta0) )
        pass
    

    from atomSel import AtomSel
    angleAtomList=[]
    for (a,b,c) in angleIDlist:
        angleAtomList.append( (AtomSel("id %d"%a,sim)[0],
                               AtomSel("id %d"%b,sim)[0],
                               AtomSel("id %d"%c,sim)[0]) )
        pass
    paramList=[]
    for (aType,bType,cType) in [(a.chemType(),b.chemType(),c.chemType())
                                for a,b,c in angleAtomList]:
        found = False
        matchParams = [params for params in angleParams
                       if params[1]==bType]
        for (a,b,c,forcec,theta0) in matchParams:
            if (a==aType and c==cType) or (a==cType and c==aType):
                paramList.append((float(forcec),float(theta0)))
                found=True
                break
            pass
        if not found:
            raise Exception("parameters not found for angle types:" +
                            "%s %s %s" % (aType,bType,cType))
        pass
    
#    print angleLines
    

    ret=[]
    for i,angleAtoms in enumerate(angleAtomList):
        ret.append( list(angleAtoms) + list(paramList[i]) )
        pass

    return ret

def deleteAngles(angles):
    """
    Given a list of atom triples, delete the associated bond angles
    from the topology.

    Side effect: all topology entries are erased.
    """
    if not angles: return

    sim=angles[0][0].simulation()
    import xplor
    runx = xplor.command
    for (a,b,c,forcec,theta0) in angles:
        patch = """
        topo
        reset
        presidue temp
        delete  angle +%s   +%s   +%s 
        end
        end""" % (a.atomName(),
                  b.atomName(),
                  c.atomName())
        runx(patch)
        runx('patch temp reference=+=(resid %d and segid "%s") end' %
             (a.residueNum(),a.segmentName()))
        pass
    pass

def delAnglesFromCenters(centerNames=[],
                         centerSels =[]):
    """Given a list of atoms, find all XPLOR bond angles from the bonding
    pattern, and delete them. The atoms can be specified as a list of atom
    names or atom selections using the centerNames or centerSels arguments,
    respectively.
    """
    from selectTools import convertToAtomSel
    centerSels = [convertToAtomSel(sel) for sel in centerSels]
    from atomSel import AtomSel
    centerSels += [AtomSel("name " + name)  for name in centerNames]
    angleList = [findAngles(sel) for sel in centerSels]

    for angle in angleList:
        deleteAngles(angle)
        pass
    return

def delAndCreateCenterAnglePots(instanceName,
                                centerNames=[],
                                centerSels =[]):
    """Given a list of atoms, delete XPLOR angle terms and add bond angle
    energy terms for all angles for which the atoms are the center atom,
    determined from the bonding pattern. The atoms can be specified as a list
    of atom names or atom selections using the centerNames or centerSels arguments, respectively.
    """

    from potList import PotList
    ret = PotList(instanceName)

    from selectTools import convertToAtomSel
    centerSels = [convertToAtomSel(sel) for sel in centerSels]
    from atomSel import AtomSel
    centerSels += [AtomSel("name " + name)  for name in centerNames]


    angleList = [findAngles(sel) for sel in centerSels]

    from anglePot import AnglePot
    for i,angles in enumerate(angleList):
        ret.append( AnglePot("%s-%s" % (instanceName,
                                        centerSels[i][0].atomName()),
                                        angles) )
        pass

    delAnglesFromCenters( centerSels=centerSels )

    return ret


