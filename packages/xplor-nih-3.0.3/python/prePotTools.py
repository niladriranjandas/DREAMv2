"""Tools to aid in setup/analysis of the Paramagnetic Relaxation Enhancement
potential term.

This module provides functions to simplify the creation, manipulation and
analysis of <m prePot>.PREPot potential terms.



"""

import protocol
protocol.addPseudoResName("TAU")

def create_PREPot(name, file=0,
                  assignType="",
                  frequency=500,
                  scaleType="obsig",
                  esim=0, restraints="",
                  eSpinQuantumNumber=2.5,
                  rlxType="r2dd",
                  tauc=3.0, taut=0.5, taui=0.05,
                  fixTau=False,
                  clockResid=None,
                  clockSegid=None,
                  paraCenter=None,
                  nucleusID=None,
                  sumWeights=[],
                  verbose=True,
                  ):
    """Create a <m prePot>.PREPot with given name, given the filename of a
    PRE assignment table and/or a string of assignments, and an
    ensemble simulation.

    The assignType argument is a a deprecated argument, and now has no effect.

    Specify the NMR Larmor frequency (in MHz) using the frequency argument.

    eSpinQuantumNumber specifies the electron spin quantum number of the
    spin label.

    rlxType is a string that specifies the type of paramagnetic relaxation
    mechanism; it can be one of 'r2dd', 'r2curie', 'r2mix', 'r1dd', 'r1curie',
    or 'r1mix'.
    
    tauc, taut, and taui specify correlation times in ns.  By default, tauc is
    optimized along with atomic coordinates in the range 0.5tc .. 2tc, where
    tc is the specified value.
    
    Changing fixTau to True will disable correlation-time optimization (clock
    pseudo atoms will not be added).

    clockResid and clockSegid can be set if clock pseudo atoms preexist.

    nucleusID can be one of  '15N', '1H', '13C', or '31P'.  If it is omitted,
    an attempt is made to determine this type automatically.

    The input restraint list is rewritten to be sum restraints if paraCenter is
    set - it specifies which selection is to be treated as the paramagnetic
    center.  For each restraint in the input table a SUM restraint is added for
    each atom which matches paraCenter.  By default the SUM terms will be given
    equal weights of 1, but the sumWeights argument can be used to change these
    values.  The SUM restraints will be in order of appearance in the PSF.
    
    """
    from prePot import PREPot
    from simulation import currentSimulation
    
    if not esim: esim = currentSimulation()

    if file: restraints += open(file).read()

    pre = PREPot(name,restraints,assignType,esim)

    if paraCenter:
        pre.mkSumRestraints(paraCenter,sumWeights)
        pass

    pre.setFreqI( frequency )

    nType=nucleusID if nucleusID else getNucleusType(pre)
    if verbose:
        print("create_PREPot: detected %s nuclei" % nType)
    
    

    #default values  -- WARNING: should match those in prePot module
    gammaI=26.752196
    gFac=2.0

    pre.setAveType("r-6")
    pre.setSclType(scaleType)
    pre.setRlxType(rlxType)
    pre.setGammaI(gyroMagneticRatios[nType])
    pre.setGfac( gFac )
    pre.setShowAllRestraints(True)
    pre.setSqn(eSpinQuantumNumber)

    if not fixTau:
        (tco,tcx,tcy,
         tto,ttx,tty,
         tio,tix,tiy) = addClockAtoms(clockResid,clockSegid,
                                      simulation=pre.simulation())

        pre.simulation().sync()
        pre.setTaucAtoms(tcx,tcy,tco)
        pre.setFixTauc(False)
        pre.setTcMin(0.5*tauc)
        pre.setTcMax(2*tauc)
        pre.setTautAtoms(ttx,tty,tto) # taut is in range 0 < tauc
        pre.setTauiAtoms(tix,tiy,tio)
        pre.setTiMin(0.5*taui)
        pre.setTiMax(2*taui)
        pass

    pre.setTauc(tauc)
    pre.setTaut(taut)
    pre.setTaui(taui)

    registerTerm(pre)

    return pre

gyroMagneticRatios={} # in  1/s 1/gauss
gyroMagneticRatios['1H']  = 26.75198
gyroMagneticRatios['15N'] =  2.7116 #sign doesn't matter here
gyroMagneticRatios['13C'] =  6.7283
gyroMagneticRatios['31P'] = 10.841

def getNucleusType(pre):
    """Determine the nucleous involved in the experiment.
    """
    import re
    ret = "1H"
    restraints = pre.restraints()
    if not restraints:
        print("getNucleusType: WARNING: no restraints: assuming type:", ret)
        return ret
        pass
    ret=''

    def isN(atom):
        return (atom and atom.atomName().startswith('N') and atom.mass() < 15)
    def isC(atom):
        return (atom and atom.atomName().startswith('C') and atom.mass() < 13)
    def isH(atom):
        return (atom and atom.atomName().startswith('H') and atom.mass() < 2)
    def isP(atom):
        return (atom and atom.atomName().startswith('P') and atom.mass() < 32)

    atoma=restraints[0].sel1[0]
    atomb=restraints[0].sel2[0]

    import selectTools
    if not atoma.residueName() in list(selectTools.residueMapProtein.values()):
        atoma=None
        pass

    if not atomb.residueName() in list(selectTools.residueMapProtein.values()):
        atomb=None
        pass


    numAssigned=0
    if (isN(atoma) or isN(atomb)):
        numAssigned += 1
        ret= '15N'

    if (isH(atoma) or isH(atomb)):
        numAssigned += 1
        ret= '1H'
        
    if (isC(atoma) or isC(atomb)):
        numAssigned += 1
        ret= '13C'

    if (isP(atoma) or isP(atomb)):
        numAssigned += 1
        ret= '31P'
        pass

    if numAssigned>1:
        raise Exception("more than one assignment type. Make sure you change "
                        "masses (via massSetup() after calling create_PREPot")

    if not ret:
        print("getNucleusType: Could not determine nucleus type.")
        print("  Please manually specify gyromagnetic ratio using setGammaI.")
        ret = "1H"
        pass
        
    return ret

def getTagResidSels(pre):
    """
    Generate a list of <m atomSel>.AtomSels for each unique tag residue. This
    can be used to conveniently configure selection pairs for non-bonded
    interactions.

    Here it is assumed that sel1 in each restraint refers to the paramagnetic
    center.

    This function returns a named tuple with sels and union members, with the
    second argument the union of the list.
    """
    from atomSel import AtomSel
    tagSel=AtomSel("")
    from atomSel import AtomSel, union
    for r in pre.restraints():
        # for now, assume the tag is specified by sel1
        tagSel = union(tagSel, r.sel1 )
        pass
    
    tagResidueSel=AtomSel("")
    for tagAtom in tagSel:
        tagResidueSel = union(tagResidueSel,
                       AtomSel('segid "%s" and resid %d' % (tagAtom.segmentName(),
                                                            tagAtom.residueNum())))
        pass
    from selectTools import getSegsResidues
    segDict = getSegsResidues(tagResidueSel)
    sels=[]
    for segid in list(segDict.keys()):
        sels += [ AtomSel('segid "%s" and resid %d' % (segid, tuple[0])) for
                  tuple in segDict[segid] ]
        pass

    from collections import namedtuple
    Ret = namedtuple('getTagResidSels',['sels','union'])
    return Ret(sels,tagResidueSel)

    

def setupSBmode(term):
    """ Configure the given potential term to use Solomon-Bloembergen equation
    to calculate PREs.
    """
    term.setEquType("sb")
    return

def setupSBMFmode(term):
    """ Configure the given potential term to use modified model-free
    Solomon-Bloembergen equation to calculate PREs.
    """
    term.setEquType("sbmf")
    term.setSbmfType("taut")
    term.setFixTauc(False)
    term.setFixTaui(True)
    term.setFixTaut(True)
    term.setTauT(0.00)
    term.setTauI(0.00)
    return

def Qfactor(term):
    """Term can be a single potential term, or a list of terms. If a list,
    the combined Q-factor is computed, not the average of the terms.
    """
    from simulationTools import potType
    if potType(term)=='PRE':
        q=term.qFactor()
    else:
        # assume we have a list of terms
        from simulationTools import flattenPotList
        prelist = flattenPotList(term)
        
        from math import sqrt

        snu = 0.0 ; sde = 0.0
        for p in prelist:
            rn=sum([1 if r.isStereoAssigned() else 2 for r in p.restraints()])
            snu += p.rms() * p.rms() * rn
            sde += p.AveSqObs(1) * rn
            pass
        q = sqrt(snu / sde)
        pass
    return q

def correlation(term):
    """Term can be a single potential term, or a list of terms, in which case
    the average correlation is returned. FIX: should return overall
    correlation instead.
    """
    from simulationTools import potType
    if potType(term)=='PRE':
        return term.correlation()
    else:
        # assume we have a list of terms
        from simulationTools import flattenPotList
        prelist = flattenPotList(term)

        if len(prelist)==0:
            return 0.

        sum=0
        for p in prelist:
            sum += p.correlation()
            pass

        return sum / len(prelist)
    return 

default_segid=""
default_resid=1200
def addClockAtoms(resid=None,segid=None,
                  tcOatom=None,tcXatom=None,tcYatom=None,
                  ttOatom=None,ttXatom=None,ttYatom=None,
                  tiOatom=None,tiXatom=None,tiYatom=None,
                  simulation=None                                    
                  ):
    """
    Create psf and initial coordinates for all clock atoms used for optimizing
    correlation times tau_c, tau_t and tau_i.
    
    """

    global default_segid, default_resid
    if not segid:
        segid = default_segid
        
    from ensembleSimulation import EnsembleSimulation_currentSimulation
    esim = simulation if simulation else EnsembleSimulation_currentSimulation()
    if not esim:
        from ensembleSimulation import EnsembleSimulation_sizeOneSimulation
#        from xplorSimulation import getXplorSimulation
        from simulation import currentSimulation
        sim = currentSimulation()
        esim = EnsembleSimulation_sizeOneSimulation(sim)
        pass

    mem0 = esim.members(0)

    from xplorSimulation import getXplorSimulation
    xSim = getXplorSimulation(mem0)

    # EnsembleSimulations and SymSimulations don't mix (yet)
    sim = mem0 if esim.size()>1 else xSim

    import protocol
    protocol.initParams("axis")

    from atomSel import AtomSel
    print('addClockAtoms: ', end=' ')
    if (resid!=None and 
        len(AtomSel('segid "%s" and resid %d' % (segid,resid),sim))):
        print('using existing atoms in segid: "%s" resid: %d' % (segid,resid))
    else:
        if resid==None: resid=default_resid
        
        while len( AtomSel('segid "%s" and resid %d' % (segid,resid),sim))>0:
            resid +=1 
            pass
        
        if len( AtomSel('segid "%s" and resid %s' % (segid,resid),sim)) > 0:
            raise Exception('segid "%s" and resid %d already exist.' %
                            (segid,resid))
        
        if esim.singleThread():
            print('adding clock atoms in segid: "%s" resid: %d' % (segid,resid))
            
            cmd = psfTemplate.replace('_n__','%-4d'%resid)
            cmd = cmd.replace('SGMT','%-4s'%segid)
            xSim.fastCommand(cmd)
            xSim.syncFrom()
            pass
        esim.multiThread()


        pass
    esim.sync()

    protocol.updatePseudoAtoms(sim)

    from vec3 import unitVec, dot, cross, Vec3, norm
    from selectTools import convertToAtom
#    if oAtom:  FIX: can't use existing clock atoms w/ this routine
#        oAtom = convertToAtom(oAtom)
#        xAtom = convertToAtom(xAtom)
#        yAtom = convertToAtom(yAtom)
#
#        p1=oAtom.pos()
#        p2=xAtom.pos()
#        p3=yAtom.pos()
#        #FIX: angle?
#        pass

    mem0 = esim.members(0)
    selPrefix = 'segid "%s" and resid %d and name ' %(segid,resid)
    tcOatom, tcXatom, tcYatom = [AtomSel(selPrefix + name, sim) for name in
                                 ('TCO', 'TCA', 'TCB')]
    ttOatom, ttXatom, ttYatom = [AtomSel(selPrefix + name, sim) for name in
                                 ('TTO', 'TTA', 'TTB')]
    tiOatom, tiXatom, tiYatom = [AtomSel(selPrefix + name, sim) for name in
                                 ('TIO', 'TIA', 'TIB')]

    if esim.singleThread():
        from atomAction import centerOfMass
        pCM = centerOfMass("known and not PSEUDO")
        pO = pCM + Vec3(20,20,20)
        pX = pO + Vec3(1,0,0)
        pY = pO + Vec3(0,1,0)

        resSel = 'segid "%s" and resid %d' % (segid,resid)

        for (oname,xname,yname) in (("TCO","TCA","TCB"),
                                    ("TTO","TTA","TTB"),
                                    ("TIO","TIA","TIB")):
            oAtom = AtomSel(resSel + " and name "+oname ,sim)[0]
            xAtom = AtomSel(resSel + " and name "+xname ,sim)[0]
            yAtom = AtomSel(resSel + " and name "+yname ,sim)[0]
            oAtom.setPos( pO )
            xAtom.setPos( pX )
            yAtom.setPos( pY )
            pass
        pass
    esim.multiThread()

    return (tcOatom,tcXatom,tcYatom,
            ttOatom,ttXatom,ttYatom,
            tiOatom,tiXatom,tiYatom)

# the string _n__ is replaced by the residue number
# the string SGMT is replaced by the segment name
psfTemplate = """
structure
PSF

       1 !NTITLE
 REMARKS   prePotTools.py: auto-generated structure parameters

       9 !NATOM
       1 SGMT _n__ TAU  TCO  OOO    0.000000E+00   10.0000           0
       2 SGMT _n__ TAU  TCA  XXX    0.000000E+00   10.0000           0
       3 SGMT _n__ TAU  TCB  YYY    0.000000E+00   10.0000           0
       4 SGMT _n__ TAU  TTO  OOO    0.000000E+00   10.0000           0
       5 SGMT _n__ TAU  TTA  XXX    0.000000E+00   10.0000           0
       6 SGMT _n__ TAU  TTB  YYY    0.000000E+00   10.0000           0
       7 SGMT _n__ TAU  TIO  OOO    0.000000E+00   10.0000           0
       8 SGMT _n__ TAU  TIA  XXX    0.000000E+00   10.0000           0
       9 SGMT _n__ TAU  TIB  YYY    0.000000E+00   10.0000           0

       6 !NBOND: bonds
       1       2       1       3       4       5       4       6
       7       8       7       9

       0 !NTHETA: angles


       0 !NPHI: dihedrals


       0 !NIMPHI: impropers


       0 !NDON: donors


       0 !NACC: acceptors


       0 !NNB

       0       0       0       0

       1       0 !NGRP
       0       0       0

end
"""

def massSetup(list=[],axisMass=300):
    """
    Appropriately setup clock masses to axisMass.

    If list is not specified, then pseudoatoms associated with
    all registered PREPot objects are configured. 
    """
    try:
        len(list)
    except TypeError:
        list = [list]
        pass

    from simulation import currentSimulation
    global registeredTerms
    if not list:
        try:
            list = registeredTerms[currentSimulation().lookupID()]
        except:
            pass

    for t in list:
        if t.fixTauc() and t.fixTaut() and t.fixTaui():
            continue
        for a in [t.tcOatom(),t.tcXatom(),t.tcYatom(),
                  t.ttOatom(),t.ttXatom(),t.ttYatom(),
                  t.tiOatom(),t.tiXatom(),t.tiYatom()]:
            a.setMass(axisMass)
            pass
        pass
    return

def topologySetup(ivm,list=[]):
    """
    Configure the given <m ivm>.IVM object's topology setup using the
    freedom string for each PREPot in list.
    This function should be called prior to
    <m ivm>.IVM.autoTorsion() or <m protocol>.torsionTopology().

    """

    try:
        len(list)
    except TypeError:
        list = [list]
        pass

    from simulation import currentSimulation
    global registeredTerms
    if not list:
        try:
            list = registeredTerms[currentSimulation().lookupID()]
        except:
            pass

    for p in list:
        for (fixed,oAtom,xAtom,yAtom) in [
            (p.fixTauc(),p.tcOatom(),p.tcXatom(),p.tcYatom()),
            (p.fixTaut(),p.ttOatom(),p.ttXatom(),p.ttYatom()),
            (p.fixTaui(),p.tiOatom(),p.tiXatom(),p.tiYatom())]:

            if not fixed:
                ivm.fix(xAtom)
                ivm.group( (oAtom,yAtom) )
                ivm.hinge("bend", oAtom,xAtom,yAtom)
            else:
                ivm.fix( (oAtom,xAtom,yAtom) )
                pass
            pass
        pass
    return

def makeTable(pre):
    """Return the assignment table (a string) corresponding to the 
    restraints associated with the specified <m prePot>.PREPot instance. 

    Does not support SUM restraints or non-stereo-assigned restraints.

    """
    ret=""
    for restraint in pre.restraints():
        ret += "assign "
        ret += "\t( %s )\n"  % restraint.sel1.string()
        ret += "\t( %s )  "  % restraint.sel2.string()
        ret += "%7.4f %7.4f\n\n" % (restraint.g(),
                                    restraint.gSigma())
        pass
    return ret

nonstereoAtomMap={}
nonstereoAtomMap['HDX'] = ["HD1","HD2"]
nonstereoAtomMap['HGX'] = ["HG1","HG2"]

def addNonstereoRestraints(pre,inString):
    """
    To the specified <m prePot>.PREPot object, add the restraints in
    inString. These can be standard restraints, or restraints in which
    geminal pairs are labelled with X/Y instead of 1/2. The optimal
    assignment can then be made using the function processNonstereoRestraints.
    """
    import potUtils
    inString = potUtils.stripBracketedComments(inString)
    #
    import re
    nonstereoList=[]
    while True:
        m=re.search(r"assi[a-z]*[^ (]",inString,re.IGNORECASE)
        if not m:
            break
        inString=inString[m.end(0):]
        #m=re.search("assi[a-z]*[^ (]",inString,re.IGNORECASE)
        #if m==None:
        #    break
        indx = m.start(0)
        start=inString.find('(')
        from parseTools import findNested, readFloat
        end = findNested('(',')',start,inString,0)
        sel1 = inString[start+1:end]
        inString=inString[end:]
        start=inString.find('(')
        end = findNested('(',')',start,inString,0)
        sel2 = inString[start+1:end]
        inString=inString[end+1:]
        val,inString=readFloat(inString)
        err,inString=readFloat(inString)
        line=inString.split('\n')[0].strip()
        comment=''
        if line and line[0]=='!':
            comment = line[1:]
            pass
        isNonstereo=False
        for name in list(nonstereoAtomMap.keys()):
            #fix: should also allow other whitespace
            if (re.search(r'name *'+name,sel2,re.IGNORECASE) or
                re.search(r'name *'+name.replace('X','Y'),sel2,re.IGNORECASE)):
                
                nonstereoList.append( (sel1,sel2,val,err,comment) )
                isNonstereo=True
                break
            pass
        if not isNonstereo:
            restraint = 'assign (%s) (%s) %f %f !%s' % (sel1,
                                                        sel2,
                                                        val,
                                                        err,comment)
            pre.addRestraints(restraint)
            pass
        pass
    partners=[]
    from atomSel import AtomSel
    pre.nonStereo=[]
    for sel1,sel2,val,err,comment in nonstereoList:
        for key,values in list(nonstereoAtomMap.items()):
            if sel2.find(key)>-1:
                sel2Partner=sel2.replace(key,key.replace('X','Y'))
                sel2=sel2.replace(key,values[0])
                comment+=" [nonstereo: %s]" % key
                restraint = 'assign (%s) (%s) %f %f !%s' % (sel1,
                                                            sel2,
                                                            val,
                                                            err,comment)
                pre.addRestraints(restraint)
                pre.nonStereo.append( [key,pre.restraints()[-1]] )
                #now look for partner
                foundPartner=False
                for psel1,psel2,pval,perr,pcomment in nonstereoList:
                    print('partner: ',sel2Partner, psel2)
                    if sel2Partner==psel2:
                        pcomment+=" [nonstereo: %s]" % key.replace('X','Y')
                        psel2=sel2.replace(values[0],values[1])
                        restraint = 'assign (%s) (%s) %f %f !%s' % (psel1,
                                                                    psel2,
                                                                    pval,
                                                                    perr,
                                                                    pcomment)
                        pre.addRestraints(restraint)
                        pre.nonStereo[-1].append(pre.restraints()[-1])
                        foundPartner=True
                        break
                    pass
                if not foundPartner:
                    print('Warning: found no partner for nonstereospecific', end=' ')
                    print("restraint", pre.restraints()[-1].name())
                    pass
                break
            pass
        pass
    return 

def processNonstereoRestraints(pre):
    """
    Swap all non-stereospecifically labeled assignments s.t. energy
    is minimized.
    """
    numSwapped=0
    from atomSel import AtomSel
    for partners in pre.nonStereo:
        e0=pre.calcEnergy()
        print(e0)
        if len(partners)==3:
            key,r1,r2 = partners
            r1.sel2, r2.sel2 = r2.sel2, r1.sel2
            e1=pre.calcEnergy()
            if e0<e1:
                #swap back
                r1.sel2, r2.sel2 = r2.sel2, r1.sel2
            else:
                print('swapped sel2 for %s and %s' % (r1.name(),
                                                      r2.name()))
                numSwapped += 1
                pass
            pass
        else: # only a single partner
            key,restraint=partners
            for key,values in list(nonstereoAtomMap.items()):
                if restraint.comment().find("nonstereo: %s" % key)>=0:
                    sel2_0 = restraint.sel2.string()
                    sel2 = sel2_0
                    if sel2.find(values[0])>=0:
                        sel2 = sel2.replace(values[0],values[1])
                    else:
                        sel2 = sel2.replace(values[1],values[0])
                        pass
                
                    restraint.sel2 = AtomSel(sel2,restraint.sel2.simulation())
                    e1=pre.calcEnergy()
                    if e0<e1:
                        #swap back
                        restraint.sel2 = AtomSel(sel2_0,restraint.sel2.simulation())
                    else:
                        numSwapped += 1
                        pass
                    pass
                pass
            pass
        pass
    print("processNonstereoRestraints: %d assignments swapped" % numSwapped)
    return
        
        
        
        
        
        
        

def spaceSeparatedToRestraint(inString,
                              paraSel,
                              atomName="HN",
                              residCol=1,
                              preCol=2,
                              errCol=None,
                              defaultErr=0.1,
                              residCharPrepend=False,
                              illegalValues=[],
                              segid=None,
                              ):
    """
    Convert string restraint table (inString) consisting of columns
    for resid and PRE values to an Xplor-NIH readable restraint
    table. Column numbers start with 1.  paraSel is an atom selection
    specifying the paramagnetic center.

    Set residCharPrepend to True if the residCol takes the form Xnum, where X
    is a chracter representing the residue type, and num is the residue
    number. The illegalValues argument is a sequence of (string) value for
    observed PRE. If these values are encountered, they will be ignored.
    """

    lines = [line for line in inString.split('\n') if not line.startswith('#')]

    ret=""
    first=True
    for line in lines:
        cols = line.split()
        if not cols:
            continue
        try:
            resid = cols[residCol-1]
            if residCharPrepend: resid = resid[1:]
            resid = int(resid)
            pre = cols[preCol-1]
            if pre in illegalValues:
                continue
            pre = float(pre)
            err = float(cols[errCol-1]) if errCol!=None else defaultErr
        except:
            if not first: # allow failure until one read successfully.
                print("Warning: could not read line: %s" % line)
                pass
            continue
            
        atomSel = "name %s and resid %d" % (atomName,resid)
        if segid != None:
            atomSel += ' and segid "%s"' % segid
            pass
        ret += "assign  (%s) (%s) %f %f ! %s\n" % (paraSel,atomSel,
                                                   pre,err,line)
        first=False
        pass
    if not ret:
        print("spaceSeparatedToRestraint: Warning: no restraints read.")
    return ret

def fitTauc(term,verbose=False):
    """Determine the value of tauc which optimizes term's energy
    """
    from minimize import linemin
    term.setFixTauc(True)
    prev_funType = term.funType()
    prev_energyScale = term.scale()
    if prev_funType=="correlation":
        term.setFunType("harmonic")
        term.setScale(1)
        pass
    def f(t):
        term.setTauc(t[0])
        ret = term.calcEnergy()
        if verbose: print('tauc: %.2g, energy: %.2f' % (term.tauc(), ret))
        return ret
    from cdsVector import CDSVector_double as RVec
    bestEnergy=1e20
    bestTauc=None
    for t0 in (0.1, 1, 4, 6, 8, 10, 15, 20):
        ret = linemin(f,RVec([t0]),RVec([1]),5)
        if term.calcEnergy()<bestEnergy:
            bestEnergy = term.calcEnergy()
            bestTauc= ret[0][0]
            if verbose: print('  best tauc: %.2g' % term.tauc())
            pass
        pass
    if verbose: print('final tauc: %.2g, Q-factor: %.3g' % (term.tauc(),
                                                      term.qFactor()))
    tauc=bestTauc
    if tauc<term.tcMin():
        print("fitTauc: Warning. Optimal tauc < tcMin.", end=' ')
        print("Using tcMin (%f)" % term.tcMin())
        tauc=term.tcMin()
        pass
    if tauc>term.tcMax():
        print("fitTauc: Warning. Optimal tauc > tcMax.", end=' ')
        print("Using optimal tauc (%f) and resetting tcMax to this val." % tauc)
        term.setTcMax( tauc )
        pass
    term.setTauc(tauc)
    term.setFunType(prev_funType)
    term.setScale( prev_energyScale )
    return ret[1]

def fitRho0(term):
    """Determine and set constant prefactor rho0 to minimize chi^2 
    """
    term.setRho0(1.)
    from cdsVector import CDSVector_double as RVec
    g=   RVec([r.calcd() for r in term.restraints()])
    gObs=RVec([r.obs() for r in term.restraints()])
    err= RVec([r.err() for r in term.restraints()])
    rho0 = sum(g*gObs/err) / sum(g**2/err)
    term.setRho0( rho0 )
    return 

def analyze(potList):
    """Perform analysis of PREPot terms and return nicely formatted summary."""

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'PRE')

    if not potList: return ret

    #set to get rid of duplicates. list for ordering.
    instanceNames = list(set([x.instanceName() for x in potList]))
    instanceNames.sort()

    ret+= "%9s  %6s  %6s  %6s %7s  %4s\n" % \
          (" " , "RMS", "Q-fac", "corr.", "tau_c", "Viols")

    for name in instanceNames:
        pre = [x for x in potList if x.instanceName()==name][0]

        print(pre.showViolations())

        print(pre.info())

        ret += "%-9s  %6.3f  %6.3f  %6.3f  %7.3f  %4d\n" % \
               (name , pre.rms(), Qfactor(pre), pre.correlation(),
                pre.tc()*1e9, pre.violations())
        pass
    
    return ret

import simulationTools
simulationTools.registerTerm(analyze,
                             "Paramagnetic Relaxation Enhancement Analysis",
                             "PRE",
r"""
For each PREPot term, print the root mean square fit of calculated to
experimental values, the Q-factor, the correlation, the overall correlation
time and the number of violations of this term.

Please see
    J. Iwahara, C.D. Schwieters and G.M. Clore. "Ensemble
    approach for NMR structure refinement against 1H paramagnetic
    relaxation enhancement data arising from a flexible
    paramagnetic group attached to a macromolecule,"
    J. Am. Chem. Soc. 126, 5879-5896 (2004).
""")

simulationTools.registerExtraStats("PRE","Q-factor",Qfactor,True)
simulationTools.registerExtraStats("PRE","correlation",correlation,True)
def tauc(term):
    from simulationTools import potType
    if potType(term)=='PRE':
        return term.tc()*1e9
    else:
        # assume we have a list of terms
        from simulationTools import flattenPotList
        prelist = flattenPotList(term)
        if len(prelist)==0:
            return 0.;
        
        from math import sqrt

        sum=0.
        for p in prelist:
            sum += p.tc()
            pass
        return sum/len(prelist) * 1e9
    return 

simulationTools.registerExtraStats("PRE","tau_c",tauc,True)

registeredTerms={}
def registerTerm(term):
    """
    Add the given PrePot object to a list associated with its Simulation.
    These objects will be automatically processed by topologySetup and
    massSetup.
    """
    global registeredTerms

    #potObj (not PotProxy) used to avoid reference cycle
    # but this should be cleaned up...
    p=term.potObj
    p.thisown=False

    sim = term.simulation()
    id = sim.lookupID()
    if not id in registeredTerms: registeredTerms[id] = []
    registeredTerms[id].append(p)
    while sim.type() in ("SymSimulation", "EnsembleSimulation"):
        if sim.type() == "SymSimulation":
            from symSimulation import fromSimulation
            sim = fromSimulation(sim).subSel().simulation()
            pass
        if sim.type() == "EnsembleSimulation":
            from ensembleSimulation import fromSimulation
            sim = fromSimulation(sim).subSim()
            pass
        id = sim.lookupID()
        if not id in registeredTerms: registeredTerms[id] = []
        registeredTerms[id].append(p)
        pass

    return
    


import pyPot 
import derivList

class twoSatePREPot(pyPot.PyPot):
    """

    """
    def __init__(self, name, popA, prePotA, prePotB):
        pyPot.PyPot.__init__(self, name)
        self.popA = popA
        self.popB = 1.0 - popA
        self.prePotA = prePotA
        self.prePotB = prePotB

        obs = []  # for original observed PRE values
        restraint_pairs = [] 
        for (a, b) in zip(self.prePotA.restraints(), self.prePotB.restraints()):
            obs.append(a.g())  # same as b.g()
            a.setG(999) # to avoid division by zero in calcEnergyAndDerivList
            b.setG(999)
            restraint_pairs.append((a, b))

        self.obs = tuple(obs)
        self.restraint_pairs = tuple(restraint_pairs)
        

    def calcEnergy(self):

        energy = 0

        # Update energy of input terms.
        (self.prePotA.calcEnergy(), self.prePotB.calcEnergy())
        
        for i in range(len(self.obs)):  # loop thru restraints

            restA = self.restraint_pairs[i][0] # restraints
            restB = self.restraint_pairs[i][1]
            
            obs = self.obs[i]
            
            calcA = restA.gamma()
            calcB = restB.gamma()
            calc = self.popA * calcA + self.popB * calcB
            
            weight = restA.weight() # same as restB.weight() 
            
            energy += weight * (obs - calc)**2
            
        return self.scale() * energy 
        
            
    def calcEnergyAndDerivList(self, derivs):

        import xplor
        sim = xplor.simulation

        # Update energy of input terms.
        (self.prePotA.calcEnergy(), self.prePotB.calcEnergy())

        for i in range(len(self.obs)):  # loop thru restraints

            restA = self.restraint_pairs[i][0] # restraints
            restB = self.restraint_pairs[i][1]

            obs = self.obs[i]

            calcA = restA.gamma()
            calcB = restB.gamma()
            calc = self.popA * calcA + self.popB * calcB

            weight = restA.weight()  # restB.weight() is the same
            diffA = restA.g() - calcA  # restA.g() set to 999 by __init__
            diffB = restB.g() - calcB  # restB.g() set to 999 by __init__

            # We are dividing by diffX below; by setting restX.g() to 999, we
            # made sure diffX isn't zero.
            
            # derivX will hold derivatives associated with input term X.
            (derivA, derivB) = (derivList.DerivList(), derivList.DerivList())
            (restA.deriv(1, derivA), restB.deriv(1, derivB))

            # Here, derivX has energy derivatives of input term X.

            derivA[sim] = derivA[sim] / (-2 * weight * diffA)
            derivB[sim] = derivB[sim] / (-2 * weight * diffB)

            # Here, derivX has calc PRE derivatives of input term X.

            derivs[sim] += -2*weight*(obs-calc) * (self.popA * derivA[sim] +
                                                   self.popB * derivB[sim])

        return self.calcEnergy()


    def q_factor(self, sel='all'):
        """Return the Q factor.

        """
        # Update energy of input terms.
        (self.prePotA.calcEnergy(), self.prePotB.calcEnergy())
        
        obs_squared = []
        diff_squared = []

        for i in range(len(self.obs)):  # loop thru restraints

            restA = self.restraint_pairs[i][0] # restraints
            restB = self.restraint_pairs[i][1]

            # Only consider restraint if at least one of its associated atom
            # selections is within input sel.
            if sel != 'all':
                import atomSel
                atoms = atomSel.AtomSel(sel)
                if not atomSel.intersect(atoms, restA.sel1) and \
                   not atomSel.intersect(atoms, restA.sel2): continue
            
            obs = self.obs[i]
            
            calcA = restA.gamma()
            calcB = restB.gamma()
            calc = self.popA * calcA + self.popB * calcB

            obs_squared.append(obs**2)
            diff_squared.append( (obs - calc)**2 )

        return (sum(diff_squared)/sum(obs_squared))**0.5
