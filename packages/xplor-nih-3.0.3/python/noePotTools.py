
"""Tools to aid in setup/analysis of the NOE potential term.

This module provides functions to simplify the creation, manipulation and
analysis of <m noePot>.NOEPot potential terms.
"""

def create_NOEPot(name,file='',esim=0,restraints='',
                  nef=None,nefRestraintName=None,
                  splitRestraints=False,deltaResidLR=5,
                  simulation2=None,
                  verbose=False):
    """Return an <m noePot>.NOEPot instance.

    The name of the potential is given by the argument name (a string).
    Restraints are specified via an NOE assignment table (whose filename is
    indicated in the file argument; a string), and/or via a string of such
    assignments (restraints argument).  The file argument can optionally be a
    sequence of assignment table filenames.  esim specifies an ensemble
    simulation (optional). The simulation2 argument is set when the second
    selection in a restraint should correspond to a different
    <m simulation>.Simulation.

    If specified, the nef argument should be an object returned by
    <m nefTools>.readNEF. If there is more than one distance restraint table
    specified in the NEF record, the desired name must be specified in the
    nefRestraintName argument.

    If splitRestraints=True, the returned potential is a <m potList>.PotList
    with <m noePot>.NOEPot instances corresponding to the following classes:

      intraresidue
      interresidue - sequential
      interresidue - short range (1 < delta resid <= deltaResidLR)
      interresidue - long range (delta resid > deltaResidLR)
      ambiguous

    The classes are arranged in the above order within the <m potList>.PotList.    

    """
    from noePot import NOEPot
    from simulation import currentSimulation

    if not esim: esim = currentSimulation()

    import types
    from builtins import type
    if type(file) in [tuple, list]:
        files=file
    elif file:
        files=[file]
    else:
        files=[]
        pass

    noe = NOEPot(name,restraints,esim,simulation2)
    if restraints and noe.numRestraints()==0:
        print("Warning: no restraints read from restraints argument")
        pass

    for file in files:
        startNum=noe.numRestraints()
        noe.addRestraints(open(file).read())
        if noe.numRestraints() == startNum:
            print("Warning: no restraints read from file", file)
            pass
        pass

    noe.setVerbose(verbose)

    #setup default values
    noe.setAveType("sum")
    noe.setPotType( "hard" )
    noe.setRSwitch( 0.5 )
    noe.setAsympSlope( 1. )
    noe.setSoftExp(1.)
    noe.setThreshold( 0.5 )

    if nef:
        import nefTools
        block = nefTools.getBlock(nef,"distance",nefRestraintName)
        readNEF(noe,block,verbose)
        pass
    elif nefRestraintName!=None:
        raise Exception("nefRestraintName specified without nef object")

    if splitRestraints:
        categories=splitTable(noe,deltaResidLR)
        from potList import PotList
        noe = PotList(name)
        types=('int','seq','sr','lr','amb')
        for i in range(len(types)):
            type=types[i]
            restraints=categories[type]
            if restraints:
                noe.append( create_NOEPot("%s-%s"%(name,type),esim=esim,
                                          restraints=restraints) )
                pass
            pass

        list(map( lambda p: p.setAveType("sum"), noe))
        list(map( lambda p: p.setPotType( "hard" ), noe))
        list(map( lambda p: p.setRSwitch( 0.5 ), noe))
        list(map( lambda p: p.setAsympSlope( 1. ), noe))
        list(map( lambda p: p.setSoftExp(1.), noe))
        list(map( lambda p: p.setThreshold( 0.5 ), noe))
        pass

    return noe


def repeated_restraints(noe, selstrings=True):
    """Find NOE restraints with similar selection pairs.

    Return a list of the type:

    [( (seli1, seli2), [(selj1, selj2), ...] ), ...]

    where (selk1, selk2) is a pair of selection strings (if selstrings=True,
    otherwise <m atomSel>.AtomSel instances) associated with a restraint of the
    input noe potential (a <m noePot>.NOEPot instance).  In the above example,
    (seli1, seli2) has similar selection pairs [(selj1, selj2), ...] within the
    noe restraint table.  Two individual selections are considered similar if
    they select at least one common atom.

    """

    def coll(s):
        """Collapse white space in string s. """
        return ' '.join(s.split())
    
    from atomSel import intersect
    
    data = {}
    data_str = {}  # string version of data
    restraints = noe.restraints()
    
    for r in restraints:

        if len(r.selPairs()) > 1:
            msg = 'Restraints with more than one selection pair are unsupported'
            raise Exception(msg) 
        
        if not data:
            data[(r.sel1, r.sel2)] = []  # first (arbitrary) seed
            data_str[(coll(r.sel1.string()), coll(r.sel2.string()))] = []
        else:
            # if (r.sel1, r.sel2) or (r.sel2, r.sel1) in data.keys()
            for (sel1, sel2) in list(data.keys()):
                if ((intersect(sel1, r.sel1) and intersect(sel2, r.sel2)) or
                    (intersect(sel1, r.sel2) and intersect(sel2, r.sel1))):
                    data[(sel1, sel2)].append((r.sel1, r.sel2))
                    data_str[(coll(sel1.string()),
                        coll(sel2.string()))].append((coll(r.sel1.string()),
                                                      coll(r.sel2.string())))
                    break
            else: 
                data[(r.sel1, r.sel2)] = []  # new (arbitrary) seed
                data_str[(coll(r.sel1.string()), coll(r.sel2.string()))] = []

    # We only want restraints with "repetitions".
    if selstrings is True:
        return [x for x in list(data_str.items()) if x[1]]
    else:
        return [x for x in list(data.items()) if x[1]]



def splitTable(noe,deltaResidLR):
    ret={}
    ret['int']=""
    ret['amb']=""
    ret['seq']=""
    ret['sr']=""
    ret['lr']=""

    for r in noe.restraints():

        distClassPrev=None
        first=True
        entry=""
        for selPair in r.selPairs():
            if first:
                entry+="assign (%s) (%s) %f %f %f !%s\n" %(selPair.a.string(),
                                                           selPair.b.string(),
                                                           r.d(),
                                                           r.dMinus(),
                                                           r.dPlus(),
                                                           r.comment())
            else:
                entry+=" OR (%s) (%s)\n" % (selPair.a.string(),
                                            selPair.b.string())
                pass
            
            deltaResid=calcDeltaResid(selPair,deltaResidLR)

            if deltaResid == -1:
                distClass='amb'
            elif deltaResid == 0:
                distClass='int'
            elif deltaResid == 1:
                distClass='seq'
            elif deltaResid <= deltaResidLR:
                distClass='sr'
            else:
                distClass='lr'
                pass
            if distClassPrev and distClassPrev!=distClass:
                print("splitTable: Warning: Restraint has multiple classifications"+\
                      " so it will be denoted ambiguous.")
                distClass='amb'
                pass
            distClassPrev=distClass
            first=False
            pass

        ret[distClass] += entry
        pass
    return ret

def calcDeltaResid(selPair,deltaResidLR):
    """
    given a selPair from an NOE restraint, determine the difference in the
    residue number of the two atom selections. If the selections have
    differing segids, deltaResid is set to deltaResidLR.

    If a selection specifies atoms with differing residue IDs, -1 is returned.
    """

    (segid1,resid1) = getSegidResid(selPair.a)
    (segid2,resid2) = getSegidResid(selPair.b)

    if resid1<0 or resid2<0:
        return -1

    if segid1!=segid2:
        return deltaResidLR

    return abs(resid1-resid2)

def getSegidResid(sel):
    """
    given an <m atomSel>.AtomSel object, return the tuple
    (segid,resid), if they are unique. If they are not unique return
    ('',-1)
    """

    (segid, resid) = (None, None)
    for atom in sel:
        if not segid: segid=atom.segmentName()
        if not resid: resid=atom.residueNum()

        if segid!=atom.segmentName() or resid!=atom.residueNum():
            return ('',-1)
        pass
    return (segid,resid)
    
def correlation(term):
    """Term can be a single potential term, or a list of terms, in which case
    the overall correlation is returned. FIX: should return overall
    correlation instead.
    """
    from simulationTools import potType
    if potType(term)=='NOEPot':
        from cdsVector import CDSVector_double, sum

        if term.numRestraints()==0:
            return 0.

        if term.numRestraints()<3:
            return 1.

        calcd = CDSVector_double([r.dist() for r in term.restraints()])
        obs   = CDSVector_double([r.d() for r in term.restraints()])
        
        ave_calcd = sum(calcd) / len(calcd)
        ave_obs   = sum(obs) / len(obs)
        num  = 0.;
        sumCalcd2=0.;
        sumObs2 = 0.;

        from math import sqrt
        for calcd,obs in zip(calcd,obs):
            diffCalcd = calcd - ave_calcd
            diffObs   = obs   - ave_obs
            num += diffCalcd * diffObs;
            sumCalcd2 += diffCalcd**2;
            sumObs2   += diffObs**2;
            pass
        denom = sqrt(sumCalcd2 * sumObs2);
        try:
            ret = num / denom;
        except ZeroDivisionError:
            return 0.
        
        return ret
    else:
        # assume we have a list of terms
        from simulationTools import flattenPotList
        prelist = flattenPotList(term)

        if len(prelist)==0:
            return 0.

        sum=0
        for p in prelist:
            sum += correlation(p)
            pass

        return sum / len(prelist)
    return 


def analyze(potList):
    "perform analysis of NOEPot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'NOEPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret += "%-9s  %6s  %6s  %6s  %6s\n" % \
           ( "", "RMS", "Devia", "Viols", "Corr.")

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        print(term.showViolations())
        print(term.info())

        ret += "%-9s  %6.3f  %6.3f  %6d  %6.4f\n" % \
               (name , term.rms(), term.deviation(), term.violations(),
                correlation(term))
        pass
    
    return ret

def makeTable(noe):
    """given an noePot, generate the associated assignment 
    table. Return this as a string.

    If the argument is a list-like object, makeTable will be called
    recursively on each element.
    """
    ret=""
    try:
        len(noe)
        # if we get here, noe is a list-like object
        for term in noe:
            ret += "\n# Restraints Generated from term: " + term.instanceName()
            ret += "\n\n"
            ret += makeTable(term)
            pass
        pass
    except TypeError:
        
        for r in noe.restraints():
            ret +="assign (%s) (%s) %7.3f %7.3f %7.3f !%s\n" %(r.sel1.string(),
                                                      r.sel2.string(),
                                                      r.d(),
                                                      r.dMinus(),
                                                      r.dPlus(),
                                                      r.comment())
            for selPair in r.selPairs()[1:]:
                ret += " OR "
                ret += "( %s ) "  % selPair.a.string()
                ret += "( %s ) \n"  % selPair.b.string()
                pass
            pass
        pass
    return ret

def removeLowerBounds(noe):
    """Given an NOEPot object, remove all lower bounds from the restraints. This
    can be appropriate during early stages of structure calculation to improve
    convergence."""
    for r in noe.restraints():
        r.setDMinus( r.d() )
        pass
    return

def genRestraints(selection,
                  selection2=None,
                  distCutoff=6.0,
                  distOffset=0.0,
                  upperBound=1.0,
                  lowerBound=1.0,
                  omitIntraResidue=True,
                  deltaSequence=0,
                  asList=False
                  ):
    """
    Given current coordinates, generate distance restraints from all
    atoms in the argument selection to all atoms in the argument
    selection2 which are within the distance specified by
    distCutoff. The nominal distance is taken from the atom selections,
    while the +/- spread (both values) are set to upperBound and
    lowerBound, respectively. The value distOffset is added to the nominal
    observed distance. If selection2 is omitted, it defaults to
    selection. If distCutoff is less than or equal to zero, there is
    no distance cutoff: restraints are generated between all atoms
    specified. 

    Intra-residue restraints are omitted if omitTranResidue is True, while
    restraints between residues whose distances in primary sequence is less
    than deltaSequence are also not generated.

    The return value is a string which contains the restraints, unless asList
    is True, in which case the return value is a list of strings, each element
    of which contains a single restraint.
    """

    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)
    selection2 = convertToAtomSel(selection2) if selection2!=None else selection

    # used to avoid doubling restraints
    seen={} if selection.simulation()==selection2.simulation() else None
    
    ret=[] if asList else ''
    from vec3 import norm
    for atom1 in selection:
        for atom2 in selection2:
            if atom1==atom2:
                continue
            dist = norm(atom1.pos()-atom2.pos())
            if distCutoff>0. and dist>distCutoff:
                continue

            index1=atom1.index()
            index2=atom2.index()
            if seen!=None and index2 in seen and index1 in seen[index2]:
                continue
            if index1 not in seen: seen[index1]=set()
            seen[index1].add(index2)

            resid1 = atom1.residueNum()
            resid2 = atom2.residueNum()
            if (omitIntraResidue and
                atom1.segmentName()==atom2.segmentName() and
                resid1 == resid2                            ):
                continue                

            sel1 = "resid %d and name %s" % (resid1,atom1.atomName())
            sel2 = "resid %d and name %s" % (resid2,atom2.atomName())
            segid1 = atom1.segmentName()
            segid2 = atom2.segmentName()

            if segid1==segid2 and abs(resid1-resid2)<deltaSequence:
                continue

            if segid1!="": sel1 += ' and segid "%s"' % segid1
            if segid2!="": sel2 += ' and segid "%s"' % segid2
            
            restraint = "ASSIgn (" + sel1 + ")\n"
            restraint += "       (" + sel2 + ") "
            dist += distOffset
            restraint += "%.3f %.3f %.3f" % (dist,lowerBound,upperBound)

            if asList:
                ret.append(restraint)
            else:
                ret += restraint + '\n';
            pass
        pass
    return ret

def readNEF(noe,block,verbose=False):
    """
    Given an NOEPot, and a saveframe from a NEF record, read in all
    specified restraints.
    """
    #   _distance_restraint_list.sf_category     distance_restraint_list
    #   _distance_restraint_list.sf_framecode    distance_restraint_list_1
    #   _distance_restraint_list.potential_type  square-well-parabolic
    #
    #   loop_
    #      _distance_restraint.restraint_index
    #      _distance_restraint.restraint_combination_id
    #      _distance_restraint.chain_code_1
    #      _distance_restraint.sequence_code_1
    #      _distance_restraint.residue_name_1
    #      _distance_restraint.atom_name_1
    #      _distance_restraint.chain_code_2
    #      _distance_restraint.sequence_code_2
    #      _distance_restraint.residue_name_2
    #      _distance_restraint.atom_name_2
    #      _distance_restraint.weight
    #      _distance_restraint.target_value
    #      _distance_restraint.target_value_uncertainty
    #      _distance_restraint.lower_limit
    #      _distance_restraint.upper_limit

    from nefTools import fromNefAtomname
    try:
        pottype = block.nef_distance_restraint_list.potential_type[0]
    except:
        pottype=None
        pass
    
    cnt=0
    restraints = block.nef_distance_restraint
    ids = [int(id) for id in restraints.restraint_id]
    idIndices={}
    for i,id in enumerate(ids):
        if id in idIndices:
            idIndices[id].append(i)
        else:
            idIndices[id] = [i]
            pass
        pass
    
    from atomSel import AtomSel
    if "XplorNIH_label" in restraints.keys():
        haveComments = True
        comments = restraints.XplorNIH_label
    else:
        haveComments = False
        pass
    for id in sorted(idIndices.keys()):
        try:
            first=True
            fromSels=[]
            toSels=[]
            for i in idIndices[id]:
                fromSegid=restraints.chain_code_1[i]
                fromResid=int(restraints.sequence_code_1[i])
                fromSel = 'resid %d' % fromResid
                if fromSegid!='.':
                    fromSel += ' and segid "%s"' % fromSegid
                    pass
                resName = AtomSel(fromSel)[0].residueName()
                fromName=fromNefAtomname(restraints.atom_name_1[i])
                fromSel += " and name %s" % fromName

                toSegid=restraints.chain_code_2[i]
                toResid=int(restraints.sequence_code_2[i])
                toSel = 'resid %d' % toResid
                if toSegid!='.': #FIX: check this
                    toSel += ' and segid "%s"' % toSegid
                    pass
                resName = AtomSel(toSel)[0].residueName()
                toName=fromNefAtomname(restraints.atom_name_2[i])
                toSel += " and name %s" % toName

                if first:
                    try:
                        weight = float(restraints.weight[i])
                    except (ValueError, IndexError):
                        weight=1
                        pass
                    if pottype=="upper-bound-parabolic":
                        dist=float(restraints.upper_limit[i])
                        upper=0
                        lower=dist
                    else:
                    
                        dist  = float(restraints.target_value[i])
                        lower = dist - float(restraints.lower_limit[i])
                        upper = float(restraints.upper_limit[i]) - dist
                        pass
                        
                    first=False
                    pass
                fromSels.append(fromSel)
                toSels.append(toSel)
                pass

            restraint = "assign (%s) (%s) %f %f %f" % (fromSels[0],
                                                           toSels[0],
                                                           dist,lower,upper)
            
            
            comment = comments[i] if haveComments else "NEF ID: %d" % id
            restraint +=  " ! " + comment + '\n'
            for fromSel,toSel in zip(fromSels[1:],toSels[1:]):
                restraint += "OR (%s) (%s)\n" % (fromSel,toSel)
                pass

            try:
                noe.addRestraints(restraint)
                noe.restraints()[-1].setWeight( float(weight) )
                pass
            except SystemError as e:
                print(e)
                print("error reading restraint id %s: %s" % (id, restraint))
                pass
            cnt+=1
        except ValueError:
            print("Error adding restraint with ID ",str(id))
            pass
        except IndexError:
            print("Error adding restraint with ID ",str(id))
            pass
        pass
    if verbose:
        print("readNEF: added %d of %d restraints" % (noe.numRestraints(),cnt))
    return

def writeNEF(noe,name):
    """Given an NOEPot term, write a NEF record with the specified saveframe
    name."""

    from nefTools import distanceRestraintHeader
    commentStr,cif = distanceRestraintHeader(name,
                                  "\nRestraints from NOEPot term: %s\n"%
                                  noe.instanceName())
    
    cat = cif[name]["nef_distance_restraint"]

    from nefTools import toNefAtomname
    id=0
    index=0
    for r in noe.restraints():
        d=r.d()
        lower = d - r.dMinus()
        upper = d + r.dPlus()

        r.selPairs(),
        selPairs = [(p.a,p.b) for p in r.selPairs()]

        from nefTools import addOneDistanceRestraint
        (id,index) = addOneDistanceRestraint(cat,selPairs,
                                             index,id,
                                             d,lower,upper,
                                             r.weight(),
                                             #r.comment()
                                             )
        pass
    return commentStr + cif.asString()

def spaceSeparatedToRestraint(lines,
                              atom1Name="HN",
                              atom2Name="HN",
                              resid1Col=1,
                              resid2Col=2,
                              distCol=2,
                              errCol=None,
                              defaultErr=1,
                              segid1=None,
                              segid2=None,
                              ):
    """
    Convert input columnar data in condensed form to an Xplor-NIH
    readable restraint table. The first argument is either a multiline
    string or sequence of strings. Lines beginning with a hash
    character (#) are skipped. A single atom name is specified for
    each selection using atom1Name and atom2Name. The residue number
    is specified by columns specified by resid1Col and
    resid2Col. Distance and +/- errors are specified by distCol and
    errCol, respectively. If errCol=None, the values of defaultErr is
    used instead. The arguments segid1 and segid2 specify segmentNames
    for the two selections. Column numbers start with 1.

    To assign distances to a single nucleus resid2Col can be set to
    -1, and then atom2Name is taken as a complete atom selection,
    ignoring resid2Col and segid2.
    """

    if type(lines)==type(str("")):
        lines = [line for line in lines.split('\n') 
                 if not line.startswith('#')]
        pass

    ret=""
    for line in lines:
        cols = line.split()
        if not cols:
            continue
        try:
            resid1 = int(cols[resid1Col-1])
            resid2 = int(cols[resid2Col-1]) if resid2Col>=0 else resid2Col
            dist= float(cols[distCol-1])
            err = float(cols[errCol-1]) if errCol!=None else defaultErr
        except IndexError:
            print("Warning: could not read line: %s" % line)
            continue
            
        atom1SelString = "name %s and resid %d" % (atom1Name,resid1)
        if segid1 != None:
            atom1SelString += ' and segid "%s"' % segid1
            pass
        if resid2<0:
            atom2SelString = atom2Name
        else:
            atom2SelString = "name %s and resid %d" % (atom2Name,resid2)
            if segid2 != None:
                atom2SelString += ' and segid "%s"' % segid2
                pass
            pass
        ret += "assign (%s) (%s)  %f %f %f ! %s\n" % (atom1SelString,
                                                      atom2SelString,
                                                      dist,err,err,line)
        pass
    return ret

nonstereoAtomMap={}
nonstereoAtomMap['HDX'] = ["HD1","HD2"]
nonstereoAtomMap['HGX'] = ["HG1","HG2"]

def addNonstereoRestraints(term,inString):
    """
    To the specified <m noePot>.NOEPot object, add the restraints in
    inString. These can be standard restraints, or restraints in which
    geminal pairs are labelled with X/Y instead of 1/2. The optimal
    assignment can then be made using the function
    processNonstereoRestraints.

    This function does not yet support restraints using the ASSIgn
    ... OR idiom. Nor does it support restraints in which both
    selections contain unassigned geminal pairs.
    
    """
    import potUtils
    inString = potUtils.stripBracketedComments(inString)
    #
    import re
    nonstereoList=[]
    while True:
        m=re.search(r"^ *assi[a-z]*[^ (]",inString,re.IGNORECASE|re.MULTILINE)
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
        errMinus,inString=readFloat(inString)
        errPlus,inString=readFloat(inString)
        line=inString.split('\n')[0].strip()
        comment=''
        if line and line[0]=='!':
            comment = line[1:]
            pass
        isNonstereo=False
        for name in list(nonstereoAtomMap.keys()):
            #fix: should also allow other whitespace
            if (re.search(r'name *'+name,sel2,re.IGNORECASE) or
                re.search(r'name *'+name.replace('X','Y'),sel2,
                          re.IGNORECASE)                        or
                re.search(r'name *'+name,sel1,re.IGNORECASE) or
                re.search(r'name *'+name.replace('X','Y'),sel1,
                          re.IGNORECASE)):
                nonstereoList.append((sel1,sel2,val,errMinus,errPlus,comment))
                isNonstereo=True
                break
            pass
        if not isNonstereo:
            restraint = 'assign (%s) (%s) %f %f %f !%s' % (sel1,
                                                           sel2,
                                                           val,
                                                           errMinus,
                                                           errPlus,
                                                           comment)
            term.addRestraints(restraint)
            pass
        pass
    partners=[]
    from atomSel import AtomSel
    term.nonStereo=[]
    for sel1,sel2,val,errMinus,errPlus,comment in nonstereoList:
        for key,values in list(nonstereoAtomMap.items()):
            if sel2.find(key)>-1:
                sel2Partner=sel2.replace(key,key.replace('X','Y'))
                sel2=sel2.replace(key,values[0])
                comment+=" [nonstereo: %s]" % key
                restraint = 'assign (%s) (%s) %f %f %f !%s' % (sel1,
                                                            sel2,
                                                            val,
                                                            errMinus,
                                                            errPlus,comment)
                term.addRestraints(restraint)
                term.nonStereo.append( [key,term.restraints()[-1]] )
                #now look for partner
                foundPartner=False
                for (psel1,psel2,pval,
                     perrMinus,perrPlus,pcomment) in nonstereoList:
                    print('partner: ',sel2Partner, psel2)
                    if sel1==psel1 and sel2Partner==psel2:
                        pcomment+=" [nonstereo: %s]" % key.replace('X','Y')
                        psel2=sel2.replace(values[0],values[1])
                        restraint = 'assign (%s) (%s) %f %f %f !%s' % (
                            psel1,
                            psel2,
                            pval,
                            perrMinus,
                            perrPlus,
                            pcomment)
                        term.addRestraints(restraint)
                        term.nonStereo[-1].append(term.restraints()[-1])
                        foundPartner=True
                        break
                    pass
                if not foundPartner:
                    print('Warning: found no partner for nonstereospecific', end=' ')
                    print("restraint", term.restraints()[-1].name())
                    pass
                break
            if sel1.find(key)>-1:
                sel1Partner=sel1.replace(key,key.replace('X','Y'))
                sel1=sel1.replace(key,values[0])
                comment+=" [nonstereo: %s]" % key
                restraint = 'assign (%s) (%s) %f %f %f !%s' % (sel1,
                                                            sel2,
                                                            val,
                                                            errMinus,
                                                            errPlus,comment)
                term.addRestraints(restraint)
                term.nonStereo.append( [key,term.restraints()[-1]] )
                #now look for partner
                foundPartner=False
                for (psel1,psel2,pval,
                     perrMinus,perrPlus,pcomment) in nonstereoList:
                    print('partner: ',sel1Partner, psel1)
                    if sel2==psel2 and sel1Partner==psel1:
                        pcomment+=" [nonstereo: %s]" % key.replace('X','Y')
                        psel1=sel1.replace(values[0],values[1])
                        restraint = 'assign (%s) (%s) %f %f %f !%s' % (
                            psel1,
                            psel2,
                            pval,
                            perrMinus,
                            perrPlus,
                            pcomment)
                        term.addRestraints(restraint)
                        term.nonStereo[-1].append(term.restraints()[-1])
                        foundPartner=True
                        break
                    pass
                if not foundPartner:
                    print('Warning: found no partner for nonstereospecific', end=' ')
                    print("restraint", term.restraints()[-1].name())
                    pass
                break
            pass
        pass
    return 

def processNonstereoRestraints(term):
    """
    Swap all non-stereospecifically labeled assignments s.t. energy
    is minimized.
    """
    numSwapped=0
    from atomSel import AtomSel
    for partners in term.nonStereo:
        e0=term.calcEnergy()
        print(e0)
        if len(partners)==3:
            key,r1,r2 = partners
            r1.selPairs_, r2.selPairs_ = r2.selPairs(), r1.selPairs()
 

            e1=term.calcEnergy()
            if e0<e1:
                #swap back
                r1.selPairs_, r2.selPairs_ = r2.selPairs(), r1.selPairs()
            else:
                print('swapped %s and %s' % (r1.name(),r2.name()))
                numSwapped += 1
                pass
            pass
        else: # only a single partner
#            #FIX: need to be able to swap atom selections for this
#            raise Exception("single partner not supported")
            key,restraint=partners
            from noePot import CDSList_NOESelPair
            for key,values in list(nonstereoAtomMap.items()):
                if restraint.comment().find("nonstereo: %s" % key)>=0:
                    selPairs = restraint.selPairs_
                    selPairs0 = CDSList_NOESelPair(selPairs)
                    for i in range(len(selPairs)):
                        sel2_0 = selPairs[i].b.string()
                        sel2 = sel2_0
                        if sel2.find(values[0])>=0:
                            sel2 = sel2.replace(values[0],values[1])
                        else:
                            sel2 = sel2.replace(values[1],values[0])
                            pass
                        sel1_0 = selPairs[i].a.string()
                        sel1 = sel1_0
                        if sel1.find(values[0])>=0:
                            sel1 = sel1.replace(values[0],values[1])
                        else:
                            sel1 = sel1.replace(values[1],values[0])
                            pass
                        restraint.selPairs_[i].a = AtomSel(sel1,restraint.
                                                           sel1.simulation())
                        restraint.selPairs_[i].b = AtomSel(sel2,restraint.
                                                           sel2.simulation())
                        pass
                    e1=term.calcEnergy()
                    if e0<e1:
                        #swap back
                        restraint.selPairs_ = selPairs0
                    else:
                        numSwapped += 1
                        pass
                    pass
                pass
            pass
        pass
    term.calcEnergy()
    print("processNonstereoRestraints: %d assignments swapped" % numSwapped)
    return
        

from simulationTools import registerTerm, registerExtraStats
registerTerm(analyze,"NOE terms","NOE",
r"""
For each term, the following are reported::

 RMS    - root mean square deivation between calculated and observed
 Devia  - deviation between ensemble members
 Viols  - number of violated terms
 Corr.  - correlation between calculated and observed distances
""")

registerExtraStats("NOEPot","Corr.",correlation)
