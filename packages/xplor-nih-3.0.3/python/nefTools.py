"""
FIX
"""

catPrefixes={ "shifts"   : "nef_chemical_shift_list",
              "spectrum" : "nef_nmr_spectrum",
              "distance" : "nef_distance_restraint_list",
              "dihedral" : "nef_dihedral_restraint_list"
              }

shiftsPrefix="nef_chemical_shift_list_"  # FIX: remove
peaksPrefix ="nef_nmr_spectrum_"         #FIX: remove


def readNEF(filename,
            initPSF=True
            ):
    """
    Read NEF data from the file given by filename into an internal data
    structure. The returned object is a <m cif>.Cif object. If the initPSF
    argument is True, all PSF information is generated from the NEF data.
    """
    dataString=open(filename).read()
    import cif
    cifData = cif.Cif()
    cifData.parse(dataString)

    if not initPSF:
        return cifData
    


    import psfGen
    seqInfo=seqFromNEF(cifData=cifData)
    for (startResid,segid,seq,nterm,cterm,cycle) in seqInfo.seqs:
        seqType="protein"
        singleChar=False
        if len(seq[0])<3:
            seqType='dna'
            singleChar=True
            if len(seq[0])<2:
                seqType='rna'
                pass
            pass
        psfGen.seqToPSF(seq,seqType=seqType,
                        ntermPatch=None if nterm else '',
                        ctermPatch=None if cterm else '',
                        startResid=int(startResid),segName=segid,
                        singleChar=singleChar)
        if cycle:
            if seqType!='protein':
                raise Exception("readNEF: cyclic nucleic acids not yet supported")
            from xplorSimulation import getXplorSimulation
            xplor=getXplorSimulation()
            endResid=startResid+len(seq)-1
            patch='PEPP' if seq[0]=='PRO' else 'PEPT'
            xplor.command('''
             patch %s
               reference=-=(segid "%s" and resid %d)
               reference=+=(segid "%s" and resid %d)
             end
             ''' % (patch,segid,endResid,segid,startResid))
            pass
        pass
    for segid,resid in seqInfo.cisPeptides:
        psfGen.cisPeptide(int(resid),segid)
        pass
    return cifData


def seqFromNEF(pdbRecord=None,
               cifData=None,
               useSeqres=False,
               processBiomt=False,
               failIfInsertion=-1,
               ):
    """Return a list of list of sequences for an input record in mmCIF format-
    True.  """


    if pdbRecord!=None:
        from cif import Cif
        cifData = Cif()
        cifData.setNumDatablocksToRead(1)
        cifData.setModelToRead(1)
        cifData.parse(pdbRecord)
    elif cifData!=None:
        cifData=cifData
    else:
        raise Exception("must specify one of pdbRecord or cif arguments")

    if "nef_molecular_system" not in cifData.keys():
        raise Exception("could not locate the nef_molecular_system saveset")

    block=cifData["nef_molecular_system"]

    #optional entries
    try:
        cis_peptide = block.nef_sequence.cis_peptide
    except:
        cis_peptide = ['.']*len(block.nef_sequence.chain_code)
        pass
    try: 
        linking = block.nef_sequence.linking
    except:
        linking = ['.']*len(block.nef_sequence.chain_code)
        pass
    try: 
        residue_variant = block.nef_sequence.residue_variant
    except:
        residue_variant = ['.']*len(block.nef_sequence.chain_code)
        pass
    
    seqs=[]
    prevSegid=None
    newSeq=True
    cisPeptides=[]
    for (segid,
         resid,
         resName,
         resVariant, #not yet used
         linking,
         cisPeptide) in zip(block.nef_sequence.chain_code,
                            block.nef_sequence.sequence_code,
                            block.nef_sequence.residue_name,
                            residue_variant,
                            linking,
                            cis_peptide):

        if not linking in ('start','middle','end','cyclic','single','dummy'):
            print('seqFromNEF: linking %s not supported.' % linking, end=' ')
            print('Skipping segid=%s, resid=%s, resname=%s' %(segid,resid,resName))
            newSeq=True
            continue

        try:
            resid=int(resid)
        except ValueError:
            print('seqFromNEF: resid %s is not an integer. Skipping')
            continue

        if resVariant and resVariant!='.':
            print("seqFromNEF: residue_variant (%s) not yet supported" % \
                  resVariant)
            pass
        
        if segid==".": segid=""
        if segid:
            segid += ' '*(4-len(segid)) # this for truncated records
            pass
        
        if prevSegid and segid!=prevSegid:
            seqs.append(( beginResid,segid,seq,nterm,False,cycle) )
            newSeq=True
            pass
        
        if newSeq:
            newSeq=False
            cterm=True #normal C-termination
            beginResid=resid
            seq=[]
            #normal N-termination
            nterm=True if (linking=='start' or linking=='single') else False
            cycle=True if linking=='cyclic' else False
            pass

        seq.append(resName)
        if cisPeptide and cisPeptide!='.': cisPeptides.append( (segid,resid) )

        if linking in ['end','single','cycle']:
            if linking=='cycle': cterm=False
            seqs.append(( beginResid,segid,seq,nterm,cterm,cycle) )
            newSeq=True
            pass
        pass
    if not seqs: raise pdbRecordError
    
    class PDBToSeqClass:
        pass
    ret=PDBToSeqClass()
    ret.seqs=seqs
    ret.cisPeptides = cisPeptides
    ret.biomt=None

    return ret
    


def genMetadata():
    import time
    date=time.strftime("%Y-%m-%d_%H:%M:%S")
    program="Xplor-NIH"

    import random
    uuid = program + "_" + date + "_%X" % int(random.uniform(0,16**8))

    import xplor
    metadata="""
save_nef_nmr_meta_data
    _nef_nmr_meta_data.sf_category     nef_nmr_meta_data
    _nef_nmr_meta_data.sf_framecode    nef_nmr_meta_data
    _nef_nmr_meta_data.format_name     nmr_exchange_format
    _nef_nmr_meta_data.format_version  1.1
    _nef_nmr_meta_data.program_name    %s
    _nef_nmr_meta_data.program_version %s
    _nef_nmr_meta_data.creation_date   %s
    _nef_nmr_meta_data.uuid            %s
save_

""" % (program,xplor.version, date, uuid)

    return metadata

def genMolecularSystem(selection="all"):
    """
    Return a string containing a nef_molecular_system saveframe valid for the
    specified selection.
    """
    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection)

    from cif import Cif, CifDatablock, CifCategory
    block = CifDatablock()

    name="nef_molecular_system"
    cat = CifCategory()
    cat["sf_category"] = "nef_molecular_system"
    cat["sf_framecode"] = name
    block["nef_molecular_system"] = cat

    cat = CifCategory()
    for key in (
        "index",
        "chain_code",
        "sequence_code",
        "residue_name",
        "residue_variant",
        "linking",
        "cis_peptide",
        ):
        cat.addKey(key)
        pass
        
    from psfGen import residueTypes
    
    segids={}
    segidList=[]
    curSegid=None
    curResid=None
    for i,atom in enumerate(selection):
        segid=atom.segmentName()
        resid=int(atom.residueNum())
        resname=atom.residueName()
        if segid!=curSegid or resid!=curResid:
            curSegid=segid
            curResid=resid
            if not segid in list(segids.keys()):
                segids[segid]=[]
                segidList.append(segid)
                pass
            segids[segid].append( (resid,resname) )
            pass
        pass
        
    lineNo=0
    for segid in segidList:
        linking="start"
        for i,(resid,resname) in enumerate(segids[segid]):
            if i==len(segids[segid])-1 or resid+1!=segids[segid][i+1][0]:
                linking="end"
                pass
            from atomSel import AtomSel
            if len(AtomSel('segid "%s" and resid %d and pseudo' %(segid,
                                                                  resid))):
                linking='dummy'
                pass
            lineNo += 1
            cat.addValue("index",str(lineNo))
            cat.addValue("chain_code",segid)# if segid else ".")
            cat.addValue("sequence_code",str(resid))
            cat.addValue("residue_name",resname)
            cat.addValue("residue_variant",".")
            cat.addValue("linking",linking)
            cat.addValue("cis_peptide",".")

            linking="middle"
            pass
        pass
    
    block["nef_sequence"] = cat

    cif=Cif()
    block.setIsSaveframe(True)
    cif[name] = block
    cif.setUseTrailingStop(True)
    cif.setUseTrailingSave(True)
    return cif.asString()

def shiftsFromList(list,
                   name):
    """
    Return a formatted NEF record for the chemical shifts in the input list
    as a string with the specified saveframe name.

    The list argument should be a sequence of (val, sel) or (val,err,sel),
    where val is the chemical shift value, err is the chemical shift error,
    and sel is an <m atomSel>.AtomSel specifying which atoms the shift value
    corresponds to.

    The resulting NEF record is sorted by the atoms order in the PSF.
    """
    from cif import Cif, CifDatablock, CifCategory
    block = CifDatablock()

    cat = CifCategory()
    cat["sf_category"] = "nef_chemical_shift_list"
    cat["sf_framecode"] = "nef_chemical_shift_list_"+name

    block["nef_chemical_shift_list"] = cat

    cat = CifCategory()
    for key in (
        "chain_code",
        "sequence_code",    
        "residue_name",     
        "atom_name",        
        "value",            
        "value_uncertainty",
        ):
        cat.addKey(key)
        pass
    

    
    entries=[]
    def addNEFentry(cat,chain_code,
                    sequence_code,
                    residue_name,
                    atom_name,
                    value,
                    var,
                    index):
        """
        """
        entries.append( (
            chain_code,
            sequence_code,
            residue_name,
            atom_name,
            str(value),
            str(var) if var else '.',
            index) )
        return

    from nefTools import toNefAtomname
    degenerate=[]
    from collections import namedtuple
    ShiftEntry = namedtuple('ShiftEntry',['val','err','sel'])

    nlist=[]
    cnt=0
    for entry in list:
        if entry in nlist:
            cnt += 1
            continue
        nlist.append( entry )
        pass
    print('shiftsFromList: removed %d duplicate entries:' % cnt)
    list = nlist

    for entry in list:
        #print len(entry), entry
        if len(entry)==3:
            val,err,selStr = entry
        else:
            val,selStr = entry[:2]
            err=None
            pass

        if type(selStr)==type([]) and len(selStr)==1:
            selStr = selStr[0]
            pass            

        from selectTools import convertToAtomSel
        sel = convertToAtomSel(selStr)
        
        if len(sel)==0:
            print("string %s selects no atoms" % selStr)
            continue

        if len(sel)==1:
            atom=sel[0]
            addNEFentry(cat,atom.segmentName(), str(atom.residueNum()),
                        atom.residueName(),atom.atomName(),
                        val,err,atom.index())
        else:
            degenerate.append( ShiftEntry(val,err,sel) )
            pass
        pass

    #NEF wildcards:"
    # '%' for 'any sequence of digits', equivalent to the regular expression
    # "[0-9]+".

    for i in range(len(degenerate)):
        if i>= len(degenerate):
            break
        foundMatch=False
        wildcard=False
        degi = degenerate[i]
        atom = degi.sel[0]
        for j in range(i+1,len(degenerate)):
            if j>= len(degenerate):
                break
            degj = degenerate[j]
            if degi.sel == degj.sel:
                foundMatch=True
                #print 'degenerate',i,j, degi.val,degj.val
                degenerate = degenerate[:j] + degenerate[j+1:]
                if degi.val==degj.val:
                    foundMatch=True
                    wildcard=True
                    continue
                degenType=None
                if len(degi.sel)==2: #methylene
                    degenType="methylene"
                    pass
                elif len(degi.sel)==6: #2 methyl groups
                    degenType="methyl"
                    pass
                else:
                    print("can not handle entries with selection:", end=' ')
                    print(degi.sel.string())
                    print("  values:", degi.val,degj.val, i,j)
                    continue
                if degi.val>degj.val:
                    suf1='x'; suf2='y'
                else:
                    suf1='y'; suf2='x'
                    pass

                #print 'adding xy'
                if degenType=="methylene":
                    addNEFentry(cat,atom.segmentName(), str(atom.residueNum()),
                                atom.residueName(),atom.atomName()[:-1]+suf1,
                                degi.val,degi.err,atom.index())
                    addNEFentry(cat,atom.segmentName(), str(atom.residueNum()),
                                atom.residueName(),atom.atomName()[:-1]+suf2,
                                degj.val,degj.err,atom.index())
                elif degenType=="methyl":
                    addNEFentry(cat,atom.segmentName(),
                                str(atom.residueNum()),
                                atom.residueName(),
                                atom.atomName()[:-2]+suf1+'%',
                                degi.val,degi.err,atom.index())
                    addNEFentry(cat,atom.segmentName(),
                                str(atom.residueNum()),
                                atom.residueName(),
                                atom.atomName()[:-2]+suf2+'%',
                                degj.val,degj.err,atom.index())
                    pass
                    
                continue
            pass
        if foundMatch:
            continue

        if len( degi.sel ) in [2,3]:
            foundMatch=True
            wildcard=True
            addNEFentry(cat,atom.segmentName(), str(atom.residueNum()),
                        atom.residueName(),atom.atomName()[:-1]+'%',
                        degi.val,degi.err,atom.index())
            pass
        if not foundMatch:
            print('still have sel:',degi.sel.string(), len(degi.sel))
        pass

    #sort entries by PSF order
    entries.sort(key= lambda e1: e1[-1])
    
    for (chain_code,
         sequence_code,
         residue_name,
         atom_name,
         value,
         var,
         index) in entries:
        cat.addValue("chain_code",chain_code)
        cat.addValue("sequence_code",sequence_code)
        cat.addValue("residue_name",residue_name)
        cat.addValue("atom_name",atom_name)
        cat.addValue("value",value)
        cat.addValue("value_uncertainty",var)
        pass

    block["nef_chemical_shift"] = cat

    cif=Cif()
    block.setIsSaveframe(True)
    cif['nef_chemical_shift_list_'+name] = block
    cif.setUseTrailingStop(True)
    cif.setUseTrailingSave(True)
    return cif.asString()
    

    return nefString

def getBlockNames(nef,blockType):
    """
    Return a list of names of the given blockType.
    """
    from nefTools import catPrefixes
    prefix = catPrefixes[blockType]
    
    blocks = [nef[key] for key in nef.keys() if key.startswith(prefix)]

    if len(blocks)==0:
        raise Exception("getBlock: could find no saveset with prefix: " +
                        prefix)

    ret = [ block[prefix].sf_framecode[0][len(prefix)+1:]
            for block in blocks]
    return ret

def getBlock(nef,
             blockType="",
             nefRestraintName=None):
    """
    Given an object nef, returned by <m nefTools>.readNEF, get the block
    of given type named nefRestraintName. If there's a
    single block of this type, its name needn't be specified.
    """

    if not blockType:
        raise Exception("blockType not specified")

    from nefTools import catPrefixes
    prefix = catPrefixes[blockType]

    blocks = [nef[key] for key in nef.keys() if key.startswith(prefix)]

    if len(blocks)==0:
        raise Exception("getBlock: could find no saveset with prefix: " +
                        prefix)

    blockNames = getBlockNames(nef,blockType)

    if len(blocks)==1:
        return blocks[0]
    elif not nefRestraintName:
        mess = "getBlock: found more than one block with prefix: "+ prefix
        mess += "\n  Please choose one of: " + " ".join(blockNames)
        raise Exception(mess)
    else:
        for i,blockName in enumerate(blockNames):
            if nefRestraintName==blockName:
                return blocks[i]
            pass
        raise Exception("getBlock: restraint named %s not found" %
                        nefRestraintName)
        pass
    return None

    
            

def genHeader(datablockName = "test1",
              cif=None):
    """
    Generate boilerplate header fields for NEF entries. Return as a valid
    NEF string.
    """

    if not datablockName.startswith("nef_"):
        datablockName = "nef_" + datablockName
        pass

    nefString = "data_" + datablockName + '\n'

    nefString += genMetadata()

    nefString += genMolecularSystem()

    if cif: cif.parse(nefString)

    return nefString

def makeCif(datablockName = "test1"):

    import cif
    cifData = cif.Cif()
    cifData.setUseTrailingStop(True)
    cifData.setUseTrailingSave(True)
    

    genHeader(datablockName, cif=cifData)
    return cifData
    
def fromNefAtomname(atomname):
    if not set(atomname).isdisjoint( set("@") ):
        raise ValueError("illegal character in atom name")
    atomname = atomname.replace('%','*')
    if atomname[-1]=='x' or atomname[-1]=='y':
        atomname=atomname[:-1]+'#'
        pass
    if len(atomname)>2 and (atomname[-2]=='x' or atomname[-2]=='y'):
        atomname=atomname[:-2]+'#'+atomname[-1]
        pass
    return atomname
    
def toNefAtomname(atomname):
    atomname = atomname.replace('#','%')
    atomname = atomname.replace('*','%')
    atomname = atomname.replace('+','%')
    return atomname
    
def nefComment(oneOrMoreComment):
    """Return a string properly formatted for one or more lines of comment
    in the NEF format (lines should begin with the '# ' string.
    """
    ret = '# ' + oneOrMoreComment.replace("\n","\n# ") + "\n"
    return ret

def distanceRestraintHeader(name,
                            comments=None,
                            origin=None):
    """
    Return a tuple of (commentStr,, <m cif>.Cif) appropriate for a distance
    restraint list with the specified name, where commentStr is a string
    appropriate to place in a NEF file, and origin is a description of the
    source of these distance restraints.
    """
    from nefTools import nefComment
    commentStr=nefComment(comments)

    from cif import Cif, CifDatablock, CifCategory
    cif=Cif()
    block = CifDatablock()
    cat = CifCategory()
    cat["sf_category"] = "nef_distance_restraint_list"
    cat["sf_framecode"] = name
    cat["potential_type"] = "undefined"
    if origin!=None:
        cat["restraint_origin"] = origin
        pass

    block["nef_distance_restraint_list"] = cat
    

    
    cat = CifCategory()
    for key in ("index",
                "restraint_id",
                "chain_code_1",
                "sequence_code_1",
                "residue_name_1",
                "atom_name_1",
                "chain_code_2",
                "sequence_code_2",
                "residue_name_2",
                "atom_name_2",
                "target_value",
                "lower_limit",
                "upper_limit",
                "weight",
                "XplorNIH_label",
                ):
        cat.addKey(key)
        pass

    block["nef_distance_restraint"] = cat
    block.setIsSaveframe(True)
    cif[name] = block
    cif.setUseTrailingStop(True)
    cif.setUseTrailingSave(True)
    return commentStr,cif
    
def addOneDistanceRestraint(cat,
                            selPairs,
                            index,id,
                            target_value,
                            lower_limit,
                            upper_limit,
                            weight=1,
                            XplorNIH_label=""):
    """
    Given a list of pairs of atom selections, add a restraint to cat,
    an input CifCategory.
    Returns a tuple of (id,index,cat)
    """

    
    class FakeAtom:
        def __init__(s,segid,resid,resname,atomname) :
            s.segid=segid
            s.resid=resid
            s.resname=resname
            s.atomname=atomname
            return
        def segmentName(s): return s.segid
        def residueNum(s): return s.resid
        def residueName(s): return s.resname
        def atomName(s): return s.atomname
        pass

    def getWildcard(sel):
        """Given an AtomSel, determine if all atoms can be replaced by a
        single atomname wildcard."""
        if len(sel)<2:
            return None
        segid = sel[0].segmentName()
        resid = sel[0].residueNum()
        resname = sel[0].residueName()
        name = sel[0].atomName()

        #check that all atoms start with same character
        for atom in [atom for atom in sel][1:]:
            if atom.segmentName() != segid or atom.residueNum() != resid:
                return None
            namei = atom.atomName()
            if namei[0]!=name[0]:
                return None
            pass

        wcIndex=1
        wcGood=None
        while wcIndex<len(name):

            # check that all atoms in wildcard selection are in sel, and that
            # there are no leftovers
            wc=name[:-wcIndex] + "*"
            from atomSel import AtomSel
            wcSel=AtomSel('segid "%s" and resid %d and name %s' %
                          (segid,resid,wc))

            #debuggery
            #print sel.string(),wcIndex,len(wcSel)

            for atom in wcSel:
                wcGood=wc
                if not atom in sel:
                    wcGood=False
                    break
                pass
            for atom in sel:
                if not atom in wcSel:
                    wcGood=None
                    break
                pass
            if wcGood: break
            wcIndex += 1
            pass


        if wcGood:
            return FakeAtom(segid,resid,resname,wcGood)
        else:
            return None
    
    seen=[]
    id += 1
    for selPair in selPairs:
        if selPair in seen:
            continue
        seen.append( selPair )

        #FIX: figure out if we can combine the selection into a
        #NEF atomname wildcard
        aWC = getWildcard(selPair[0])
        bWC = getWildcard(selPair[1]) 
        aSel=[aWC] if aWC else selPair[0]
        bSel=[bWC] if bWC else selPair[1]
        for a in aSel:
            for b in bSel:
                index += 1
                cat.addValue("index",str(index))
                cat.addValue("restraint_id",str(id))
                cat.addValue("chain_code_1",a.segmentName())
                cat.addValue("sequence_code_1",str(a.residueNum()))
                cat.addValue("residue_name_1",a.residueName())
                cat.addValue("atom_name_1",toNefAtomname(a.atomName()))
                cat.addValue("chain_code_2",b.segmentName())
                cat.addValue("sequence_code_2",str(b.residueNum()))
                cat.addValue("residue_name_2",b.residueName())
                cat.addValue("atom_name_2",toNefAtomname(b.atomName()))
                cat.addValue("target_value",str(target_value))
                cat.addValue("lower_limit",str(lower_limit))
                cat.addValue("upper_limit",str(upper_limit))
                cat.addValue("weight", str(weight) )
                cat.addValue("XplorNIH_label", XplorNIH_label )
                pass
            pass
        pass
    return (id,index)
