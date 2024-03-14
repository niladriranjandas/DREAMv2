
"""Tools for converting from variance chemical shift formats to that used by
the Xplor-NIH potential terms.
"""

def convertRestraints(restraints,
                      format,
                      saveSet=None,
                      segid=None,
                      ambiguousGlyHA=False,
                      verbose=False):
    """Convert restraints to Xplor-NIH native format.
    The format of this string can be 'plain' (internal Xplor-NIH format),
    'TALOS', 'PIPP', or 'NMRSTAR'. The saveSet argument can be specified to
    choose a specific set of chemical shifts from a NMRSTAR table.

    The segid argument can be used to add a segment name to each restraint.

    If ambiguousGlyHA is True, HA1 and HA2 stereo specific chemical
    shifts will be made ambiguous.
    """
    ret=""

    if format.lower()=="plain":
        ret = restraints
    elif format.lower()=="talos":
        ret = convertFromTalos(restraints,
                               ambiguousGlyHA=ambiguousGlyHA,
                               segid=segid,
                               verbose=verbose)
    elif format.lower()=="nmrstar":
        from pasd import readSTAR, starShifts
        from simulationTools import mktemp
        tmpfilename = mktemp()
        open(tmpfilename,'w').write(restraints)
        readSTAR(tmpfilename)
        import os
        os.unlink(tmpfilename)
        shiftData = starShifts(segmentName=segid,
                               saveSet=saveSet,
                               tclOutput=False)
        ret = ""
        #FIX: ambiguousGlyHA?
        for val,sels in shiftData:
            fullSel = " or ".join(["(%s)" % sel for sel in sels])
            bad=False
            from atomSel import AtomSel
            for atom in AtomSel(fullSel):
                if (verbose and
                    not atom.atomName() in "C CA CB HA HN N".split()):
                    print("Warning. Sparta does not support atom", end=' ')
                    print(atom.atomName())
                    bad=True
                    break
                pass
            if not bad:
                ret += "assi (%s) %f\n" % (fullSel,float(val))
                pass
            pass
        pass
    elif format.lower()=="pipp":
        ret = convertFromPipp(restraints,
                              ambiguousGlyHA=ambiguousGlyHA,
                              verbose=verbose)
        pass
    elif format.lower()!='plain':
        raise Exception("unsupported format: %s" % format)
    
    return ret

def convertFromTalos(restraints,
                     segid=None,
                     ambiguousGlyHA=False,
                     verbose=False):
    """
    Parse the given string restraint table and split it into tables for each
    supported atom type present. If segid is not specified, then the
    selection omits a segid selection.
    """


    start=False
    tables={}
    for line in restraints.split('\n'):
        line=line.strip()
        if not line or line[0]=='#' or line.startswith('REMARK'):
            continue
        if line.startswith('FORMAT'):
            start=True
            continue
        if not start:
            continue
        (resid,resname,name,val) = line.split()

        from atomSel import AtomSel
        if segid!=None:
            selstr='segid "%s" and ' % segid
        else:
            selstr=''
        selstr +='resid %d'%int(resid)
            
        tagSel=AtomSel(selstr)
        if len(tagSel)==0:
            if verbose:
                print("convertFromTalos: WARNING:", end=' ')
                print("tagSel contains no atoms: " +tagSel.string())
                pass
            continue
        if len(resname)==1:
            import selectTools
            resname=selectTools.renameResidues([resname])[0]
            pass
        if verbose and resname!=tagSel[0].residueName():
            print("convertFromTalos: WARNING: residue name mismatch:")
            print("   %4s != %4s" % (resname,tagSel[0].residueName()), end=' ')
            if segid:
                print(" segid: %4s" % segid, end=' ')
                pass
            print(" resid: %d" %int(resid))
            pass
                                     
        (nameType,name)=nameConvertType(name)
        if not nameType:
            continue
        if not nameType in list(tables.keys()):
            tables[nameType]=''
            pass
        if resname=="GLY" and name.startswith("HA") and ambiguousGlyHA:
            name="HA#"
            pass
        if name:
            tables[nameType]+= 'assi (%s and name %s) %s\n' % (selstr,
                                                               name,val)
            pass
        pass
    return tables

def convertFromPipp(restraints,
                    segid=None,
                    ambiguousGlyHA=False,
                    verbose=False):
    """
    Parse the given string restraint table and split it into tables for each
    supported atom type present. If segid is not specified, then the
    selection omits a segid selection.
    """

    if ambiguousGlyHA: raise Exception("not yet supported")

    inDef=False
    tables={}
    for line in restraints.split('\n'):
        line=line.strip()
        if not line or line[0]=='#':
            continue
        if line.startswith('RES_ID'):
            (text,resid) = line.split()
            continue
        if line.startswith('RES_TYPE'):
            (text,resname) = line.split()
            inDef=True
            continue
        if line.startswith('SPIN_SYSTEM_ID'):
            continue
        if line.startswith('END_RES_DEF'):
            inDef=False
            continue

        if not inDef:
            continue

        if segid!=None:
            selstr='segid "%s" and ' % segid
        else:
            selstr=''
        selstr +='resid %d'%int(resid)

        from atomSel import AtomSel
        tagSel=AtomSel(selstr)
        if len(tagSel)==0:
            if verbose:
                print("convertFromPipp: WARNING:", end=' ')
                print("tagSel contains no atoms: " +tagSel.string())
                pass
            continue
            pass

        if verbose and resname!=tagSel[0].residueName():
            print("convertFromPipp: WARNING: residue name mismatch:")
            print("   %4s != %4s" % (resname,tagSel[0].residueName()), end=' ')
            if segid:
                print(" segid: %4s" %segid, end=' ')
            print(" resid: %d" %int(resid))
            pass

        (name,val) = line.split()
        
        if "|" in name:
            names=name.split('|')
            name=names[0]
            for n in names[1:]:
                name += ' or name ' + n
                pass
            pass

#        if "#" in name: # ignore these entries for now
#            continue

        if not name:
            continue

        atomType,name = nameConvertType(name)

        if not name:
            continue

        if not atomType in list(tables.keys()):
            tables[atomType]=''
            pass
        if name:
            tables[atomType]+= '''assi (%s and (name %s))
                                   %s\n''' % (selstr,
                                              name,
                                              val)
            pass
        pass
    return tables

def nameConvertType(name): #,swapHA23=False):
    """
    Given an atom name, return a tuple consisting of the corresponding
    Sparta atom type and a possibly converted atom name
    """
    ret=(None,None)
    if name in ("N","HA","C","CA","CB","HN"):
        return (name,name)

    if name.startswith("HA"):
        ret = ("HA",name)
        pass

    if  name=="HA3":
        ret = ("HA","HA1")
        pass
    return ret
    
