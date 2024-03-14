
def Xplor_readNEF(block,
                  nefRestraintName=None,
                  xplorSim=None,
                  verbose=True):
    """
    Read Torsion angle restraints from the saveframe in the block argument
    into the specified <m xplorSimulation>.XplorSimulation process (the
    current XplorSimulation, by default).

    Alternatively, the block argument could be an object returned by
    <m nefTools>.readNEF(). In this case, the dihedral restraint given by
    nefRestraintName is read, or if that argument is None, and there is a
    single dihedral restraint table in the NEF data block, that is read, while
    if there is more than one such table an exception is raised.

    """

#      _dihedral_restraint.restraint_index
#      _dihedral_restraint.restraint_combination_index
#      _dihedral_restraint.chain_code_1
#      _dihedral_restraint.sequence_code_1
#      _dihedral_restraint.residue_type_1
#      _dihedral_restraint.atom_name_1
#      _dihedral_restraint.chain_code_2
#      _dihedral_restraint.sequence_code_2
#      _dihedral_restraint.residue_type_2
#      _dihedral_restraint.atom_name_2
#      _dihedral_restraint.chain_code_3
#      _dihedral_restraint.sequence_code_3
#      _dihedral_restraint.residue_type_3
#      _dihedral_restraint.atom_name_3
#      _dihedral_restraint.chain_code_4
#      _dihedral_restraint.sequence_code_4
#      _dihedral_restraint.residue_type_4
#      _dihedral_restraint.atom_name_4
#      _dihedral_restraint.weight
#      _dihedral_restraint.target_value
#      _dihedral_restraint.target_value_uncertainty
#      _dihedral_restraint.lower_limit
#      _dihedral_restraint.upper_limit

    import nefTools
    if not hasattr(block,"nef_dihedral_restraint"):
        block = nefTools.getBlock(block,"dihedral",nefRestraintName)
        pass

    restraints=block.nef_dihedral_restraint
    ids = [int(id) for id in restraints.restraint_id]
    ids_set = set(ids)
    for id in ids:
        if not id in ids_set:
            continue
        if ids.count(id)>1:
            print("Xplor_readNEF: restraint id %d has >1 entry. Not supported."%id)
            ids_set.remove(id)
            pass
        pass
    
    from nefTools import fromNefAtomname
    import xplorSimulation
    xplorSim = xplorSimulation.getXplorSimulation(xplorSim)
    
    restraintString=""
    nullChars=('.','?','')
    cnt=0
    from atomSel import AtomSel
    for (id,
         segid1,resid1,name1,
         segid2,resid2,name2,
         segid3,resid3,name3,
         segid4,resid4,name4,
         weight,target,lower,upper) in zip(
        restraints.restraint_id,
        restraints.chain_code_1,restraints.sequence_code_1,restraints.atom_name_1,
        restraints.chain_code_2,restraints.sequence_code_2,restraints.atom_name_2,
        restraints.chain_code_3,restraints.sequence_code_3,restraints.atom_name_3,
        restraints.chain_code_4,restraints.sequence_code_4,restraints.atom_name_4,
        [float(v) for v in restraints.weight],
        [float(v) for v in restraints.target_value],
        [float(v) for v in restraints.lower_limit],
        [float(v) for v in restraints.upper_limit]):

        id=int(id)
        if not id in ids_set:
            continue

        try:
            sel1 = 'resid %d and name %s' % (int(resid1),name1)
            if not segid1 in nullChars: sel1 += ' and segid "%s"' % segid1
            
            sel2 = 'resid %d and name %s' % (int(resid2),name2)
            if not segid2 in nullChars: sel2 += ' and segid "%s"' % segid2
    
            sel3 = 'resid %d and name %s' % (int(resid3),name3)
            if not segid3 in nullChars: sel3 += ' and segid "%s"' % segid3
    
            sel4 = 'resid %d and name %s' % (int(resid4),name4)
            if not segid4 in nullChars: sel4 += ' and segid "%s"' % segid4
    
            val = target
            delta = 0.5 * abs(upper - lower)

            for sel in [AtomSel(str,xplorSim) for str in (sel1,sel2,sel3,sel4)]:
                if len(sel)!=1:
                    mess ="selection (%s) selects %d atoms" % (sel.string(),
                                                               len(sel)) 
                    print(mess)
                    raise Exception(mess)
                pass

            restraint = """assign (%s) (%s)
            (%s) (%s) %f %f %f 2""" % (sel1,sel2,sel3,sel4,
                                       weight,val,delta)
            restraint +=  " ! NEF ID: %d\n" % id
            restraintString += restraint
            cnt += 1
        except:
            print("Xplor_readNEF: error reading restraint",id)
            pass
        pass
    if verbose:
        print("readNEF: adding %d restraints" % cnt)
        pass
    import protocol
    protocol.initDihedrals(string=restraintString,simulation=xplorSim,
                           useDefaults=False)
    return

def readSegResName(line):
    """return tuple of (segid, resid, resname, atomname) parsed from the
    input string"""
    fields=line.split()
    if len(fields)==4:
        segid=fields[0]
        resid=int(fields[1])
        resname=fields[2]
        name=fields[3]
    elif len(fields)==3:
        segid=""
        resid=int(fields[0])
        resname=fields[1]
        name=fields[2]
    else:
        raise Exception("unable to parse atom record: line")
    return (segid,resid,resname,name)

def Xplor_writeNEF(name,
                   xplorSim=None,
                   origin=None):
    """
    Return a formatted NEF record for the torsion angle restraints in the
    specified <m xplorSimulation>.XplorSimulation as a string with the
    specified saveframe name. The optional original argument can be used to
    specify the source of the restraints.
    """

    from xplorSimulation import currentSimulation, getXplorSimulation
    if not xplorSim:
        xplorSim = getXplorSimulation(currentSimulation())
        pass

    from cif import Cif, CifDatablock, CifCategory
    block = CifDatablock()

    cat = CifCategory()
    cat["sf_category"] = "nef_dihedral_restraint_list"
    cat["sf_framecode"] = name
    cat["potential_type"] = "square-well-parabolic"
    if origin!=None:
        cat["origin"] = origin
        pass
    

    block["nef_dihedral_restraint_list"] = cat

    cat = CifCategory()
    for key in (
        "index",
        "restraint_id",
#        "restraint_combination_id",
        "chain_code_1",
        "sequence_code_1",
        "residue_name_1",
        "atom_name_1",
        "chain_code_2",
        "sequence_code_2",
        "residue_name_2",
        "atom_name_2",
        "chain_code_3",
        "sequence_code_3",
        "residue_name_3",
        "atom_name_3",
        "chain_code_4",
        "sequence_code_4",
        "residue_name_4",
        "atom_name_4",
        "weight",
        "target_value",
        "lower_limit",
        "upper_limit",):
        cat.addKey(key)
        pass
    

    outputState=xplorSim.disableOutput()

    from simulationTools import mktemp
    cmd=r'"print thres=-1 CDIH%s"'
    xplorFilename=mktemp('xplor-violations')
    xplorSim.fastCommand("set print %s end" % xplorFilename)
    xplorSim.fastCommand(eval(cmd))
    xplorSim.fastCommand("set print $prev_print_file end" )
    xplorSim.fastCommand("close %s end" % xplorFilename)
    xplorSim.enableOutput(outputState)
    
    lines=[l[:-1] for l in open(xplorFilename).readlines()]
    import os
    from nefTools import toNefAtomname
    os.unlink(xplorFilename)
    cnt=0
            

    id=0
    index=0
    while cnt < len(lines):
        id += 1
        if '==============' in lines[cnt]:
            atom1 = readSegResName(lines[cnt+1])
            atom2 = readSegResName(lines[cnt+2])
            atom3 = readSegResName(lines[cnt+3])
            atom4 = readSegResName(lines[cnt+4])
            weight = float(lines[cnt+5].split()[5])
            target = float(lines[cnt+5].split()[7])
            delta = float(lines[cnt+6].split()[1])
            cnt += 6
            
            index += 1
            cat.addValue("index",str(index))
            cat.addValue("restraint_id",str(id))
            cat.addValue("chain_code_1",atom1[0])
            cat.addValue("sequence_code_1",str(atom1[1]))
            cat.addValue("residue_name_1",atom1[2])
            cat.addValue("atom_name_1",toNefAtomname(atom1[3]))
            cat.addValue("chain_code_2",atom2[0])
            cat.addValue("sequence_code_2",str(atom2[1]))
            cat.addValue("residue_name_2",atom2[2])
            cat.addValue("atom_name_2",toNefAtomname(atom2[3]))
            cat.addValue("chain_code_3",atom3[0])
            cat.addValue("sequence_code_3",str(atom3[1]))
            cat.addValue("residue_name_3",atom3[2])
            cat.addValue("atom_name_3",toNefAtomname(atom3[3]))
            cat.addValue("chain_code_4",atom4[0])
            cat.addValue("sequence_code_4",str(atom4[1]))
            cat.addValue("residue_name_4",atom4[2])
            cat.addValue("atom_name_4",toNefAtomname(atom4[3]))
            cat.addValue("weight",str(weight))
            cat.addValue("target_value",str(target))
            cat.addValue("lower_limit",str(target-delta))
            cat.addValue("upper_limit",str(target+delta))
            
            pass
        cnt += 1
        pass
    block["nef_dihedral_restraint"] = cat

    cif=Cif()
    block.setIsSaveframe(True)
    cif[name] = block
    cif.setUseTrailingStop(True)
    cif.setUseTrailingSave(True)
    return cif.asString()

    
def XPLOR_makeTable(xplorSim=None):
    """
    Return a XPLOR-formatted restraint table as a string for the 
    specified <m xplorSimulation>.XplorSimulation.
    """

    from xplorSimulation import currentSimulation, getXplorSimulation
    if not xplorSim:
        xplorSim = getXplorSimulation(currentSimulation())
        pass

    outputState=xplorSim.disableOutput()

    from simulationTools import mktemp
    cmd=r'"print thres=-1 CDIH%s"'
    xplorFilename=mktemp('xplor-violations')
    xplorSim.fastCommand("set print %s end" % xplorFilename)
    xplorSim.fastCommand(eval(cmd))
    xplorSim.fastCommand("set print $prev_print_file end" )
    xplorSim.fastCommand("close %s end" % xplorFilename)
    xplorSim.enableOutput(outputState)
    
    lines=[l[:-1] for l in open(xplorFilename).readlines()]
            

    ret="! restraint table created by Xplor_makeTable\n"
    cnt=0
    while cnt < len(lines):
        if '==============' in lines[cnt]:
            atoms=[]
            for i in range(4):
                atoms.append( readSegResName(lines[cnt+i+1]) )
                pass
            weight = lines[cnt+5].split()[5]
            target = lines[cnt+5].split()[7]
            delta  = lines[cnt+6].split()[1]
            cnt += 6
            sels=[]
            for i in range(4):
                sels.append( "resid %s and name %s" % (atoms[i][1],
                                                       atoms[i][3]))
                if atoms[i][0]!='':
                    sels[-1] += ' and segid "%s"' % atoms[i][0]
                    pass
                pass

            ret += "assign (%s) (%s) (%s) (%s) %s %s %s 2\n" \
                   % (tuple(sels)+(weight,target,delta))
            pass
        cnt += 1
        pass
    return ret

    
