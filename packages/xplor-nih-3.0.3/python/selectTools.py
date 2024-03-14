"""tools for manipulating atom selections
"""
#
# rigidProteinSelections:
#  list of rigid regions of protein sidechains
#
# """resname PRO and 
#    (name N or name CA or name CB or name CG or name CD)""",
rigidProteinSelections = {}
for (res,sel) in \
[("PHE","name CG or name CD1 or name CD2 or name CE1 or name CE2 or name CZ"),
("HIS","name CG or name ND1 or name CD2 or name CE1 or name NE2"),
("TYR","name CG or name CD1 or name CD2 or name CE1 or name CE2 or name CZ"),
("TRP","""name CG or name CD1 or name CD2 or name NE1 or name CE2 or
       name CE3 or name CZ2 or name CZ3 or name CH2"""),
("ASN","name CG or name OD1 or name ND2"),
("ASP","name CG or name OD*"),
("GLN","name CD or name OE1 or name NE2"),
("GLU","name CD or name OE*"),
("ARG","name NE or name CZ or name NH*"),
("CED","name MN or name OY* or name CY* or name NY*"),
("CYSP","name " +
        " or name ".join("CS6 CS7 CS5 NS1 OS1 CS2 CS3 CS4 CS8 CS9".split()))
]:
    rigidProteinSelections[res] = sel


#
# rigidNucleicSelections:
#  list of rigid regions of nucleic acid bases
#
rigidNucleicSelections = {}
for (res, sel) in \
[("CYT", """name N1 or name C6 or name C5 or name C4 or name N3 or name C2 or
        name N4 or name H5 or name H6 or name O2 or name H4#"""), 
 
 ("GUA", """name N9 or name C4 or name N3 or name C2 or name N1 or name N2 or
       name C6 or name C5 or name N7 or name C8 or name O6
       or name H1 or name H8 or name H2#"""), 
 
 ("ADE", """name N9 or name C4 or name N3 or name C2 or name N1 or
       name C6 or name C5 or name N7 or name C8 or name N6
       or name H2 or name H8 or name H6#"""), 

 ("AP7", """name N9 or name C4 or name N3 or name C2 or name N1 or
       name C6 or name C5 or name N7 or name C8 or name N6
       or name H2 or name H8 or name H1 or name H6#"""), 
 
 ("TED", """name N1 or name C6 or name C5 or name C4 or name N3 or name C2
        or name H6 or name C7 or name O4 or name H3 or name O2"""),
        # add H7, H8, C8, C9 ?
 
 ("THY", """name N1 or name C6 or name C5 or name C4 or name C5A
       or name N3 or name C2 or name O4 or name O2 or name H3 or name H6"""),
 
 ("URI", """name N1 or name C6 or name C5 or name C4 
       or name N3 or name C2 or name O4 or name O2 or name H3 or name H5
       or name H6""")]:
    
    rigidNucleicSelections[res] = sel



def groupRigidSideChains(sel="not PSEUDO",
                         sim=None):
    """return a list of groups atoms which are members of fixed
    regions of protein sidechains.
    The second argument is a selection in sim in which the groupings are
    to be made."""
    return groupRigidSel(rigidProteinSelections,sel,sim)

def groupRigidProtein(sel="not PSEUDO",
                      sim=None):
    """return a list of groups atoms which are members of fixed
    regions of protein backbone and sidechains.
    The second argument is a selection in sim in which the groupings are
    to be made."""
    import simulation
    if not sim:     sim = simulation.currentSimulation()
    if len(sel)==0: sel = sim.select("name ca")
    resIDs = []
    for i in sel:
        id = sim.residueNum(i)
        if not id in resIDs:
            resIDs.append(id)
        pass
    groupList = groupRigidSideChains(sel,sim)
    for resID in resIDs:
        prev = resID - 1
        sel = sim.select("""(resid %d and (name N or name HN)) or
        (resid %d and (name C or name O))""" % (resID,prev))
        if len(sel):
            groupList.append( sel )
            pass
        pass
    return groupList


def groupRigidNucleicAcid(sel="not PSEUDO",
                          sim=None):
    """return a list of groups atoms which are members of fixed
    regions of nucleic acid bases.
    The second argument is a selection in sim in which the groupings are
    to be made."""
    import simulation
    if not sim: sim = simulation.currentSimulation()
    return groupRigidSel(rigidNucleicSelections,sel,sim)


def getSegsResidues(sel="tag",
                    sim=0):
    """given an atom selection return a dictionary whose keys are segment
    names. The value of each dictionary entry is a list of named tuples of
    (resid, resname) for each residue contained therein.
    """
    sel = convertToAtomSel(sel,sim)
    sim = sel.simulation()
    

    segments={}
    from collections import namedtuple
    ResidResname = namedtuple('ResidResname',['resid','resname'])
    for i in sel.indices():
        seg = sim.segmentName(i)
        if seg not in segments: segments[seg]=[]
        res = sim.residueNum(i)
        resName = sim.residueName(i)
        if not (res,resName) in segments[seg]:
            segments[seg].append( ResidResname(res,resName) )
            pass
        pass
    return segments


def groupRigidSel(rigidSelList,
                  sel="not PSEUDO",
                  sim=None):
    """Return a list of groups atoms. The groups are formed by iterating
    through all residues in the selection, and making groups from each
    group selection string.
    The last argument is a selection in sim in which the groupings are
    to be made."""
    import simulation
    if not sim:
        sim = simulation.currentSimulation()
        pass
    sel = convertToAtomSel(sel,sim)
    
    segments = getSegsResidues(sel,sim)

    groupList=[]
    for segID in list(segments.keys()):
        for (resID, resName) in segments[segID]:
            if not resName in list(rigidSelList.keys()):
                continue
            selStr = rigidSelList[resName]
            sel = sim.select('segid "%s" and resid %d and (%s)' %
                             (segID,resID,selStr))
            if len(sel):
                groupList.append( sel )
                pass
            pass
        pass
    return groupList

def getResids(sel=0):
    """return a list of unique residue ids. If an atom selection sel is
    specified, return the unique residue ids in the selection.
    """
    from atomSel import AtomSel
    if type(sel)==type("string"): sel = AtomSel(sel)
    
    if not sel: sel = AtomSel("tag")
    ret = []
    for atom in sel:
        if not atom.residueNum() in ret: 
            ret.append( atom.residueNum() )
            pass
        pass
    return ret

def getSegids(sel="tag"):
    """return a list of unique segment ids. If an atom selection sel is
    specified, return the unique segment ids in the selection.
    """
    sel=convertToAtomSel(sel)

    ret = []
    for atom in sel:
        if not atom.segmentName() in ret: 
            ret.append( atom.segmentName() )
            pass
        pass
    return ret

def maxResid(sel=0,sim=0):
    """return the largest residue id
    """
    from atomSel import AtomSel
    if type(sel)==type("string"): sel = AtomSel(sel)
    
    import simulation
    if not sel: sel = AtomSel("tag")
    if not sim: sim = simulation.currentSimulation()
    ret = 0
    for atom in sel:
        if atom.residueNum()>ret: ret = atom.residueNum()
        pass
    return ret

def minResid(sel=0,sim=0):
    """return the smallest residue id
    """
    from atomSel import AtomSel
    if type(sel)==type("string"): sel = AtomSel(sel)
    
    import simulation
    if not sel: sel = AtomSel("tag")
    if not sim: sim = simulation.currentSimulation()
    ret = 100000
    for atom in sel:
        if atom.residueNum()<ret: ret = atom.residueNum()
        pass
    return ret

def numResidues(sel="not PSEUDO"):
    """return the number of residues in the specified <m atomSel>.AtomSel
    """
    from atomSel import AtomSel
    if type(sel)==type("string"): sel = AtomSel(sel)

    ret = 0
    segResidList = []
    for atom in sel:
        segResid = atom.segmentName() + "~%d" % atom.residueNum()
        if not segResid in segResidList:
            ret += 1
            segResidList.append(segResid)
            pass
        pass
    return ret

cisTol=7 #tolerance in degrees
def isCisPeptide(resid,segid=""):
    """
    Return True if the resid/segid specifies a protein residue
    and the omega angle between the specified residue and its predecessor is
    in the range -cisTol ... +cisTol.
    """
    from dihedral import Dihedral
    try:
        d = Dihedral('name CA and resid %d and segid "%s"' % (resid-1,segid),
                     'name C  and resid %d and segid "%s"' % (resid-1,segid),
                     'name N  and resid %d and segid "%s"' % (resid,segid),
                     'name CA and resid %d and segid "%s"' % (resid,segid)   )
        from math import pi
        val = d.value() * 180./pi
        if abs(val)<cisTol:
            return True
        pass
    except:
        pass
    return False


def IVM_groupRigidSidechain(ivm,
                            sel="not PSEUDO"):
    """Group rigid sidechain regions in selection sel for the <m ivm>.IVM object

    Default protein and nucleic acid rigid regions are grouped, as specified in
    rigidProteinSelections and rigidNucleicSelections.
    """
    sim = ivm.simulation
    sel = convertToAtomSel(sel,sim)

    #don't consider sidechains in rigid groups
    from atomSel import union
    rigidRegions = union( *ivm.groups() )
    from atomSel import intersection, notSelection
    sel =  intersection(sel, notSelection(rigidRegions))
    
    #combine nucleic and protein rigid group info
    d={}
    for key in list(rigidProteinSelections.keys()):
        d[key] = rigidProteinSelections[key]
        pass
    for key in list(rigidNucleicSelections.keys()):
        d[key] = rigidNucleicSelections[key]
        pass
    gl = ivm.groupList()

    for group in groupRigidSel(d,sel,sim):
        gl.append( group )
        pass
    ivm.setGroupList(gl)
    return

def IVM_groupRigidBackbone(ivm,
                           sel="name CA"):
    """For the <m ivm>.IVM object, group omega angles in selection sel. Due
    to implementation details, the dihedral angle which is constrained is
    currently actually the HN-N-C-O angle, so that the bond angle term is
    essential to maintain peptide bond planarity. 
    """
    from atomSel import AtomSel
    sim = ivm.simulation
        
    sel = convertToAtomSel(sel,sim)

    segments = getSegsResidues(sel,sim)

    from atomSel import AtomSel
    for segID in list(segments.keys()):
        for (resID, resName) in segments[segID]:

            prev = resID - 1
            aSel = AtomSel('''
            (segid "%s" and resid %d and name N HN) or
            (segid "%s" and resid %d and name C O)''' % (segID,resID,
                                                          segID,prev))
            if len(aSel): ivm.group( aSel )
            pass
        pass
    return

def IVM_breakProlines(ivm,
                      sel=0,
                      breakSelStr='name CD or name CG'):
    """In an <m ivm>.IVM object, break proline rings at a specified position.

    The position of the break is specified by the argument breakSelStr, a string
    with an XPLOR-style selection of two covalently bonded ring atoms.

    """
    # WARNING: default value of breakSelStr should match that of breakProlineSel
    # in protocol.torsionTopology.
    
    sim = ivm.simulation

    if not sel:
        sel = sim.select("tag")
        pass
    sel = convertToSelectionIndices(sel,sim)

    segments = getSegsResidues(sel,sim)
    
    from atomSel import AtomSel
    for segID in list(segments.keys()):
        for (resID, resName) in segments[segID]:
            if resName!='PRO':
                continue
            aSel = AtomSel('''segid "%s" and resid %d and
                           (%s)''' % (segID, resID, breakSelStr))
            if len(aSel)==2:
                ivm.breakBond(aSel)
                pass
            pass
        pass
    return
    
def IVM_breakRiboses(ivm,
                     sel="tag",
                     breakSelStr="name C4' or name O4'"):
    """In an <m ivm>.IVM object, break ribose rings at a specified position.

    The position of the break is specified by the argument breakSelStr, a string
    with an XPLOR-style selection of two covalently bonded ring atoms.    

    """
    # WARNING: default value of breakSelStr should match that of breakRiboseSel
    # in protocol.torsionTopology.
    
    sim = ivm.simulation
    from selectTools import convertToAtomSel
    sel = convertToAtomSel(sel,sim)

    #don't consider riboses in rigid groups
    from atomSel import union
    rigidRegions = union( *ivm.groups() )
    from atomSel import intersection, notSelection
    sel =  intersection(sel, notSelection(rigidRegions))

    segments = getSegsResidues(sel,sim)
    
    from atomSel import AtomSel
    for segID in list(segments.keys()):
        for (resID, resName) in segments[segID]:
            # should check here that the residue is a nucleic acid
            aSel = AtomSel('segid "%s" and resid %d and '  % (segID,resID) +
                           "(%s)"% breakSelStr)
            if len(aSel)==2:
                ivm.breakBond(aSel)
                ivm.constrainBond(aSel)
                pass
            pass
        pass
    return

def IVM_breakDisulfides(ivm,
                        sel=0):
    """for the <m ivm>.IVM object, break disulfide bonds in given
    selection (default: all disulfides)"""

    sim = ivm.simulation
    if not sel:
        sel = sim.select("tag")
        pass
    sel = convertToSelectionIndices(sel,sim)

    segments = getSegsResidues(sel,sim)
    

    from atomSel import AtomSel
    seen=[]
    for segID in list(segments.keys()):
        for (resID, resName) in segments[segID]:
            if resName!='CYS':
                continue
            segRes = "%s_%d" %(segID,resID)
            if segRes in seen:
                continue
            seen.append(segRes)
            aSel = AtomSel('segid "%s" and resid %d and name SG' %
                           (segID,resID))
            bSel=0
            if len(aSel)>1:
                raise "IVM_breakDisulfides: CYS with more than one sulfur??"
            elif aSel:
                bSel = AtomSel("name SG and bondedto id %d" %
                               (aSel[0].index()+1) )
                pass

            if bSel:
                seen.append( "%s_%d" % (bSel[0].segmentName(),
                                        bSel[0].residueNum()) )
                ivm.breakBond(AtomSel("id %d or id %d" %
                                      (aSel[0].index()+1, bSel[0].index()+1) ))
                
                pass
            pass
        pass
    return

def IVM_breakIsopeptides(ivm,
                         sel=0):
    """for the <m ivm>.IVM object, break isopeptide bonds in given
    selection (default: all isopeptides)"""

    sim = ivm.simulation
    if not sel:
        sel = sim.select("tag")
        pass
    sel = convertToSelectionIndices(sel,sim)

    segments = getSegsResidues(sel,sim)
    

    from atomSel import AtomSel
    seen=[]
    for segID in list(segments.keys()):
        for (resID, resName) in segments[segID]:
            if resName!='LYS':
                continue
            segRes = "%s_%d" %(segID,resID)
            if segRes in seen:
                continue
            seen.append(segRes)
            aSel = AtomSel('segid "%s" and resid %d and name NZ' %
                           (segID,resID))
            bSel=0
            if len(aSel)>1:
                raise "IVM_breakIsopeptides: LYS with more than one NZ??"
            elif aSel:
                bSel = AtomSel("name C and bondedto id %d" %
                               (aSel[0].index()+1) )
                pass

            if bSel:
                seen.append( "%s_%d" % (bSel[0].segmentName(),
                                        bSel[0].residueNum()) )
                ivm.breakBond(AtomSel("id %d or id %d" %
                                      (aSel[0].index()+1, bSel[0].index()+1) ))
                
                pass
            pass
        pass
    return

def IVM_flexibilizeRiboses(ivm, sel="all"):
    """Make all endocyclic bond angles within ribose rings flexible.

    A subset of rings can be made flexible by specifying the sel argument
    to be an <m atomSel>.AtomSel.
    
    """
    sim = ivm.simulation
    sel = convertToSelectionIndices(sel, sim)
    segments = getSegsResidues(sel, sim)
    
    from atomSel import AtomSel
    for segid in list(segments.keys()):
        for (resid, resname) in segments[segid]:
            if not len( AtomSel('segid "%s" and resid %d and name C3\'' %
                                (segid, resid)) ):
                continue
                            
            ivm.hinge('bendtorsion',
                      'segid "%s" and resid %d and name C3\'' % (segid, resid),
                      'segid "%s" and resid %d and name C2\'' % (segid, resid),
                      'segid "%s" and resid %d and name C4\'' % (segid, resid))
            ivm.hinge('bendtorsion',
                      'segid "%s" and resid %d and name C2\'' % (segid, resid),
                      'segid "%s" and resid %d and name C1\'' % (segid, resid),
                      'segid "%s" and resid %d and name C3\'' % (segid, resid))
            ivm.hinge('bendtorsion',
                      'segid "%s" and resid %d and name C1\'' % (segid, resid),
                      'segid "%s" and resid %d and name O4\'' % (segid, resid),
                      'segid "%s" and resid %d and name C2\'' % (segid, resid)) 
    return

    
residueMapProtein={     'G' : 'GUA' ,   #single character residue names
                        'A' : 'ALA',
                        'R' : 'ARG',
                        'N' : 'ASN',
                        'D' : 'ASP',
                        'C' : 'CYS',
                        'c' : 'cys',  #Bax group convention
                        'Q' : 'GLN',
                        'E' : 'GLU',
                        'G' : 'GLY',
                        'H' : 'HIS',
                        '#' : 'HIH',  #Bax group convention
                        'h' : 'his',  #Bax group convention
                        'I' : 'ILE',
                        'L' : 'LEU',
                        'K' : 'LYS',
                        'M' : 'MET',
                        'F' : 'PHE',
                        'P' : 'PRO',
                        'S' : 'SER',
                        'T' : 'THR',
                        'W' : 'TRP',
                        'Y' : 'TYR',
                        'V' : 'VAL',
                        }

def renameResidues(seq,type='protein'):
    """ Convert single letter to three letter residue names. The type
    argument may be 'protein' (default) or 'nucleic'.
    """
    if type=='protein':
        residueMap = residueMapProtein
        pass
    elif type=='nucleic':
        import psfGen
        residueMap = psfGen.residueMap

    ret = []
    for r in seq:
        if r.isspace():
            continue
        r=r.strip()
        if r.upper() in list(residueMap.keys()):
            ret.append( residueMap[r] )
        else:
            ret.append(r)
            pass
        pass
    return ret

def threeToOne(name):
    """
    Convert a three character residue name to the corresponding one
    character version.
    """
    name  = name.upper()
    for ( one,three ) in list(residueMapProtein.items()):
        if three==name:
            return one
        pass
    import psfGen
    residueMapNucleic = psfGen.residueMap
    for ( one,three ) in list(residueMapNucleic.items()):
        if three==name:
            return one
        pass
    return None


def oneToThree(name):
    """
    Convert a one character residue name to the corresponding three
    character version.
    """
    if name==" ": return name
    name  = name.upper()
    for ( one,three ) in list(residueMapProtein.items()):
        if one==name:
            return three
        pass
    import psfGen
    residueMapNucleic = psfGen.residueMap
    for ( one,three ) in list(residueMapNucleic.items()):
        if one==name:
            return three
        pass
    return None


def convertToAtom(arg,sim=0):
    """
    given an <m atomSel>.AtomSel or a string convert to atom selection
    and check that just a single <m atom>.Atom is selected. Return that Atom.
    If the argument is an Atom, just return it, if its simulation matches the
    sim argument, otherwise, try to find an atom with the same segid, resid,
    and name in simulation sim (or currentSimulation if sim is not specified).
    """

    import simulation
    from atomSel import AtomSel
    if not sim:
        sim = simulation.currentSimulation()
        pass

    if (repr(arg).startswith(r"<atom.Atom") or
        repr(arg).endswith(r"p_Atom>")       ):
        if arg.simulation()==sim:
            return arg
        else:
            return convertToAtom('atom "%s" %d %s' % (arg.segmentName(),
                                                      arg.residueNum(),
                                                      arg.atomName()),
                                 sim)
        pass
    
    if type(arg)==type(""):
        arg = AtomSel(arg,sim)
        pass
        
    if len(arg)!=1 :
        raise Exception("selection does not select a single atom: %s" %
                        arg.string())

    return arg[0]

    

def convertToAtomSel(arg,
                     sim=None,
                     ordered=None):
    """
    convert arg to an <m atomSel>.AtomSel, if it is not already one. The
    argument can be an AtomSel, a selection string, a sequence of atoms, a
    sequence of atom indices, or a raw SWIG pointer string.

    If sim is specified it is used instead of <m simulation>.currentSimulation.
    If sim is specified, arg is an AtomSel, and the simulations don't match,
    a new AtomSel is created using the string from arg and the sim argument.

    If specified, the ordered argument will specify the value of ordered()
    for the returned AtomSel.
    """

    import simulation
    from atomSel import AtomSel

    if (repr(arg).startswith("<atomSel.AtomSel") or
        repr(arg).endswith(r"p_AtomSel>")          ):
        if sim and arg.simulation().lookupID() != sim.lookupID():
            if ordered==None: ordered=arg.ordered()
            arg = AtomSel(arg.string(), sim, ordered)
            pass
        pass
    elif repr(arg).endswith("_p_AtomSel'"):
        #convert from raw pointer
        from atomSel import AtomSelPtr
        from atomSel import AtomSel_fromPtr
        #want arg to be owned by Python interface: copy it
        tmp = AtomSel_fromPtr(arg)
        arg = AtomSel(tmp.string(),tmp.simulation(),tmp.ordered())
        if sim and arg.simulation().lookupID() != sim.lookupID():
            if ordered==None: ordered=arg.ordered()
            arg = AtomSel(arg.string(), sim, ordered)
            pass
        pass
    elif type(arg)==type("string"):
        if not sim: sim = simulation.currentSimulation()
        if ordered==None: ordered=False
        arg = AtomSel(arg,sim,ordered=ordered)
        pass
    elif hasattr(arg,"__len__"):
        if len(arg)>0 and type(arg[0])!=type(1):
            if not sim: sim = arg[0].simulation()
            arg = AtomSel([atom.index() for atom in arg],sim)
        else:
            if not sim: sim = simulation.currentSimulation()
            arg = AtomSel(arg,sim)
            pass
        pass
    else:
        raise TypeError("could not convert >%s< to an AtomSel" % str(arg))

    return arg

    

def convertToSelectionIndices(sel,
                              sim=0):
    """allowed input: string, <m atom>.Atom, <m atomSel>.AtomSel, single index,
        or list of indices
    output: list of indices"""
    import types
    import re
    import simulation
    from atomSel import AtomSel
    if not sim:
        sim = simulation.currentSimulation()

    if type(sel) == tuple: sel = list(sel)
    if repr(sel).endswith(r"p_CDSListTAtom_t>"): sel = list(sel)
    if type(sel) == int:
        return [ sel ]
    if type(sel) == list:
        for i in range( len(sel) ):
            if not type(sel[i])==int:
                try:
                    #assume the list element is an atom
                    sel[i] = sel[i].index()
                except:
                    raise Exception("bad element in array (index %d): %s" %
                                    (i, repr(sel[i])))
                pass
            pass
        return sel
        pass

    if type(sel) == bytes or type(sel) == type("string"):
        return AtomSel(sel,sim).indices()
    

    if (repr(sel).startswith(r"<atom.Atom") or
        repr(sel).endswith(r"p_Atom>")  ):
        return (sel.index(),)
    if (repr(sel).endswith(r"p_AtomSel>") or
        repr(sel).startswith("<atomSel.AtomSel")):
        return sel.indices()

    raise Exception("error converting sel to a list of indices: " + repr(sel))


planarityDefs = {'ADE': {1: 'name N1', 3: 'name N1 or name C6 or name C2'},
                 'URI': {1: 'name N3', 3: 'name N3 or name C2 or name C4'},
                 
                 'GUA': {1: 'name N1', 3: 'name N1 or name C6 or name C2'},
                 'CYT': {1: 'name N3', 3: 'name N3 or name C2 or name C4'},

                 'PRF': {1: 'name N1', 3: 'name N1 or name C6 or name C2'}}


planarityRestraint = r'''
group
   select = ((resid %i and (%s)) or
             (resid %i and (%s)))
   weight = %s
end

'''

def genPlanarityRestraints(basePairs,xSim):
    """Given a sequence of base pairs, and an
    <m xplorSimulation>.XplorSimulation, generate planarity restraints for
    nucleic acid base pairs. Each base pair can be represented as a pair of
    residue numbers, or as a pair of tuples specifying (segid,resid).

    The restraints are returned as a string.
    """
    ret = ""
    import types
    for base1,base2 in basePairs:

        if type(base1) is int:
            sel1 = "resid %d" % base1
        else:
            sel1 = 'segid "%s" and resid %d' % base1
            pass
        if type(base2) is int:
            sel2 = "resid %d" % base2
        else:
            sel2 = 'segid "%s" and resid %d' % base2
            pass
        
        from atomSel import AtomSel
        sel1 = AtomSel(sel1,xSim)
        resname1 = sel1[0].residueName()
        resid1 = sel1[0].residueNum()

        from atomSel import AtomSel
        sel2 = AtomSel(sel2,xSim)
        resname2 = sel2[0].residueName()
        resid2 = sel2[0].residueNum()

        remark = '{%s %i   -   %s %i}' % (resname1, resid1, resname2, resid2)
        ret += remark + '\n'

        ret += planarityRestraint % (resid1,
                                     planarityDefs[resname1][1],
                                     resid2,
                                     planarityDefs[resname2][3], '30.0')

        ret += planarityRestraint % (resid1,
                                     planarityDefs[resname1][3],
                                     resid2,
                                     planarityDefs[resname2][1], '30.0')
        pass
    return ret

    
symmetricSidechains = {}

from math import pi
for (res, chiRange, chiSels, symPairs) in \
    (("PHE", (0,pi)      ,("ca", "cb", "cg", "cd1"), (("cd1", "cd2"),
                                                      ("ce1", "ce2"),
                                         ("hd1","hd2"),("he1","he2"))),
     ("TYR", (0,pi)      ,("ca", "cb", "cg", "cd1"), (("cd1", "cd2"),
                                                      ("ce1", "ce2"),
                                         ("hd1","hd2"),("he1","he2"))),
     ("ARG", (-pi/2,pi/2),("cd", "ne", "cz", "nh1"), (("nh1", "nh2"),
                                         ("hh11","hh21"),("hh12","hh22"))),
     ("ASP", (-pi/2,pi/2),("ca", "cb", "cg", "od1"), (("od1", "od2"),)),
     ("GLU", (-pi/2,pi/2),("cb", "cg", "cd", "oe1"), (("oe1", "oe2"),)),
     ):
    symmetricSidechains[res] = (chiRange, chiSels, symPairs)
    pass
    
def correctSymmetricSidechains(sel=None,
                               sim=0):
    """Correct PHE, TYR, ASP, ARG, and GLU sidechains if outside angle range.

    Step through a selection (specified via the sel argument, a selection
    string), looking for PHE, TYR, ASP, ARG, and GLU residues.  Make sure that
    each of those sidechains has its chi2 (chi3 for GLU and chi5 for ARG) angle
    in a specific range.  If not, exchange the positions of the atoms in
    symmetric locations to force that angle to be in the correct range.  For
    ARG, ASP and GLU, the angle range is taken from IUPAC rules [Biochemistry
    9:3471-3479 (1970), section 2.3.2.].  For PHE and TYR, the rules makes no
    sense as the most populated rotamer is right at the edge of the range, so
    the range 0..180 degrees is used instead.

    After the correction is made a group may still be outside the range if it's
    not completely planar.  In such case, the orientation that minimizes the
    difference between the angle and the closest range bound is kept.

    Note that protein topology/parameter sets before 2.0 have the wrong 
    definition for ARG, such that it is treated as a special case.
    """

    import dihedral
    import simulation
    from atomSel import AtomSel
    from vec3 import Vec3

    def exchange(segID, resID, nameA, nameB, sim):
        """Exchange atom positions (A <-> B)."""
        selA = AtomSel('segid "%s" and resid %d and name %s' %
               (segID,resID,nameA), sim)
        selB = AtomSel('segid "%s" and resid %d and name %s' %
               (segID,resID,nameB), sim)
        atomA = convertToAtom(selA)
        atomB = convertToAtom(selB)
        temp = Vec3(atomA.pos())
        atomA.setPos(atomB.pos())
        atomB.setPos(temp)


    if not sim:
        sim = simulation.currentSimulation()
        pass

    from xplorSimulation import getXplorSimulation    
    xSim = getXplorSimulation(sim)
    #fixup for mistake in traditional XPLOR topology for ARG
    import protocol
    
    outputState = xSim.disableOutput()
    protocol.initTopology('protein',simulation=sim)
    xSim.enableOutput(outputState)

    if ('protein' in protocol.topVersion and
        protocol.topVersion['protein']=="1.0"     ):
        symmetricSidechains['ARG'] = ((-pi/2,pi/2),
                                      ("cd", "ne", "cz", "nh2"),
                                      (("nh1", "nh2"),
                                       ("hh11","hh21"),
                                       ("hh12","hh22")))
        pass
    


    import simulation
    saveSimulation=simulation.currentSimulation()
    simulation.makeCurrent(sim)
    
    if not sel:
        sel = AtomSel("tag and (%s)"%
                      " or ".join(["resname %s" %resname for resname in
                                   list(symmetricSidechains.keys())]),
                      sim)
        pass

    segments = getSegsResidues(sel,sim)

    for segID in list(segments.keys()):
        for (resID, resName) in segments[segID]:
            if not resName in list(symmetricSidechains.keys()):
                continue
            (chiRange, chiSels, symPairs) = symmetricSidechains[resName]
            
            sel0 = 'segid "%s" and resid %d and name %s' % (segID,resID,
                                                            chiSels[0])
            sel1 = 'segid "%s" and resid %d and name %s' % (segID,resID,
                                                            chiSels[1])
            sel2 = 'segid "%s" and resid %d and name %s' % (segID,resID,
                                                            chiSels[2])
            sel3 = 'segid "%s" and resid %d and name %s' % (segID,resID,
                                                            chiSels[3])
            try:
                d =  dihedral.Dihedral(sel0, sel1, sel2, sel3)
                chi = d.value()
                if chi <= chiRange[0] or chi > chiRange[1] :
                    if chi <= chiRange[0]: diff = abs(chi - chiRange[0])
                    if chi > chiRange[1]: diff = abs(chi - chiRange[1])
                    for (nameA, nameB) in symPairs:
                        try:
                            exchange(segID, resID, nameA, nameB, sim)
                            
                            # Account for non-planar group.
                            d = dihedral.Dihedral(sel0, sel1, sel2, sel3)
                            chi = d.value()
                            if chi <= chiRange[0] or chi > chiRange[1]:
                                if chi <= chiRange[0]: newdiff =abs(chi-chiRange[0])
                                if chi > chiRange[1]: newdiff = abs(chi-chiRange[1])
                                if newdiff > diff:  # revert
                                    exchange(segID, resID, nameA, nameB, sim)
                                    pass
                                pass
                            pass
                        except:
                            #this allows for missing atoms
                            pass
                        pass
                    pass
                pass
            except IndexError:
                #missing atoms
                pass
            pass
        pass
    simulation.makeCurrent(saveSimulation)
    return 

def makeAtomSelStringFromResidList(resList):
    """
    Make an Xplor atom selection string from a list of residue numbers.
    """
    resList.sort()

    listLen=len(resList)

    s=''
    ret=''
    if listLen > 0:
     r_s=resList[0]
     r_st=r_s
     k=1
     while k < listLen:
        r_c=resList[k]
        if r_c > r_st+1:
           if r_s==r_st:
              s += ' resid %i'%r_s
           else:
              s += ' resid %i'%r_s+':%i'%r_st
              pass
           r_s=r_c
           r_st=r_c
        else:
           r_st=r_c
           pass
        k+=1
        pass

     if r_s==r_st:
        s += ' resid %i'%r_s
     else:
        s += ' resid %i'%r_s+':%i'%r_st
        pass

     ret=s[1]

     count=0    
     k=2
     pos_in=0
     while k < len(s):   
        if s[k]=='r':
          count+=1
          st=ret[0:len(ret)-1]
          ret=st+')'
          ret += ' or ('
          pass
        ret += s[k]
        k+=1
        pass
 
     if count > 0:
       ret='('+ret+')'

    return ret

def breakAtomSelString(selString,
                       lineLen=80,
                       spacer="",
                       initialShift=0):
    """
    Break the line obtained from makeAtomSelStringFromResidList into
       several lines.
    """

    k=0
    n_lines=-1
    while k < len(selString):
        if (selString[k]=='o' and
            ((len(selString[0:k])-
              n_lines+initialShift)/lineLen)-n_lines > 1):
            selString=selString[0:k-1]+'\n'+spacer+selString[k:len(selString)]
            k += len(spacer)+1
            n_lines += 1
            pass
        k=k+1
        pass  
    
    return selString

def fixMethylImpropers(sel="all",verbose=False):
    """
    Swap proton positions on methyls such that improper energies are lowered.

    Swapped coordinates include any C, N or S bound to three protons.
    """
    sel=convertToAtomSel(sel)
    from atomSel import intersection, AtomSel
    from xplorPot import XplorPot
    simulation=sel.simulation()
    impr=XplorPot("IMPR",simulation)
    Eimpr=impr.calcEnergy()
    for heavy in intersection(sel,AtomSel("name C* or name S* or name N*",
                                          simulation)):
        protons=AtomSel('''bondedto ATOM "%s" %d %s
                        and name H*''' %(heavy.segmentName(),
                       heavy.residueNum(),
                       heavy.atomName()),
                        simulation)
        if len(protons)==3:
            coords=[atom.pos() for atom in protons]
            protons[0].setPos(coords[1])
            protons[1].setPos(coords[0])
            new=impr.calcEnergy()
            if new<Eimpr:
                if verbose:
                    print("fixMethylImpropers: swapped  %s and %s" % \
                          (protons[0].string(), protons[1].string()))
                Eimpr=new
            else:
                protons[0].setPos(coords[0])
                protons[1].setPos(coords[1])
                pass
            pass
        pass
    return

def fixProtonImpropers(sel="all",verbose=False):
    """
    Swap proton positions on methyls or methylenes such that improper
    energies are lowered. 

    Swapped coordinates include any C, N or S bound to two or three protons.
    """
    sel=convertToAtomSel(sel)
    from atomSel import intersection, AtomSel
    from xplorPot import XplorPot
    simulation=sel.simulation()
    impr=XplorPot("IMPR",simulation)
    Eimpr=impr.calcEnergy()
    for heavy in intersection(sel,AtomSel("name C* or name S* or name N*",
                                          simulation)):
        protons=AtomSel('''bondedto ATOM "%s" %d %s
                        and name H*''' %(heavy.segmentName(),
                       heavy.residueNum(),
                       heavy.atomName()),
                        simulation)
        if len(protons)==2 or len(protons)==3:
            coords=[atom.pos() for atom in protons]
            protons[0].setPos(coords[1])
            protons[1].setPos(coords[0])
            new=impr.calcEnergy()
            if new<Eimpr:
                if verbose:
                    print("fixMethylImpropers: swapped  %s and %s" % \
                          (protons[0].string(), protons[1].string()))
                Eimpr=new
            else:
                protons[0].setPos(coords[0])
                protons[1].setPos(coords[1])
                pass
            pass
        pass
    return
            


            
def toAtomSelString(atom):
    """
    Given an atom, return a selection string.
    """
    ret = 'atom "%s" %d %s' % (atom.segmentName(),
                               atom.residueNum(),
                               atom.atomName()   )
    return ret
