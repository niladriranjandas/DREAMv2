"""
The residue affinity potential encodes the arrangement of hydrophobic, polar
and charged residues in a simple potential term.

It can be used as a scoring function, or in refinement.

This term is described in
  Y. Ryabov, J.-Y.  Suh, A. Grishaev,  G.M. Clore and C.D.
  Schwieters, ``Using the experimentally determined
  components of the overall rotational diffusion tensor to restrain
  molecular shape and size in NMR structure determination of globular
  proteins and protein-protein complexes}{J. Am. Chem. Soc. 131, 9522-9531
  (2009)
"""

def create_ResidueAffPot(name,sel="known",
                         sequentialCutoff=3,
                         intradomainContacts=True,
                         interdomainContacts=False,
                         selectionPairs=[],
                         potentialName="Thomas"):
    """
    The residue affinity potential encodes the arrangement of hydrophobic,
    polar and charged residues in a simple potential term.

    If an atom selection is given, only the specified residues will be
    included in the term.

    sequentialCutoff specifies the neighborhood in sequence about which
    the potential term is not evaluated.

    This potential term is implemented as a <m selNBPot>.SelNBPot.

    The flags interdomainContacts and intradomainContacts control whether or
    not interactions between domains are computed. Domain identity is
    determined by segid.

    Alternatively, active interactions can be specified using the selectionPairs
    argument, a list of pairs of atom selections, the interactions
    between which are not nonzero. If selectionPairs is specified,
    intradomainContacts and interdomainContacts arguments are ignored.

    By default the 5-class contact potential parameters are used from
    Thomas and Dill,
    Proc. Natl. Acad. Sci. USA 93, 11628-11633 (1996).

    Also available is the potential from
    S. Miyazawa and R.L. Jernigan, PROTEINS 34, 49-68 (1999).
    This can be used by specifying potentialName='Miyazawa' 

    """
    from atomSel import AtomSel
    if type(sel)==type(""):
        sel = AtomSel(sel)
        pass


    import selNBPot
    groups = []
    classes = []

    try:
        residue   = residues[potentialName]
        potential = potentials[potentialName]
    except KeyError:
        raise Exception("invalid potentialName specified")

    segidResidMapping = {}
    for atom in AtomSel("tag and (%s)" % sel.string(),
                        sel.simulation()             ):
        resid = atom.residueNum()
        segid = atom.segmentName()
        resName = atom.residueName()
        if not resName in list(residue.keys()):
            print("create_residueAffPot: residue %s not known" % resName)
            continue
        sstring = 'resid %d and segid "%-4s" and (' % (resid,segid)
        (atoms,classID) = residue[resName]
        sstring += "name " + atoms[0]
        if len(atoms)>1:
            for a in atoms[1:]:
                sstring += " or name " + a
                pass
            pass
        sstring += ")"
        if segid not in segidResidMapping: segidResidMapping[segid] = {}
        segidResidMapping[segid][resid] = len(groups)
        groups.append(AtomSel(sstring,sel.simulation()))
        if len(groups[-1])==0:
            raise Exception("no atoms for selection >%s<" % sstring)
        classes.append(classID)
        pass

    
    from atomSel import intersection
    from selectTools import convertToAtomSel
    for i,(a,b) in enumerate(selectionPairs):
        a,b = (convertToAtomSel(a), convertToAtomSel(b))
        a = [segidResidMapping[atom.segmentName()][atom.residueNum()]
             for atom in intersection("tag",a)]
        b = [segidResidMapping[atom.segmentName()][atom.residueNum()]
             for atom in intersection("tag",b)]
        selectionPairs[i] = (a,b)  # converted to group indices
        pass        

    def getPot(i,j):
        return potential[max(i,j)][min(i,j)]

    from cdsMatrix import SymMatrix_double as SymMatrix
    intMat=SymMatrix(len(groups))
    intMat.set(0.)
    
    from atomSel import minResidDiff
    if len(selectionPairs):
        for aGroups,bGroups in selectionPairs:
            for a in aGroups:
                for b in bGroups:
                    groupA = groups[a]
                    groupB = groups[b]
                    if minResidDiff(groupA,groupB)<= sequentialCutoff:
                        continue
                    intMat[a,b] = getPot(classes[a],classes[b])
                    pass
                pass
            pass
        pass
    else:
        for i in range(len(groups)):
            for j in range(i):
                atom_i = groups[i][0]
                atom_j = groups[j][0]

                isInteracting  = ((atom_i.segmentName()==atom_j.segmentName() and
                                   intradomainContacts)
                                  or
                                  (atom_i.segmentName()!=atom_j.segmentName() and
                                   interdomainContacts))
            
                if isInteracting and \
                   abs(atom_i.residueNum()-
                       atom_j.residueNum() ) > sequentialCutoff:
                    intMat[i,j] = getPot(classes[i],classes[j])
                    pass
                pass
            pass
        pass
    

    pot = selNBPot.SelNBPot(name,groups,intMat)

    

    return pot

residues={}
potentials={}

# 5-class
# contact potential parameters from Thomas and Dill, Proc. Natl. Acad.
# Sci. USA 93, 11628-11633 (1996).
#

residue={}
#residue          active atoms                    class
residue['GLY'] = (("CA",)  ,                        1)       #   G
residue['ALA'] = (("CB",)  ,                        0)       #   A
residue['VAL'] = (("CG*",) ,                        0)       #   V
residue['PHE'] = (("CD*","CE*","CZ") ,              0)       #   F
residue['PRO'] = (("CB","CG","CD") ,                1)       #   P
residue['MET'] = (("SD","CE") ,                     0)       #   M
residue['ILE'] = (("CG2","CD1") ,                   0)       #   I
residue['LEU'] = (("CD*",) ,                        0)       #   L
residue['ASP'] = (("OD*",) ,                        3)       #   D
residue['GLU'] = (("OE*",) ,                        3)       #   E
residue['LYS'] = (("NZ",) ,                         4)       #   K
residue['ARG'] = (("NH*",) ,                        4)       #   R
residue['SER'] = (("OG*",) ,                        1)       #   S
residue['THR'] = (("OG1","CG2",) ,                  1)       #   T
residue['TYR'] = (("CD*","CE*","CZ","OH") ,         0)       #   Y
residue['HIS'] = (("ND1","CD2","CE1","NE2") ,       1)       #   H
residue['CYS'] = (("SG",) ,                         2)       #   C
residue['ASN'] = (("OD1","ND2") ,                   1)       #   N
residue['GLN'] = (("OE1","NE2") ,                   1)       #   Q
residue['TRP'] = (("CD*","NE1","CE*","CZ*","CH2") , 0)       #   W

#potential
potential=[[-0.65,                          ],
           [ 0.04,  0.05,                   ],
           [-0.68, -0.26, -3.25,            ],
           [ 0.32,  0.12,  1.42,  0.61,     ],
           [-0.17,  0.17, -0.09, -0.43, 1.91]]

from cdsMatrix import CDSMatrix_double, transpose
residues["Thomas"]   = residue
potentials["Thomas"] = potential

#
# 20 group contact energies from S. Miyazawa and R.L. Jernigan,
# PROTEINS 34, 49-68 (1999).
#

res={}
#  residue          active atoms                                      class
res['GLY'] = (("CA",)                        , 9)
res['ALA'] = (("CB",)                        , 8)
res['VAL'] = (("CG*",)                       , 5)
res['PHE'] = (("CD*","CE*","CZ")             , 2)
res['PRO'] = (("CB","CG","CD")               ,19)
res['MET'] = (("SD","CE")                    , 1)
res['ILE'] = (("CG2","CD1")                  , 3)
res['LEU'] = (("CD*",)                       , 4)
res['ASP'] = (("OD*",)                       ,15)
res['GLU'] = (("OE*",)                       ,14)
res['LYS'] = (("NZ",)                        ,18)
res['ARG'] = (("NH*",)                       ,17)
res['SER'] = (("OG*",)                       ,11)
res['THR'] = (("OG1","CG2",)                 ,10)
res['TYR'] = (("CD*","CE*","CZ","OH")        , 7)
res['HIS'] = (("ND1","CD2","CE1","NE2")      ,16)
res['CYS'] = (("SG",)                        , 0)
res['ASN'] = (("OD1","ND2")                  ,13)
res['GLN'] = (("OE1","NE2")                  ,12)
res['TRP'] = (("CD*","NE1","CE*","CZ*","CH2"), 6)

# these values are in units of J/mol (RT units) in the paper, but are used
# in kcal/mol here, so an implicit scale factor has been introduced.

# lower triangle is the contact energy which includes solvent (it has been transposed from table 5 in the paper)
#
# CYS ,  MET ,  PHE ,  ILE ,  LEU ,  VAL ,  TRP ,  TYR ,  ALA ,  GLY ,  THR ,  SER ,  GLN ,  ASN ,  GLU ,  ASP ,  HIS ,  ARG ,  LYS ,  PRO
potential_Miyazawa=[
[-1.19, -0.54, -0.03, -0.02, -0.03,  0.00, -0.03, -0.06,  0.16,  0.01, -0.03,  0.12,  0.08,  0.10,  0.21,  0.31,  0.26, -0.01,  0.33,  0.35],  # CYS
[-0.61, -0.7 , -0.19, -0.24, -0.13, -0.12, -0.03, -0.21, -0.08,  0.00, -0.10,  0.09,  0.19,  0.04,  0.19,  0.16,  0.37, -0.02,  0.21,  0.24],  # MET
[-0.67, -0.83, -0.88, -0.22, -0.11, -0.14, -0.11, -0.08, -0.02, -0.02, -0.16,  0.13,  0.10,  0.07,  0.22,  0.26,  0.33,  0.00,  0.20,  0.22],  # PHE
[-0.64, -0.66, -0.73, -0.74, -0.18, -0.20, -0.16, -0.05,  0.02, -0.07, -0.17,  0.08,  0.20,  0.12,  0.32,  0.24,  0.32,  0.18,  0.21,  0.23],  # ILE
[-0.65, -0.7 , -0.80, -0.81, -0.84, -0.19, -0.19, -0.03,  0.00, -0.04, -0.18,  0.13,  0.20,  0.14,  0.26,  0.28,  0.41,  0.17,  0.21,  0.24],  # LEU
[-0.59, -0.51, -0.67, -0.67, -0.74, -0.65, -0.20, -0.02,  0.08, -0.07, -0.09,  0.10,  0.16,  0.15,  0.25,  0.27,  0.41,  0.19,  0.22,  0.22],  # VAL
[-0.66, -0.73, -0.68, -0.60, -0.62, -0.51, -0.64, -0.11,  0.01,  0.01, -0.03,  0.19,  0.15,  0.09,  0.06,  0.05,  0.15, -0.08, -0.03,  0.05],  # TRP
[-0.39, -0.56, -0.58, -0.49, -0.55, -0.38, -0.49, -0.45,  0.00,  0.05, -0.03,  0.08,  0.04, -0.06,  0.02, -0.07, -0.03, -0.05, -0.10, -0.12],  # TYR
[-0.33, -0.27, -0.36, -0.37, -0.38, -0.32, -0.27, -0.20, -0.12, -0.09, -0.05,  0.00,  0.00,  0.08,  0.07,  0.18,  0.10,  0.10,  0.18,  0.13],  # ALA
[-0.31, -0.17, -0.19, -0.13, -0.16, -0.15, -0.25, -0.22, -0.08, -0.29, -0.26, -0.08, -0.10, -0.01, -0.10,  0.13, -0.06,  0.04,  0.02,  0.01],  # GLY
[-0.15, -0.11, -0.15, -0.15, -0.15, -0.07, -0.02, -0.09,  0.04, -0.04, -0.03, -0.08, -0.13, -0.09, -0.12, -0.11, -0.13, -0.06, -0.03, -0.02],  # THR
[-0.13,  0.05, -0.12,  0.03, -0.02,  0.04, -0.01, -0.08,  0.10, -0.01, -0.04,  0.05, -0.17, -0.04, -0.12, -0.14, -0.20, -0.06, -0.04, -0.05],  # SER
[-0.07, -0.06, -0.11, -0.01, -0.04,  0.08, -0.02, -0.14,  0.22,  0.13, -0.12,  0.22,  0.20, -0.10, -0.20, -0.10, -0.10,  0.02, -0.14, -0.17],  # GLN
[-0.01,  0.04, -0.01,  0.14,  0.04,  0.12, -0.10, -0.11,  0.15, -0.01, -0.04,  0.09,  0.06, -0.06, -0.27, -0.19, -0.27, -0.08, -0.08, -0.18],  # ASN
[ 0.20,  0.12,  0.14,  0.17,  0.17,  0.26,  0.00, -0.08,  0.38,  0.32, -0.16,  0.18,  0.27,  0.12,  0.46,  0.03,  0.04, -0.19, -0.52, -0.58],  # GLU
[ 0.12,  0.3 ,  0.18,  0.22,  0.27,  0.36,  0.07, -0.07,  0.27,  0.11, -0.11,  0.10,  0.24,  0.02,  0.44,  0.29, -0.08, -0.26, -0.51, -0.49],  # ASP
[-0.36, -0.29, -0.34, -0.13, -0.18, -0.06, -0.37, -0.30,  0.07,  0.00, -0.03,  0.04,  0.15,  0.00,  0.00, -0.10, -0.40, -0.36, -0.02,  0.10],  # HIS
[ 0.08,  0.03, -0.05,  0.00, -0.04,  0.08, -0.21, -0.25,  0.24,  0.09, -0.11,  0.16,  0.09,  0.10, -0.22, -0.24,  0.05,  0.19,  0.03,  0.28],  # ARG
[ 0.33,  0.29,  0.19,  0.24,  0.22,  0.29,  0.09, -0.05,  0.41,  0.29, -0.33,  0.36,  0.28,  0.22, -0.06, -0.01,  0.38,  0.66,  0.76,  0.16],  # LYS
[-0.18, -0.13, -0.19, -0.05, -0.12, -0.05, -0.37, -0.25,  0.15,  0.02, -0.13,  0.20,  0.17,  0.18,  0.37,  0.33,  0.01,  0.17,  0.47,  0.11]]  # PRO


residues["Miyazawa"] = res
potentials["Miyazawa"] = potential_Miyazawa

def analyze(potList):
    """perform analysis of ResidueAffPot terms and return nicely formatted
    summary"""

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'SelNBPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()


    ret += "%-9s  %6s  %6s   %6s \n" % \
           ( "", "contacts", "attractive", "repulsive")

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]
#
#        print term.showViolations()
        print(term.info())
        #
        d=term.distances()
        groups=term.groups()
        intMat=term.interactionMat()
        contacts=0
        attractive=0
        repulsive=0
        for i in range(len(groups)):
            for j in range(i):
                if d[i,j]>=0:
                    contacts += 1
                    if intMat[i,j]<0: attractive += 1
                    if intMat[i,j]>0: repulsive += 1
                    pass
                pass
            pass
        ret += "%-9s  %6d     %6d     %6d\n" % \
               (name , contacts, attractive, repulsive  )
        pass
    
    return ret


from simulationTools import registerTerm
registerTerm(analyze,"Residue Affinity terms","ResidueAffinity",
r"""
For each term, report the total number of residue-residue contacts,
and the number of these which are attractive and repulsive. See
<m residueAffPotTools>.
""")
    
