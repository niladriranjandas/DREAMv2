
"""
  Tools to help create/analyze potential terms based on the difference between
  a glob position computed from atomic coordinates and that supplied by a
  bead model calculation.
"""

def create_GlobDiffPot(name,
                       selection="not PSEUDO and not name H*",
                       globSelection=None,
                       globFilename=None,
                       psf=None,
                       ):
    """Create a GlobDiffPot term based on the specified atom selection.

    The name of the GlobDiffPot term is given by the name argument (string).
    If globFilename is specified, the comparison is made with the coordinates
    in that file (using globSelection).  It is an error to omit both the
    globSelection and globFilename arguments.  For this potential term, the
    order of globs in globSelection must match those generated from selection.

    This potential term is based on the <m posSymmPot>.PosSymmPot potential,
    where glob coordinates are generated from atomic positions by weighting
    each position by the corresponding atom's number of electrons.

    If globFile (a string with a PDB filename) is specified, by default the
    associated coordinates are loaded using <m protocol>.loadPDB (with
    deleteUnknownAtoms=True); this assumes a standard topology.  Non-standard
    topology for the input PDB file can be accomodated with the psf argument,
    which can be a string with the contents of a PSF file, a string with a PSF
    filename, or a sequence with one of more PSF filenames.  If psf is
    specified, the coordinates are loaded with <m protocol>.initCoords, using
    deleteUnknownAtoms=True.

    The returned potential term has four members beyond those in a
     <m posSymmPot>.PosSymmPot term:

        globFilename  - the value of the globFilename argument.
        selection     - an <m atomSel>.AtomSel corresponding to the selection
                        argument.
        globSelection - an <m atomSel>.AtomSel corresponding to the
                        globSelection argument.
        xsim          - a <m xplorSimulation>.XplorSimulation generated for
                        containing glob position information if globFile is
                        specified.
    
    """

    if not globSelection: globSelection = selection
    
    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)

    sim = selection.simulation()

    xsim=None
    if globFilename:
        from xplorSimulation import XplorSimulation, getXplorSimulation
        from simulation import currentSimulation, makeCurrent

        #hack s.t. currentSimulation is an XplorSimulation at creation time
        cursim = currentSimulation()
        xsim1=getXplorSimulation(cursim)
        makeCurrent(xsim1)
        xsim = XplorSimulation()
        makeCurrent(cursim)
        
        import protocol
        if not psf:  # assume standard topology
            protocol.loadPDB(globFilename,simulation=xsim,deleteUnknownAtoms=True)
        else:  # read input psf
            from atomSel import AtomSel
            protocol.initStruct(psf, simulation=xsim)
            protocol.initCoords(globFilename, deleteUnknownAtoms=True,
                                selection=AtomSel('all', xsim))

        globSelection = convertToAtomSel(globSelection,xsim)
        pass
    else:
        globSelection = convertToAtomSel(globSelection)
        pass

    print('num globs:', len(globSelection))

    if len(selection) < 1 :
        print("Warning: empty selection. String:", selection.string())
        #raise Exception("atomselection does not select any atoms")
    
    from posSymmPot import PosSymmPot

    pot = PosSymmPot(name)
    pot.resetPotName( "GlobDiffPot" )

    pot.selection = selection
    pot.globSelection = globSelection
    pot.globFilename = globFilename
    pot.xsim = xsim #need to keep this simulation around

    from selectTools import getSegsResidues
    segsResidues = getSegsResidues(selection)
    globSegsResidues = getSegsResidues(globSelection)


    sels = []
    from atomSel import AtomSel, intersection
    for segid in segsResidues:
        for resid,resname in segsResidues[segid]:
            sels.append( intersection(selection,
                                      AtomSel('segid "%s" and resid %d' %
                                              (segid,resid),
                                              selection.simulation()    )) )
            pass
        pass

    globSels = []
    for segid in globSegsResidues:
        for resid,resname in globSegsResidues[segid]:
            sel = AtomSel('segid "%s" and resid %d' % (segid,resid),
                          globSelection.simulation()                )
            if len(sel)!=1:
                raise Exception("create_GlobDiffPot: bad number of globs in"+
                                " segid %s resid %d" % (segid,resid))
            globSels.append( sel )
            pass
        pass

    if len(sels) != len(globSels):
        raise Exception("size mismatch between residues in selection (%d)" %
                        len(sels) + " and number of globs (%d)" %
                        len(globSels))

    for sel,globSel in zip(sels,globSels):
        resName = sel[0].residueName()
        weights=[]
        sum=0
        for atom in sel:
            electrons = electronCounts[resName][atom.atomName()][2]
            sum += electrons
            weights.append( electrons )
            pass
        from cdsVector import CDSVector_double
        weights = CDSVector_double(weights)
        weights.scale(1./sum)
        pot.addEquivAtomSelPair(sel,globSel,weights)
        pass

    pot.setThreshold(0.01) #default threshold value

    return pot

def analyze(potList):
    """perform analysis of GlobDiffPot terms and return nicely formatted
    summary"""

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'GlobDiffPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    
    hasPDB=False
    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        if term.globFilename: hasPDB=True
        pass
        

    ret += "%-9s  %6s   %8s" % \
           ( "", "RMSD", "Delta Variance")
    ret += '\n'

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        term.calcEnergy()

        
        print("term:", term.instanceName())
        print(" lowerBound:", term.lowerBound())
        print(" upperBound:", term.upperBound())
        print(" scale:     ", term.scale())
        if term.globFilename:
            print("comparison PDB:", term.globFilename)
            pass
        print("selection string:", term.selection.string())
        print(" (Simulation: %s)" % term.selection.simulation().name())
        print("globSelection string:", end=' ') 
        if term.selection.string() == term.globSelection.string():
            print("same")
        else:
            print(term.globSelection.string())
            pass
        print(" (Simulation: %s)" % term.globSelection.simulation().name())
        print("calculated RMSD:", term.rms())

        import mat3
        print("fit rotation:")
        print(term.rotMat())
        print("group centroid positions:")
        print("  0: ", term.centroid(0), "   1: ", term.centroid(1))
 
        ret += "%-9s  %7.3f  %8.3f" % \
               (name , term.rms(), term.deltaVariance())
        ret += '\n'
        pass
    
    return ret


from simulationTools import registerTerm
registerTerm(analyze,"Glob Difference Potential","GlobDiff",
r"""
For each term, report the root mean square difference between calculated and
observed value, and the matrix norm of the difference in the variance tensor.
Please see <m posSymmPot>.
""")
    

        
    
        

electronCounts={
    'ADE': {'P': ('P', 'P', 15.0, 'A', 'NUC'), 
            'O1P': ('OP1', 'O', 8.0, 'A', 'NUC'), 
            'O2P': ('OP2', 'O', 8.0, 'A', 'NUC'), 
            "O5'": ("O5'", 'O', 8.0, 'A', 'NUC'), 
            "C5'": ("C5'", 'CH2', 8.0, 'A', 'NUC'),
            "C4'": ("C4'", 'CH', 7.0, 'A', 'NUC'),
            "O4'": ("O4'", 'O', 8.0, 'A', 'NUC'), 
            "C3'": ("C3'", 'CH', 7.0, 'A', 'NUC'), 
            "O3'": ("O3'", 'O', 8.0, 'A', 'NUC'), 
            "C2'": ("C2'", 'CH', 7.0, 'A', 'NUC'),
            "O2'": ("O2'", 'OH', 9.0, 'A', 'NUC'), 
            "C1'": ("C1'", 'CH', 7.0, 'A', 'NUC'), 
            "N9": ("N9", 'N', 7.0, 'A', 'NUC'), 
            "C8": ("C8", 'CH', 7.0, 'A', 'NUC'),
            "N7": ("N7", 'N', 7.0, 'A', 'NUC'),
            "C5": ("C5", 'C', 6.0, 'A', 'NUC'), 
            "C6": ("C6", 'C', 6.0, 'A', 'NUC'), 
            "N6": ("N6", 'NH2', 9.0, 'A', 'NUC'), 
            "N1": ("N1", 'N', 7.0, 'A', 'NUC'),
            "C2": ("C2", 'CH', 7.0, 'A', 'NUC'), 
            "N3": ("N3", 'N', 7.0, 'A', 'NUC'), 
            "C4": ("C4", 'C', 6.0, 'A', 'NUC')},

    'CYT': {'P': ('P', 'P', 15.0, 'C', 'NUC'), 
            'O1P': ('OP1', 'O', 8.0, 'C', 'NUC'), 
            'O2P': ('OP2', 'O', 8.0, 'C', 'NUC'), 
            "O5'": ("O5'", 'O', 8.0, 'C', 'NUC'), 
            "C5'": ("C5'", 'CH2', 8.0, 'C', 'NUC'),
            "C4'": ("C4'", 'CH', 7.0, 'C', 'NUC'),
            "O4'": ("O4'", 'O', 8.0, 'C', 'NUC'), 
            "C3'": ("C3'", 'CH', 7.0, 'C', 'NUC'), 
            "O3'": ("O3'", 'O', 8.0, 'C', 'NUC'), 
            "C2'": ("C2'", 'CH', 7.0, 'C', 'NUC'),
            "O2'": ("O2'", 'OH', 9.0, 'C', 'NUC'), 
            "C1'": ("C1'", 'CH', 7.0, 'C', 'NUC'), 
            "N1": ("N1", 'N', 7.0, 'C', 'NUC'), 
            "C2": ("C2", 'C', 6.0, 'C', 'NUC'),
            "O2": ("O2", 'O', 8.0, 'C', 'NUC'),
            "N3": ("N3", 'N', 7.0, 'C', 'NUC'), 
            "C4": ("C4", 'C', 6.0, 'C', 'NUC'), 
            "N4": ("N4", 'NH2', 9.0, 'C', 'NUC'), 
            "C5": ("C5", 'CH', 7.0, 'C', 'NUC'),
            "C6": ("C6", 'CH', 7.0, 'C', 'NUC')}, 
 
    'GUA': {'P': ('P', 'P', 15.0, 'G', 'NUC'), 
            'O1P': ('OP1', 'O', 8.0, 'G', 'NUC'), 
            'O2P': ('OP2', 'O', 8.0, 'G', 'NUC'), 
            'OP3': ('OP3', 'O', 8.0, 'G', 'NUC'), 
            "O5'": ("O5'", 'O', 8.0, 'G', 'NUC'), 
            "C5'": ("C5'", 'CH2', 8.0, 'G', 'NUC'),
            "C4'": ("C4'", 'CH', 7.0, 'G', 'NUC'),
            "O4'": ("O4'", 'O', 8.0, 'G', 'NUC'), 
            "C3'": ("C3'", 'CH', 7.0, 'G', 'NUC'), 
            "O3'": ("O3'", 'O', 8.0, 'G', 'NUC'), 
            "C2'": ("C2'", 'CH', 7.0, 'G', 'NUC'),
            "O2'": ("O2'", 'OH', 9.0, 'G', 'NUC'), 
            "C1'": ("C1'", 'CH', 7.0, 'G', 'NUC'), 
            "N9": ("N9", 'N', 7.0, 'G', 'NUC'), 
            "C8": ("C8", 'CH', 7.0, 'G', 'NUC'),
            "N7": ("N7", 'N', 7.0, 'G', 'NUC'),
            "C5": ("C5", 'C', 6.0, 'G', 'NUC'), 
            "C6": ("C6", 'C', 6.0, 'G', 'NUC'), 
            "O6": ("O6", 'O', 8.0, 'G', 'NUC'), 
            "N1": ("N1", 'NH', 8.0, 'G', 'NUC'),
            "C2": ("C2", 'C', 6.0, 'G', 'NUC'), 
            "N2": ("N2", 'NH2', 9.0, 'G', 'NUC'), 
            "N3": ("N3", 'N', 7.0, 'G', 'NUC'),
            "C4": ("C4", 'C', 6.0, 'G', 'NUC')},
    
    'URI': {'P': ('P', 'P', 15.0, 'U', 'NUC'), 
            'O1P': ('OP1', 'O', 8.0, 'U', 'NUC'), 
            'O2P': ('OP2', 'O', 8.0, 'U', 'NUC'), 
            "O5'": ("O5'", 'O', 8.0, 'U', 'NUC'), 
            "C5'": ("C5'", 'CH2', 8.0, 'U', 'NUC'),
            "C4'": ("C4'", 'CH', 7.0, 'U', 'NUC'),
            "O4'": ("O4'", 'O', 8.0, 'U', 'NUC'), 
            "C3'": ("C3'", 'CH', 7.0, 'U', 'NUC'), 
            "O3'": ("O3'", 'O', 8.0, 'U', 'NUC'), 
            "C2'": ("C2'", 'CH', 7.0, 'U', 'NUC'),
            "O2'": ("O2'", 'OH', 9.0, 'U', 'NUC'), 
            "C1'": ("C1'", 'CH', 7.0, 'U', 'NUC'), 
            "N1": ("N1", 'N', 7.0, 'U', 'NUC'), 
            "C2": ("C2", 'C', 6.0, 'U', 'NUC'),
            "O2": ("O2", 'O', 8.0, 'U', 'NUC'),
            "N3": ("N3", 'NH', 8.0, 'U', 'NUC'), 
            "C4": ("C4", 'C', 6.0, 'U', 'NUC'), 
            "O4": ("O4", 'O', 8.0, 'U', 'NUC'), 
            "C5": ("C5", 'CH', 7.0, 'U', 'NUC'),
            "C6": ("C6", 'CH', 7.0, 'U', 'NUC')}, 
    }

