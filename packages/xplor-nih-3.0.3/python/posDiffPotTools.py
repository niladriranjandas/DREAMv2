
"""
  Tools to help create/analyze potential terms based on atomic coordinate
  difference
"""

def create_PosDiffPot(name,
                      selection,
                      selection2=None,
                      pdbFile=None,
                      psf=None,
                      absoluteFit=False,
                      cmpSel=None):
    """Create a PosDiffPot term based on the specified atom selection.

    The name of the PosDiffPot term is given by the name argument (string).
    If selection2 is not specified, it is taken to be the same as selection.
    If pdbFile is specified, the comparison is made with the coordinates in that
    file (using selection2).  It is an error to omit both the selection2 and
    pdbFile arguments.  For this potential term, the order of atom indices in
    selection2 must match that in selection. If selection is an ordered
    <m atomSel>.AtomSel object, then selection2 will also be ordered if omitted
    or specified by string.

    This potential term is based on the <m posSymmPot>.PosSymmPot potential.
    In this case each atom is placed in its own equivalent set (as in the
    NCS potential).

    If absoluteFit=True, selection and selection2 will be physically overlaid,
    as rigid body translation and rotation will not be performed before the
    fit.

    cmpSel is a selection or pair of selections (the later is required if
    selection2 is supplied) to use for an additional rmsd computation.  It is
    not used in an energy calculation.

    If pdbFile (a string with a PDB filename) is specified, by default the
    associated coordinates are loaded using <m protocol>.loadPDB (with
    deleteUnknownAtoms=True); this assumes a standard topology.  A non-standard
    topology for the input PDB file can be accomodated with the psf argument,
    which can be a string with the contents of a PSF file, a string with a PSF
    filename, or a sequence with one of more PSF filenames.  If psf is
    specified, the coordinates are loaded with <m protocol>.initCoords, using
    deleteUnknownAtoms=True.

    The returned potential term has four members beyond those in a
     <m posSymmPot>.PosSymmPot term:

        pdbFile     - the value of the pdbFile argument.
        selection   - an <m atomSel>.AtomSel corresponding to the selection
                      argument.
        selection2  - an <m atomSel>.AtomSel corresponding to the selection2
                      argument.
        xsim        - a <m xplorSimulation>.XplorSimulation generated for
                      containing atomic information if pdbFile is
                      specified.
    
    """

    if not selection2: selection2 = selection
    
    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)

    sim = selection.simulation()

    xsim=None
    if pdbFile:
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
            protocol.loadPDB(pdbFile,simulation=xsim,deleteUnknownAtoms=True)
        else:  # read input psf
            from atomSel import AtomSel
            protocol.initStruct(psf, simulation=xsim)
            protocol.initCoords(pdbFile, deleteUnknownAtoms=True,
                                selection=AtomSel('all', xsim))
            pass

        selection2 = convertToAtomSel(selection2,xsim,
                                      selection.ordered())
        pass
    else:
        selection2 = convertToAtomSel(selection2,
                                      ordered=selection.ordered())
        pass

    if len(selection) < 1 :
        print("Warning: empty selection. String:", selection.string())
        #raise Exception("atomselection does not select any atoms")
    
    if len(selection) != len(selection2):
        print("measure selections of different size: %d != %d" %(
            len(selection), len(selection2)))
        for i,atom in enumerate(selection):
            j=i
            try:
                atom2=selection2[j]
            except IndexError:
                print("no matching atom in selection2:", atom.string())
                break
            if (atom.atomName()    != atom2.atomName() #or
#                atom.residueName() != atom2.residueName()
                ):
                print("  possible mismatch: %s != %s" % (atom.string(),
                                                         atom2.string()))
                break
            pass
        raise Exception("len(selection) [%d] != len(selection2) [%d]" %
                        (len(selection), len(selection2)))

    from posSymmPot import PosSymmPot

    pot = PosSymmPot(name)
    pot.resetPotName( "PosDiffPot" )
    pot.setAbsoluteFit( absoluteFit )

    pot.selection = selection
    pot.selection2 = selection2
    pot.pdbFile = pdbFile
    pot.xsim = xsim #need to keep this simulation around

    if cmpSel:
        try:
            cmpSel=[convertToAtomSel(cmpSel,selection.simulation()),
                    convertToAtomSel(cmpSel,selection2.simulation())]
        except TypeError:
            cmpSel=[convertToAtomSel(cmpSel[0],selection.simulation()),
                    convertToAtomSel(cmpSel[1],selection2.simulation())]
            pass
        if len(cmpSel[0]) != len(cmpSel[1]):
            print("cmpSel: selections are of unequal length:")
            print("  0: %s -- %d atoms" % (cmpSel[0].string(),len(cmpSel[0])))
            print("  1: %s -- %d atoms" % (cmpSel[1].string(),len(cmpSel[1])))
            raise Exception("cmpSel: selections are of unequal length")
        pass
    pot.cmpSel = cmpSel

    

    from ensembleSimulation \
         import EnsembleSimulation_getEnsembleSimulation as \
         getEnsembleSimulation
    
    if (pot.simulation().size()==1 or
        (getEnsembleSimulation(selection.simulation()).size()==1 and
         getEnsembleSimulation(selection2.simulation()).size()==1 )):
        for i in range(len(selection)):
            pot.addEquivAtomSelPair( selection[i], selection2[i] )
            pass
        pass
    else:
        if (pot.simulation() != selection.simulation() or
            pot.simulation() != selection2.simulation() ):
            raise Exception("all selections must be EnsembleSimulations!")
        pot.addEquivAtoms(selection.string(),
                          selection2.string())
        pass

    pot.setThreshold(0.01) #default threshold value

    return pot

def analyze(potList):
    """perform analysis of PosDiffPot terms and return nicely formatted
    summary"""

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'PosDiffPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    
    hasPDB=False
    hasCMPSel=False
    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        if term.pdbFile: hasPDB=True
        if term.cmpSel: hasCMPSel=True
        pass
        

    ret += "%-9s  %6s   %8s" % \
           ( "", "RMSD", "Delta Variance")
    if hasCMPSel: ret += " Comparison RMSD"
    ret += '\n'

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        term.calcEnergy()

        
        print("term:", term.instanceName())
        print(" lowerBound:", term.lowerBound())
        print(" upperBound:", term.upperBound())
        print(" scale:     ", term.scale())
        if term.pdbFile:
            print("comparison PDB:", term.pdbFile)
            pass
        print("selection string:", term.selection.string())
        print(" (Simulation: %s)" % term.selection.simulation().name())
        print("selection2 string:", end=' ') 
        if term.selection.string() == term.selection2.string():
            print("same")
        else:
            print(term.selection2.string())
            pass
        print(" (Simulation: %s)" % term.selection2.simulation().name())
        if term.cmpSel:
            print("Comparison selections")
            print("  0: %s" % term.cmpSel[0].string())
            print("   (Simulation: %s)" % term.cmpSel[0].simulation().name())
            print("  1: %s" % term.cmpSel[1].string())
            print("   (Simulation: %s)" % term.cmpSel[1].simulation().name())
            pass
        if term.planarFit():
            print("planarFit with u=",term.u())

        print("calculated RMSD:", term.rms())

        import mat3
        print("fit rotation:")
        print(term.rotMat())
        print("group centroid positions:")
        print("  0: ", term.centroid(0), "   1: ", term.centroid(1))
 
        ret += "%-9s  %7.3f  %8.3f" % \
               (name , term.rms(), term.deltaVariance())
        if term.cmpSel:

            # this could be faster...
            import mat3
            posA=[]
            for atom in term.cmpSel[0]:
                posA.append( atom.pos() - term.centroid(0) )
                pass
            posB=[]
            for atom in term.cmpSel[1]:
                posB.append( term.rotMat()*(atom.pos() - term.centroid(1)) )
                pass

            from vec3 import norm
            rmsd=0
            for i in range(len(posA)):
                rmsd += norm(posA[i]-posB[i])**2
                pass
            
            from math import sqrt
            if rmsd:
                rmsd /= len(posA)
                rmsd=sqrt(rmsd)
                pass
            
            ret += "      %6.3f" % rmsd
            pass
        ret += '\n'
        pass
    
    return ret


from simulationTools import registerTerm
registerTerm(analyze,"Position Difference Potential","PosDiff",
r"""
For each term the following are reported:
  RMSD            - root mean square deviation between calculated and target
                    values.
  comparison RMSD - root mean square deviation between calculated and target
                    values corresponding to the cmpSel argument, if specified.
  Delta Variance  - Measure of difference in ensemble spread for different
                    equivalent sets. Defined only for <m ensembleSimualtion>
                    calculations. See <m posSymmPot>.
""")
    

        
    
        
