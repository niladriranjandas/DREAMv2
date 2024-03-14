"""
Functions for creating and analyzing the quartic repulsive
<m repelPot>.RepelPot.
"""

def create_RepelPot(name,
                    selection="not PSEUDO",
                    moveTol=0.5,
                    repel=None,
                    use14=False,
                    selPairs=[],
                    extraRadii=[],
                    suppressException=False,
                    verbose=False):

    """Return an instance of <m repelPot>.RepelPot.

    The selection argument specifies which atoms are included in the
    nonbonded calculation.  moveTol specifies the maximum atom
    displacement before the nonbonded list is recalculated.  repel is
    a scale factor for the atomic radii which is set by the function
    initRepel.  The selPairs argument is an optional sequence of pairs
    of atom selection strings which specify particular pairs of atoms
    to compute interactions between.  The actual selections used will
    be the intersections of the selection argument and each selection
    pair entry.  The use14 argument controls whether or not 1-4
    interactions are enabled.  Setting use14=False is equivalent to
    nbxmod=4 in the <l
    https://nmr.cit.nih.gov/xplor-nih/xplorMan/vdwFunction.html older
    XPLOR VDW term>.  The default setting of use14=False is
    appropriate when using TorsionDBPot created by 
    <m torsionDBPotTools>.create_TorsionDBPot.  The suppressException argument
    should be set to True only if a nonstandard repel value is specified, and
    should only be used for testing purposes.

    The extraRadii argument contains a list of pairs (chemType,radius) which
    can be used to specify atomic radii for chemical types not present in
    the XPLOR PSF.
    """
    from selectTools import convertToAtomSel
    sel = convertToAtomSel(selection)

    sim = sel.simulation()

    if len(selPairs):
        from atomSel import union, AtomSel
        pairs=[]
        for sel1,sel2 in selPairs:
            sel1 = convertToAtomSel(sel1,sim)
            sel2 = convertToAtomSel(sel2,sim)
            pairs.append( ( sel1 , sel2 ) ) 
            pass
        pass
                          

    from nonBondTools import getNBParamsFromXplor
    nbInfo = getNBParamsFromXplor(sim)

    if verbose:
        print("Calculating Repel interactions between %d atoms." % sel.size())


    #
    # now convert these index-based exclusions to per-residue type
    # exclusions where atoms are referenced by name
    #

    nbParams = nbInfo.nbParams
    for key in list(nbParams.keys()):
        params = nbParams[key]
        if params[1]==0. : params[1]=1
        if params[3]==0. : params[3]=1
        nbParams[key] = params
        pass

    nbfixParams = nbInfo.nbfixParams
    for key in list(nbfixParams.keys()):
        params = nbfixParams[key]
        if params[1]==0. : params[1]=1
        if params[3]==0. : params[3]=1
        nbfixParams[key] = params
        pass

    for chemType,radius in extraRadii:
        A=radius**6
        B=1
        # why the factor of 2?
        nbParams[chemType] = (1,2*radius,1,2*radius)
        pass


##=======
    #now do something with nbfixParams: they are  A, B, A14, B14
    #indexed by the string T1*T2, where T1 and T2 are chemical types
    #nonbonded parameters are generated using combination rules from
    # entries given by NONBon specifications, except for those pairs
    # explicitly specified by NBFIx    

##====================

    from atomSel import AtomSel, intersection
    if len(intersection( AtomSel("not known",sel.simulation()), sel)):
        raise Exception("selection includes atoms with uninitialized coordinates")
                               
    from repelPot import RepelPot
    pot = RepelPot(name,nbParams,nbfixParams,nbInfo.resExclude,sel)

    pot.setVerbose(verbose)
    pot.setE14Factor(1)

    pot.savedSelPairs=[]
    initRepel(pot,moveTol=moveTol,repel=repel,
              use14=use14,suppressException=suppressException)

    if selPairs:
        pot.deleteSelectionPair(0)
        for pair in pairs:
            pot.addSelectionPair( *pair )
            pot.savedSelPairs.append( pair )
            pass
        pass

    return pot

def checkGetRepel(repel=None,
                  suppressException=False):
    """
    Given a repel value, check that it is consistent with the loaded Xplor-NIH
    parameter set, and throw an exception if it is inconsistent, unless the
    suppressException argument is set to True. A check is also made that
    the parameter sets for proteins and nucleic acids take the same repel
    value.

    If the default Xplor-NIH parameter set is not loaded, do not modify repel.
    """
    
    cutVersion=3.1 # parameter versions greater or equal to this use
                   # molprobity radii and require repel values of 0.9
                   # while earlier versions can have smaller values (default 
                   # is 0.8).
    vers=0.
    prevVers=None
    from protocol import parVersion
    foundVersion=False
    for system in ('nucleic','protein'):
        if system in parVersion:
            foundVersion=True
            try:
                vers = float(parVersion[system])
            except ValueError:
                #version is more complicated than simple floating point number
                # so ignore it.
                break
            pass
        if repel and vers>=cutVersion and repel<0.9:
            message= "Repel value of %f is inconsistent with " % repel + \
                     "%s parameter version %s" % (system,vers)
            if suppressException:
                print("initRepel: Warning! %s" %message)
            else:
                raise Exception(message)
            pass
        if prevVers and (prevVers<cutVersion  and vers>=cutVersion or
                         prevVers>=cutVersion and vers<cutVersion  ):
            message= "Mixing protein and nucleic parameter sets with " + \
                     "inconsistent atomic radii. " + \
                     "[Nucleic version %.1f, Protein version %.1f]" % (prevVers,
                                                                       vers)
            if suppressException:
                print("initRepel: Warning! %s" %message)
            else:
                raise Exception(message)
            pass
        prevVers=vers
        pass

    #non-standard parameters- use default repel value
    if not foundVersion: return repel

    if not repel: repel = 0.8 if vers<cutVersion else 0.9

    return repel
               

def initRepel(term,
              scale=1,
              use14=None,
              moveTol=0.5,
              repel=None,
              interactingAtoms=None,
              suppressException=False,
              ):
    """
    A helper routine for resetting a RepelPot terms to default values.
    FIX: rationalize default values vs. None
    Additionally, a simple selection pair can be specified by specifying a
    single atom selection string as the interactingAtoms argument. If the
    interactingAtoms argument is omitted, the selection pairs will be reset to
    those specified in create_RepelPot.

    The default value of repel is 0.8 for pre-3.1 version Xplor-NIH protein
    and nucleic acid parameters, and 0.9 for later versions. If a repel value
    less than 0.9 is specified to be used with post-3.0 parameter sets, an
    exception will be raised, unless suppressExceptions is set to True.

    """

    repel = checkGetRepel(repel,suppressException)
    
    if "AvePot" in str(type(term)): term = term.subPot()
    term.setScale(scale)
    term.setMoveTol(moveTol)
    if repel!=None: term.setRepel(repel)
    if use14!=None: term.setUse14( use14 )
    sim = term.simulation()
    from selectTools import convertToAtomSel
    term.resetSelectionPairs()
    if interactingAtoms:
        term.deleteSelectionPair(0)
        term.addSelectionPair(  convertToAtomSel(interactingAtoms,sim),
                                convertToAtomSel(interactingAtoms,sim) )
    else:
        if term.savedSelPairs:
            term.deleteSelectionPair(0)
            for pair in term.savedSelPairs:
                term.addSelectionPair( *pair )
                pass
            pass
        pass
    return

def analyze(potList):
    """perform analysis of RepelPot terms and return nicely formatted summary.

    """

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'RepelPot')

    if not potList: return ret

    potList.sort(key=lambda x: x.instanceName())

    for pot in potList:
        name = pot.instanceName()

        print()
        print("  Nonbonded RepelPot Violations")
        print("        threshold: %f" % pot.threshold())
        print("        scale: %.2f"   % pot.scale())
        print()
        print("%26s %26s %6s  %7s" % ("    ", "     ",
                                      "actual", "allowed"))
        print("%26s %26s %6s  %6s  %6s %6s" % ("atom1     ", "atom2     ",
                                               " dist", " dist ", "diff  ",
                                               "energy"))
        print("_"*78)

        selPairs = pot.selectionPairs()
        from math import log10, floor
        for ipair, selPair in enumerate( selPairs ):
            aString = selPair.a.string()
            aSize = selPair.a.size()
            bString = selPair.b.string()
            bSize = selPair.b.size()
            strLen = max(len(aString),len(bString))
            numLen = int(floor(log10(max(aSize,bSize,1))))
            print("""  selection pair %2d: (%-*s) [ %*d atoms ]
                     (%-*s) [ %*d atoms ]""" % (ipair,
                                                strLen,aString,
                                                numLen,aSize,
                                                strLen,bString,
                                                numLen,bSize))
            bumps = pot.bumps(ipair)

            for b in bumps:
                print("%26s %26s %6.2f %6.2f %6.2f %6.3f" % (b.atomi.string(),
                                                             b.atomj.string(),
                                                             b.dist, b.dist0,
                                                             b.dist0-b.dist,
                                                             b.energy))
                pass
            pass

        print()
        print(pot.info())

        pass
    
    return ret


from simulationTools import registerTerm
registerTerm(analyze,"RepelPot terms","RepelPot",
r"""
No header information is reported for this term.
""")
