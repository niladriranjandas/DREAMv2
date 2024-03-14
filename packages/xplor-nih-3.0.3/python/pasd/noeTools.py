"""
Output PASD data as resrtaint lists.
"""


def writeXplorAssignments(filename=None,
                          remarksList=[],
                          noePot=None,
                          peaks=[],
                          removeIntersected=True):
    """
    Routine to write XPLOR assignment statements for PASD peaks.

    Given a list of <m pasdPeak>.Peaks or a <m pasdPot>.PASDPot, write
    out appropriate XPLOR-formatted restraints from the associated assignments.
    """

    import pasdPot
    if not noePot and not peaks:
        raise Exception("Neither noePot nor peakist is specified.")

    #make sure its a Python list
    peaks = list( peaks )
    if noePot:
        peaks += list( noePot.peaks() )
        pass

    peaks.sort( key=lambda a: a.name() )

    outfile = open(filename,"w") if filename else sys.stout

    #
    # write remarks
    #
    if len(remarksList):
        outfile.write("!\n")
        for rem in remarksList:
            outfile.write("! %s\n"%rem)
            pass
        outfile.write("!\n")
        pass

    #
    # write each peak
    #
    for peak in peaks:
        writeOneXplorAssignment(outfile,peak,removeIntersected)
        pass
    return

def writeOneXplorAssignment(outfile,
                            peak,
                            removeIntersected):
    """
    Internal routine to be called by writeXplorAssignments.
    """

    outfile.write("! \n")
    outfile.write("! peak name %s\n" % peak.name())
    outfile.write("! bounds %f %f\n" % (peak.upBound(),peak.lowBound()))
    outfile.write("! likelihood %f\n"% peak.prevLikelihood())
    outfile.write("! intensity %f\n" % peak.intensity())
    outfile.write("! from proton shift %f\n" % peak.fromProtonShift())
    outfile.write("! to proton shift %f\n"   % peak.toProtonShift())
    outfile.write("! from heavyatom shift %f\n" % peak.fromHeavyatomShift())
    outfile.write("! to heavyatom shift %f\n"   % peak.toHeavyatomShift())
    
    if peak.note():
        for note in peak.note().split('\n'):
            outfile.write("! note %s\n" % note)
            pass
        pass
    outfile.write("! \n")


    #
    # eliminate peakAssignments that have non-disjoint from and to proton
    # selections, since classic xplor can't handle them
    #
    if removeIntersected:
        numRemoved = removeIntersectedPeakAssignments(peak)
        if numRemoved>0:
            outfile.write(
                "! Removed %d peakAssignments whose from and to proton"
                " selections were not disjoint.\n" %  numRemoved)
            pass
        pass

    pas = list( peak.peakAssignments() )

    if len(pas)==0:
        outfile.write("! NO ASSIGNMENTS FOR THIS PEAK\n")
        outfile.write(" \n")
        return

    pas.sort( key=lambda a: a.name() )

    #
    # run through the peakAssigns, finding the maximal up- and
    # lowBoundCorrections These will be applied to all PAs, since
    # the classic xplor format doesn't allow individual PAs to
    # have individual bound corrections
    #

    maxUBC=max([pa.upBoundCorrection()  for pa in pas])
    maxLBC=max([pa.lowBoundCorrection() for pa in pas])

    upBound  = peak.upBound()  + maxUBC
    lowBound = peak.lowBound() - maxLBC

    middleDist = lowBound + 0.5*(upBound - lowBound)
    dminus     = middleDist - lowBound
    dplus      = upBound - middleDist

    selPairs = [(pa.fromProtonSelection().string() ,
                 pa.toProtonSelection().string()   ) for pa in pas]
    
    comment = peak.name()
    outfile.write("assign (%s) (%s) %f %f %f ! %s\n" %
                  (selPairs[0][0], selPairs[0][1], 
                   middleDist, dminus, dplus,comment))
    
    seen=[selPairs[0]]
    for selPair in selPairs[1:]:
        if selPair in seen:
            continue
        seen.append( selPair )
        outfile.write("    or (%s) (%s) \n" % tuple(selPair))
        pass
    outfile.write("\n")
    return

def removeIntersectedPeakAssignments(peak):
    """
    Removes peakAssignments whose from and to shiftAssignments have 
    intersecting proton selections.
    
    To be used before writing out in classic xplor table format, 
    since the classic xplor NOE potential can't handle non-disjoint selection
    pairs.
    """

    numRemoved=0

    from atomSel import intersection
    for pa in peak.peakAssignments():
	    
        fromSel = pa.fromProtonSelection()
        toSel   = pa.toProtonSelection()

        if intersection( fromSel, toSel ):
            numRemoved += 1
            peak.removePeakAssignmentNamed(pa.name())
            pass
        pass
    return numRemoved

def removeLowLikelihoodPeakAssignments(pot,
                                       cutoff=0.9,
                                       preventUnassignment=False):
    """Remove peak assignments with likelihood lower than the cutoff value.
    If preventUnassignment=True, in the event of a peak having no assignments
    greater than cutoff, keep that with the largest likelihood rather than
    deleting the peak entirely.

    Returns a named tuple with these members
      numRemoved       - number of removed peak assignments
      numPeaksAffected - number of peaks having some peak assignments removed
      numLRremoved     - number of removed long-range peak assignments
      numGoodRemoved   - number of removed good peak assignments
      numGoodLRremoved - number of removed good long-range peak assignments
    """

    numRemoved = 0
    numLRremoved = 0
    numGoodLRremoved = 0     
    numGoodRemoved = 0
    numPeaksAffected = 0    
    for peak in pot.peaks():

        curCutoff = max(cutoff,peak.previousLikelihood()) \
                    if preventUnassignment else cutoff

        peakAffected = False
        for pa in peak.peakAssignments():
            if pa.previousLikelihood() < curCutoff:
                numRemoved += 1
                if pa.isLongRange(6):
                    numLRremoved += 1
                    if pa.isGood(): numGoodLRremoved += 1
                    pass

                if pa.isGood(): numGoodRemoved += 1

                peakAffected = True

                peak.removePeakAssignmentNamed( pa.name() )
                pass
            pass
        if peakAffected: numPeaksAffected += 1
        pass

    # 
    # eliminate any remaining links to PAs just deleted
    #

    pot.clearMissingLinksToPeakAssignments()
    
    from collections import namedtuple
    RetType = namedtuple("removeLowLikelihoodPeakAssignmentsInfo",
                         ('numRemoved',
                          'numPeaksAffected',
                          'numLRremoved',
                          'numGoodRemoved',
                          'numGoodLRremoved'))
    return RetType(numRemoved,
                   numPeaksAffected,
                   numLRremoved,
                   numGoodRemoved,
                   numGoodLRremoved)
    
def makeNEFRestraintLinks(restraintTables,linkageInfo):
    """
    Generate a
    nef_peak_restraint_links saveset linking PASD peak assignments to
    NEF peaks.

    Return a string NEF record.
    """
    from cif import Cif, CifDatablock, CifCategory
    cif=Cif()
    block = CifDatablock()
    cat = CifCategory()
    name="nef_peak_restraint_links"
    cat["sf_category"] = name
    cat["sf_framecode"] = name

    block[name] = cat
    cat = CifCategory()
    for key in "nmr_spectrum_id peak_id restraint_list_id restraint_id".split():
        cat.addKey(key)
        pass

    for restraintTable,linkage in zip(restraintTables,linkageInfo):
        for peakTable,peakID,restraintID in linkage:
            from nefTools import catPrefixes
            cat.addValue("nmr_spectrum_id",
                         catPrefixes["spectrum"]+'_'+peakTable)
            cat.addValue("peak_id",peakID)
            cat.addValue("restraint_list_id",
                         catPrefixes["distance"]+'_'+restraintTable)
            cat.addValue("restraint_id",str(restraintID))
            pass
        pass

    block["nef_peak_restraint_link"] = cat
    block.setIsSaveframe(True)
    cif[name] = block
    cif.setUseTrailingStop(True)
    cif.setUseTrailingSave(True)
    return "\n" + cif.asString() + "\n"
    
    
    
