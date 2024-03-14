
def initLikelihoods(pot,
                    maxFiltersFailed=0):
    """set previous likelihoods to 0 or 1 based on numFilteredFailed relative
    to maxFiltersFailed."""
    for peak in pot.peaks():
        for pa in peak.peakAssignments():
            pa.setPreviousLikelihood( 0 if (pa.numFiltersFailed()>
                                            maxFiltersFailed    )
                                      else 1 )
            pass
        pass
    return
    

def jointFilter(pots,
                minPeakScore=-1,
                minExpectedScore=0.08,
                primarySeqFilter=False,
                refPDB=None,
                writeFiles=False,
                assignmentThreshold=None,
                activePAThreshold=None,
                inactiveAssignmentThreshold=None,
                deleteNonIntraPAs=False,
                filenamePrefix="NAME_pass2",
                **kwargs
                ):
    """
    Filenames are generated from filenamePrefix by substituting the PASDPot
    instanceName for the 'NAME' literal, and then appending .exceptions,
    .peaks or .shiftAssignments for the respective filenames.

    assignmentThreshold - if set, delete all peaks which have more peak
                          assignments than this value.

    inactiveAssignmentThreshold - if set, delete peaks with no active peaks,
                                  but with more assignments than this value.
                                  
    activePAThreshold   - if set, delete peaks with more active peak
                          assignments than this value.
    
    deleteNonIntraPAs   - if True, delete all non-intramolecular peak
                          assignments if there is an intramolecular peak
                          assignment.

    These additional arguments are passed through to netfilter.netFilter:
       knownContacts, printResiduePairScores, passFrac, numIters,
       initScoresFrom

    The above filters are run after the network filter.
    """

    netFilterKWArgs={}
    for argname in list(kwargs.keys()):
        if argname in """knownContacts printResiduePairScores
                       passFrac numIters initScoresFrom""".split():
            netFilterKWArgs[argname] = kwargs.pop(argname)
            pass
        pass
    if len(kwargs):
        raise Exception("unexpected keyword argument(s): " + str(kwargs))
    from pasd.netfilter import netFilter
    nfResults = netFilter(pots,verbose=2,
                          minPeakScore=minPeakScore,
                          minExpectedScore=minExpectedScore,
                          **netFilterKWArgs)
    
    print((nfResults.remarks))
    maxFiltersFailed=0 #max allow to fail
    minLikelihood=0.9
    maxLikelihood=2.0

    from tclInterp import TCLInterp
    tcl = TCLInterp()
    
    for pot in pots:
        name = pot.instanceName()
        prefix = filenamePrefix.replace("NAME",name)
        exceptFilename= prefix + ".exceptions"
        peakFilename  = prefix + ".peaks"
        saFilename    = prefix + ".shiftAssignments"
    
    
        saRemarks=[]
        peakRemarks=nfResults.remarks
    
        tclRemarks = "{"+peakRemarks+"} "
        tcl.command('set peakRemarks {%s}' % tclRemarks)

            
        if not refPDB: refPDB = ""
        tcl.command("""
     	initialShiftAssignmentAnalysis \
     	    -pot noe_%s \
     	    -referenceStructureFile "%s" \
     	    -minLikelihood 0.9 \
     	    -description "after joint network filter" \
     	    -remarksVariableName saRemarks""" % (name,refPDB))
        tcl.command("""
     	initialPeakAnalysis \
     	    -pot noe_%s \
     	    -violCutoff 0.5 \
     	    -minLikelihood 0.9 \
     	    -maxLikelihood 2.0 \
     	    -referenceStructureFile "%s" \
     	    -description "after joint network filter" \
     	    -remarksVariableName peakRemarks""" % (name,refPDB))
    
        tcl.command("""
    	addPrimarySequenceExceptions \
    	    -pot noe_%s \
    	    -remarksVariableName peakRemarks""" % name)
    
        if primarySeqFilter:
            tcl.command("""
    	    primarySequenceDistanceFilter \
    		-pot noe_%s \
    		-intraresidue \
    		-remarksVariableName peakRemarks""" % name)
    
            
            initLikelihoods(pot)
            pass
    
        #activeAssignmentThreshold=4
        removedPeaks=0
        removedPAs=0
        for peak in pot.peaks():
            #print 'peak:', peak.name(), peak.numPeakAssignments()
            if (assignmentThreshold!= None and
                peak.numPeakAssignments()>assignmentThreshold):
                pot.removePeakNamed(peak.name())
                removedPeaks+=1
                continue
            if (activePAThreshold!=None and
                peak.numActivePeakAssigns()>activePAThreshold):
                pot.removePeakNamed(peak.name())
                removedPeaks+=1
                continue
            if (inactiveAssignmentThreshold!=None and
                peak.numActivePeakAssigns()==0 and
                peak.numPeakAssignments()>inactiveAssignmentThreshold):
                pot.removePeakNamed(peak.name())
                removedPeaks+=1
                continue

            #
            # Delete all non-intraresidue assignments from intraresidue
            # peaks
            #
            if (deleteNonIntraPAs and
                peak.isIntraresidue()):
                for pa in peak.peakAssignments():
                    if not pa.isIntraresidue():
                        peak.removePeakAssignmentNamed( pa.name() )
                        removedPAs += 1
                        pass
                    pass
                pass
            pass
        
        print(("removed %d peaks and %d peak assignments" % (removedPeaks,
                                                            removedPAs)))
        if writeFiles:
            tcl.command("""
            writePASDFiles -pot noe_%s \
            -referenceStructureFile "%s" \
            -exceptionFileName %s \
            -peakFileName %s \
            -shiftAssignmentsFileName %s \
            -minLikelihood %f \
            -maxLikelihood %f \
            -peakRemarks $peakRemarks \
            -shiftAssignmentRemarks $saRemarks""" % (name,
                                                     refPDB,
                                                     exceptFilename,
                                                     peakFilename,
                                                     saFilename,
                                                     minLikelihood,
                                                     maxLikelihood))
            pass
    
        pass
    
def markGoodPeakAssigns(peaks,
                        cutoff=0.5):
    """
    set as good all peak assignments which satisfy distance restraints
    to within cutoff.
    """

    #
    # set the official "good/bad" flag for each peakAssignment
    #

    nBadPAs          = 0
    nGoodPAs         = 0
    nBadPeaks        = 0
    nGoodPeaks       = 0
    nUnassignedPeaks = 0

    for peak in peaks:
        peakGood= False
        
        for pa in peak.peakAssignments():
            v = pa.violation(peak.upBound(), peak.lowBound())
            
            if v != -1 and v < cutoff:
                pa.setGood()
                peakGood=True
                nGoodPAs += 1
            else:
                pa.setBad()
                nBadPAs += 1
                pass
            pass
        
        if not peak.isAssigned():
            nUnassignedPeaks += 1
        elif peakGood:
            nGoodPeaks += 1
        else:
            nBadPeaks += 1
            pass
        pass
    
    
    remarks = "MarkGoodPeakAssigns:  "
    remarks += "using distance violation cutoff of %f A\n" % cutoff
    remarks += "  %d peakAssignments were marked good," % nGoodPAs
    remarks += " %d were marked bad" % nBadPAs
    remarks += "  %d peaks had a good peakAssignment," % nGoodPeaks
    remarks += " %d peaks did not" % nBadPeaks
    remarks += "  %d peaks had no peakAssignments at all" % nUnassignedPeaks

    return remarks

    

def writePASDFiles(pot,
                   shifts,
                   peaks,
                   exceptions,
                   refPDB="",
                   minLikelihood=0.8,
                   maxLikelihood=1.0,
                   ):
    """
    pot is a PASDPot term
    """
    from tclInterp import TCLInterp
    tcl = TCLInterp()

    #make TCL version of PASD NOEPot
    from pyInterp import portableStringRep
    from pyInterp import portableStringRep
    tcl.command("rc_PASDPot noe_%s -this $ptr\n" % pot.instanceName(),
                ("ptr",portableStringRep(pot)))
    
    tcl.command('''
    writePASDFiles -pot noe_%s \
    -referenceStructureFile "%s" \
    -shiftAssignmentsFileName %s \
    -peakFileName %s \
    -exceptionFileName %s \
    -minLikelihood %f \
    -maxLikelihood %f  '''%(pot.instanceName(),
                            refPDB,
                            shifts, peaks, exceptions,
                            0.8,
                            1.0))

import trace
#@trace.notrace_decorate
def activateLikelyPAs(pot,
                      threshold=0.9):
    
    for peak in pot.peaks():
        for pa in peak.peakAssignments():
            if pa.previousLikelihood() >= threshold: pa.activate();
            pass
        pass
    return
                      
