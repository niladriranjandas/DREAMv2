
def create_PASDPot(name,
                   shifts=None,
                   peaks=None,
                   exceptions=None,
                   asTuple=False,
                   verbose=True):
    """
    Create a <m pasdPot>.PASDPot object with specified instanceName.
    Optionally. shift assignments, peak definitions and inverse exceptions
    can be specified in the shifts, peaks and exceptions arguments,
    respectively as either a string or a filename.

    If asTuple=True, the return value is a tuple of forward and inverse PASD
    energy terms, otherwise the combination is returned.
    """

    from tclInterp import TCLInterp
    tcl = TCLInterp()
    tcl.command("package require marvin")

    from pasdPot import NOEPot    #create PASD energy term
    ret=NOEPot(name)
    ret.setUseSingleAssignment(True)

    #make TCL version of PASD NOEPot
    from pyInterp import portableStringRep
    tcl.command("rc_PASDPot noe_%s -this $ptr\n" % name,
                ("ptr",portableStringRep(ret)))
    tcl.command("puts $errorInfo")

    import os, sys
    from simulationTools import mktemp
    if shifts:
        if os.path.exists(shifts):
            shiftsFilename = shifts
            shifts = open(shifts).read()
        else:
            shiftsFilename = "string"
            pass
        tmpfilename = mktemp('.shiftAssignments')
        open(tmpfilename,'w').write(shifts)
        if verbose:
            print('reading shift assignments from %s ...' % shiftsFilename, end=' ')
            sys.stdout.flush()
            pass
        tcl.command('readShiftAssignments -fileName "%s" -pot noe_%s' %
                    (tmpfilename,name))
        os.unlink(tmpfilename)
        if verbose:
            print(' [%d loaded] ' % len(ret.shiftAssignments()))
            pass
        pass

    if peaks:
        if os.path.exists(peaks):
            peaksFilename=peaks
            peaks = open(peaks).read()
        else:
            peaksFilename="string"
            pass
        tmpfilename = mktemp('.peaks')
        open(tmpfilename,'w').write(peaks)
        if verbose:
            print('reading peaks from %s ...' % peaksFilename, end=' ')
            sys.stdout.flush()
            pass
        tcl.command('readMarvinPeaks -fileName "%s" -pot noe_%s' %
                    (tmpfilename,name))
        os.unlink(tmpfilename)
        if verbose:
            print(' [%d loaded] ' % len(ret.peaks()))
            pass
        pass

    if exceptions:
        if os.path.exists(exceptions):
            exceptionsFilename = exceptions
            exceptions = open(exceptions).read()
        else:
            exceptionsFilename = "string"
            pass
        tmpfilename = mktemp('.exceptions')
        open(tmpfilename,'w').write(exceptions)
        if verbose:
            print('reading exceptions from ' + exceptionsFilename, end=' ')
            sys.stdout.flush()
            pass
        tcl.command('readExplicitInverseExceptions -fileName "%s" -pot %s' %
                    (tmpfilename,"noe_"+name))
        os.unlink(tmpfilename)
        if verbose:
            print(' [%d loaded] ' % len(ret.exceptions()))
            pass
        pass

    
    if asTuple:
        retVal = (ret.distPot(),ret.invPot())
        retVal[0].base = ret
        ret=retVal
        pass
    
    return ret


def linearInterp(frac,minMax):
    """
    Return interpolated value between start and stop at fraction step/steps

    step takes values between 0 and steps-1
    """
    return minMax[0] + frac * (minMax[1] - minMax[0])
    pass

def maybeShuffleAssignments(pasdTerms,
                            fracDone,
                            reshuffleProb,
                            potProperties,
                            numMCsteps=1,
                            charDeltaFac=1,
                            verbose=False,
                            ):
    import random
    if random.uniform(0,1)>= reshuffleProb:
        return

    potPropertiesVals=[]
    for name,range in potProperties:
        potPropertiesVals.append((name,linearInterp(fracDone,range)))
        pass                                 
    
    for pot in pasdTerms:
        for name,val in potPropertiesVals:
            setName = name[0].upper() + name[1:]
            eval("pot.set%s(%f)"%(setName,val))
#            setVal=eval("pot.%s()"%name)
#            print pot.instanceName(), name, val, setVal, fracDone
            pass

        from pasdPotTools import updateActivation
        updateActivation(pot,
                         numMCsteps=numMCsteps,
                         charDeltaFac=charDeltaFac,
                         verbose=verbose)
        pass
    return
    
def updateActivation(pot,
                     numMCsteps=1,
                     charDeltaFac=1.,
                     verbose=False):
    """
    Update which assignments are active using a Monte Carlo algorithm.
    """

    charDeltaDV0 = pot.characteristicDeltaDV()
    charDeltaPL0 = pot.characteristicDeltaPL()
    charDeltaNC0 = pot.characteristicDeltaNoeCompleteness()
    charDeltaPS0 = pot.characteristicDeltaScatter()
    
    cdf = charDeltaFac

    if verbose:
        print("Updating NOE activation in PASDPot %s by Monte Carlo\n" % \
              pot.instanceName())
        pass

    for i in range(numMCsteps):
	
        if verbose:
            print("MC step %d of %d:" % (i, numMCsteps))
            pass

        fracDone = float(i) / numMCsteps            

        pot.setCharacteristicDeltaDV( 
            linearInterp( fracDone, (charDeltaDV0,cdf*charDeltaDV0)) )
        pot.setCharacteristicDeltaPL(
            linearInterp( fracDone, (charDeltaPL0,cdf*charDeltaPL0)) )
        pot.setCharacteristicDeltaNoeCompleteness(
            linearInterp( fracDone, (charDeltaNC0,cdf*charDeltaNC0)) )
        pot.setCharacteristicDeltaScatter(
	    linearInterp( fracDone, (charDeltaPS0,cdf*charDeltaPS0)) )

        #
        # updateActivation returns a list of numbers that describe the 
        # move it just took. It throws an error if the number of 
        # unaccepted MC tries is too large, which usually means we're
        # at a local minimum
        #
        
        try:
            summary = pot.updateActivation()
        except:
            if verbose:
                print("Gave up after too many unaccepted MC tries. ", end=' ')
                print("Reverting to previous step's activations")
                pass
            break
        
        if verbose:
            print("   violation score           %.5f --> %.5f" \
                  % tuple(summary[:2]))
            print("   previous likelihood score %.5f --> %.5f" \
                  % tuple(summary[2:4]))
            print("   NOE completeness score    %.5f --> %.5f" \
                  % tuple(summary[4:6]))
            print("   peak posn scatter score   %.5f --> %.5f" \
                  % tuple(summary[6:8]))
            print("   normal energy             %.5f --> %.5f" \
                  % tuple(summary[8:10]))
            print("   inverse energy            %.5f --> %.5f" \
                  % tuple(summary[10:12]))
            print("   number of inactive peaks  %5d --> %5d" \
                  % tuple([int(e) for e in summary[12:14]]))
            print("   number of shift assignments in wrong neighborhoods %5d --> %5d" \
                  % tuple([int(e) for e in summary[14:16]]))
            print("   fraction of peak assignments inactive %.5f --> %.5f" \
                  % tuple(summary[16:18]))
            if summary[18] > 0 or summary[19] > 0:
                print("   frac active longrange forces that are good  %.5f --> %.5f" \
                      % tuple(summary[18:20]))
                print("   frac good longrange peak assignments active ", end=' ')
                print("%.5f --> %.5f" % tuple(summary[20:22]))
                pass
            
            
            print("   Took %d tries, overall likelihood = %.5f" \
                  % (int(summary[22]), summary[23]))
            print("\n")
            pass
        pass

    pot.setCharacteristicDeltaDV(              charDeltaDV0 )
    pot.setCharacteristicDeltaPL(              charDeltaPL0 )
    pot.setCharacteristicDeltaNoeCompleteness( charDeltaNC0 )
    pot.setCharacteristicDeltaScatter(         charDeltaPS0 )
    return

def highTemp1_updateTerms(terms       ,
                          fractionDone,
                          totSteps    ,
                          numNOEshuffles=64,
                          numMCsteps=10,
                          charDistViol=10.0,
                          charDeltaDV=0.01,
                          prevLikeliWeight=1
                          ):
    

    maybeShuffleAssignments(terms ,
                            fractionDone,
                            reshuffleProb=float(numNOEshuffles)/totSteps,
                            charDeltaFac=1,
                            numMCsteps=1,
                            verbose=True,
                        potProperties=[ #linear ramping for these
    ("characteristicNoeCompleteness"     , ( 0.0  , 0.0  )),
    ("scale"                             , ( 1    , 1    )),
    ("characteristicScatter"             , (9999999.9    , 9999999.9    )),
    ("characteristicDistanceViol"        , (charDistViol , charDistViol )),
    ("distanceViolationWeight"           , ( 0    , 0    )),
    ("previousLikelihoodWeight"          , ( 1    , 1    )),
    ("noeCompletenessWeight"             , ( 0    , 0    )),
    ("scatterWeight"                     , ( 0    , 0    )),
    ("characteristicDeltaDV"             , ( 0.01 , 0.01 )),
    ("characteristicDeltaPL"             , ( 0.1  , 0.1  )),
    ("characteristicDeltaNoeCompleteness", ( 0.05 , 0.05 )),
    ("characteristicDeltaScatter"        , ( 0.001, 0.001)),
    ]
                            )
    
    return

def highTemp2_updateTerms(terms       ,
                          fractionDone,
                          totSteps    ,
                          numNOEshuffles=64,
                          numMCsteps=10,
                          charDistViol=10.0,
                          charDeltaDV=0.1,
                          prevLikeliWeight=1
                          ):
    

    maybeShuffleAssignments(terms ,
                            fractionDone,
                            reshuffleProb=float(numNOEshuffles)/totSteps,
                            charDeltaFac=0.25,
                            numMCsteps=10,
                            verbose=True,
                        potProperties=[ #linear ramping for these
    ("characteristicNoeCompleteness"     , ( 0.05 , 0.05 )),
    ("scale"                             , ( 1    , 1    )),
    ("characteristicScatter"             , ( 1    , 1    )),
    ("characteristicDistanceViol"        , (10.0  , 10.0 )),
    ("distanceViolationWeight"           , ( 1    , 1    )),
    ("previousLikelihoodWeight"          , ( 1    , 1    )),
    ("noeCompletenessWeight"             , ( 1    , 1    )),
    ("scatterWeight"                     , ( 0    , 0    )),
    ("characteristicDeltaDV"             , ( 0.1  , 0.1  )),
    ("characteristicDeltaPL"             , ( 0.1  , 0.1  )),
    ("characteristicDeltaNoeCompleteness", ( 0.05 , 0.05 )),
    ("characteristicDeltaScatter"        , ( 0.001, 0.001)),
    ]
                            )
    
    return

def cooling_updateTerms(terms       ,
                        fractionDone,
                        totSteps    ,
                        numNOEshuffles=64,
                        numMCsteps=10,
                        distViolWeight=(1,1),
                        prevLikeliWeight=(1,0),
                        noeCompleteWeight=(1,1),
                        charDistViol=(10.0,2.0),
                        charNoeComplete=(0.05,0.05),
                        charScatter=(1,1),
                        charDeltaDV=(0.01,0.01),
                        charDeltaNoeComplete=(0.05,0.05),
                        charDeltaScatter=( 0.001, 0.001),
                        ):
    

    maybeShuffleAssignments(terms ,
                            fractionDone,
                            reshuffleProb=float(numNOEshuffles)/totSteps,
                            charDeltaFac=0.25,
                            numMCsteps=10,
                            verbose=True,
                            potProperties=[ #linear ramping for these
        ("characteristicNoeCompleteness"     , charNoeComplete),
        ("characteristicScatter"             , charScatter),
        ("characteristicDistanceViol"        , charDistViol   ),
        ("distanceViolationWeight"           , distViolWeight ),
        ("previousLikelihoodWeight"          , prevLikeliWeight),
        ("noeCompletenessWeight"             , noeCompleteWeight),
        ("scatterWeight"                     , ( 0    , 0    )),
        ("characteristicDeltaDV"             , charDeltaDV    ),
        ("characteristicDeltaPL"             , ( 0.1  , 0.1  )),
        ("characteristicDeltaNoeCompleteness", charDeltaNoeComplete),
        ("characteristicDeltaScatter"        , charDeltaScatter)
        ]

                            )
    
    return

def writeNEF(noe,name,
             removeIntersected=True):
    """Given an <m pasdPot>.NOEPot term,
    write a NEF record with the specified saveframe name.
    If removeIntersected==True eliminate peakAssignments that have
    non-disjoint from and to proton selections.

    Return the a tuple of the string associated with the generated NEF object,
    and a list of peak linkage information.
    """


    from nefTools import distanceRestraintHeader
    commentStr,cif = distanceRestraintHeader(name,
                                  "\nRestraints from PASDPot term: %s\n"%
                                  noe.instanceName(),
                                             origin="NOE")
    
    
    cat = cif[name]["nef_distance_restraint"]
    id=0
    index=0
    peaks = list(noe.peaks())
    peaks.sort( key=lambda a: a.name() )
    linkageInfo=[]
    for peak in peaks:
        id,index,str = writeOneNEFAssignment(cat,peak,
                                             removeIntersected,id,index)
        commentStr += str
        if len( peak.peakAssignments() ) == 0:
            continue
        notes=[note for note in peak.note().split('\n')
               if note.startswith('from NEF table')]
        foundLinkage=False
        for note in notes:
            noteFields=note.split()
            if len(noteFields)!=6:
                continue
            peakTable=noteFields[3]
            peakID = noteFields[5]
            linkageInfo.append( (peakTable,peakID,id) )
            foundLinkage=True
            break
        if not foundLinkage:
            print("writeNEF: Warning: no NEF peak "
                  "specified in peak with name",peak.name())
            print("  found notes:",notes)
            pass

        pass

    return (commentStr + cif.asString(), linkageInfo)

def writeOneNEFAssignment(cat,
                          peak,
                          removeIntersected,
                          id,index):
    """
    Internal routine to be called by writeNEFAssignments. Returns the comments
    associated with the particular peak
    """
    
    from nefTools import toNefAtomname, nefComment
    
    commentStr = ''
    for comment in ["",
                    "peak name %s" % peak.name(),
                    "bounds %f %f" % (peak.upBound(),peak.lowBound()),
                    "likelihood %f"% peak.prevLikelihood(),
                    "intensity %f" % peak.intensity(),
                    "from proton shift %f" % peak.fromProtonShift(),
                    "to proton shift %f"   % peak.toProtonShift(),
                    "from heavyatom shift %f" % peak.fromHeavyatomShift(),
                    "to heavyatom shift %f"   % peak.toHeavyatomShift(),]:
        commentStr += nefComment(comment)
        pass
    
    if peak.note():
        commentStr += nefComment( peak.note() )
        pass
    commentStr += nefComment( "" )


    #
    # eliminate peakAssignments that have non-disjoint from and to proton
    # selections, since classic xplor can't handle them
    #
    if removeIntersected:
        from pasd.noeTools import removeIntersectedPeakAssignments
        numRemoved = removeIntersectedPeakAssignments(peak)
        if numRemoved>0:
            commentStr += nefComment("Removed %d peakAssignments whose from "
                                     "and to proton selections were not "
                                     "disjoint.\n" %  numRemoved)
            pass
        pass
    
    pas = list( peak.peakAssignments() )

    if len(pas)==0:
        commentStr += nefComment("NO ASSIGNMENTS FOR THIS PEAK\n")
        return (id,index,commentStr)

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

    selPairs = [(pa.fromProtonSelection() ,
                 pa.toProtonSelection()   ) for pa in pas]

    from nefTools import addOneDistanceRestraint
    id,index= addOneDistanceRestraint(cat,selPairs,index,id,
                                      middleDist,lowBound,upBound,# FIX: weight?
                                      XplorNIH_label=peak.name())
    return (id,index,commentStr)


                                   
