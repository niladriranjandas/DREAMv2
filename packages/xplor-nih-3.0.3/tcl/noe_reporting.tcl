#
# procedures for reporting NOE agreements, etc
#
# JJK 4/24/02
#

package provide marvin 1.0
package require xplorpot
package require marvin 

proc xplorReport args {

    # expand this to include testing and reporting of specific
    # violations    

    set temp [create_XplorPot BOND] 
    set bondEner [$temp calcEnergy]
    rename $temp ""

    set temp [create_XplorPot ANGL] 
    set angleEner [$temp calcEnergy]
    rename $temp ""

    set temp [create_XplorPot IMPR] 
    set imprEner [$temp calcEnergy]
    rename $temp ""

    set temp [create_XplorPot VDW] 
    set vdwEner [$temp calcEnergy]
    rename $temp ""

    set temp [create_XplorPot RAMA] 
    set ramaEner [$temp calcEnergy]
    rename $temp ""

    set temp [create_XplorPot CDIH]
    set cdihEner [$temp calcEnergy]
    rename $temp ""

    return [format \
    "xplor energies:  bond: %f angle: %f impr: %f vdw: %f rama: %f cdih: %f" \
     $bondEner $angleEner $imprEner $vdwEner $ramaEner $cdihEner]

}

proc compareInstanceNames {a b} {
    return [string compare [$a instanceName] [$b instanceName]]
}

proc energyReport pots {

    # sort by instance name
    set pots [lsort -command compareInstanceNames $pots]

    set tot 0
    foreach pot $pots {
	set tot [expr $tot + [$pot calcEnergy]]
    }

    set ret [list [format "%10s %10s %10s %10s" Type Name Energy RMSD]]
    lappend ret [format "%21s %10.3f" "Total Energy" $tot]
    foreach pot $pots {
	lappend ret [format "%10s %10s %10.3f %10.3f" \
			[$pot potName] \
			[$pot instanceName] \
			[$pot calcEnergy]   \
			[$pot rms]        ]
    }

    return $ret
}


proc individualStructNoeReports {pots} {

    set ret [list [format "%10s %8s %8s %8s %8s %8s" \
		       "Marvin NOE" Num Num    "Num Non-" Num RMSD]]
    lappend ret [format "%10s %8s %8s %8s %8s %8s" \
		     Potential Peaks Violated Violated Active ""]

    foreach pot $pots {
	set peaks [$pot peaks]

	set rms [$pot rms]

	lappend ret [format "%10s %8d %8d %8d %8d %8.3f" \
			 [$pot instanceName] \
			 [$pot numPeaks] \
			 [numViolatedPeaks  -peakList $peaks] \
			 [numNonviolatedPeaks -peakList $peaks] \
			 [numActivePeaks -peakList $peaks] \
			 $rms]
    }
    return $ret
}





#
# Given a list of <filename, atomPosArray> pairs from grabPDBfiles, 
# check each one's violations of a given set of peaks
# and return a list of the frac% best <filename, atomPosArray> pairs
#
# Could be modified to do other tests
#

proc chooseBestFraction args {
    
    set files              [requiredFlagVal $args -files]
    set peakList           [requiredFlagVal $args -peakList]
    set completenessPots   [flagVal $args -completenessPots]
    set frac               [flagVal $args -frac 0.10]
    set violCutoff         [flagVal $args -violCutoff 0.5]
    set invBound           [requiredFlagVal $args -inverseBound]
    set invMeth            [requiredFlagVal $args -inverseMethylCorrection]
    set completenessWeight [requiredFlagVal $args -completenessWeight]
    set useNonviolCount    [flagExists $args -useNonviolatedCount]
    set numTalos           [flagVal $args -numTalosRestraints 0]

    foreach pot $completenessPots {
	$pot setInverseBound $invBound
	$pot setInverseMethylCorrection $invMeth
    }


    set nViols    [list]
    set NCs       [list]
    set violRepts [list]

    foreach file $files {
	
	xplorSim setAtomPosArr [lindex $file 1]
	
	if {$useNonviolCount} {
	    set curNViols [numNonviolatedPeaks \
			       -peakList $peakList -violCutoff $violCutoff]
	} else {
	    set curNViols [numViolatedPeaks \
			       -peakList $peakList -violCutoff $violCutoff]
	}

	if {$numTalos > 0} {
	    XplorCommand "print threshold 0.01 cdih"
	    set curTalosNonViol \
		[expr $numTalos - [XplorVariableNamed violations]]
	} else {
	    set curTalosNonViol 0
	}

	set curNViols [expr $curNViols + $curTalosNonViol]

	set totClose 0
	set totCloseAcc 0

	foreach pot $completenessPots {
	    $pot updateActivation
	
	    set totClose \
		[expr $totClose    + [$pot numCloseShiftAssignmentPairs]]
	    set totCloseAcc \
		[expr $totCloseAcc + \
		     [$pot numCloseAccountedShiftAssignmentPairs]]
	}

	set curNC [expr $totCloseAcc / double($totClose)]

	lappend violRepts [list $file $curNViols $curNC]
	lappend nViols $curNViols
	lappend NCs $curNC
    }
    
    set minNViols [min $nViols]
    set maxNViols [max $nViols]

    set nViolRange [expr $maxNViols - $minNViols]
    if {$nViolRange == 0} {
	set nViolRange 1
    }

    set minNC     [min $NCs]
    set maxNC     [max $NCs]

    set ncRange [expr $maxNC - $minNC]
    if {$ncRange == 0} {
	set ncRange 1
    }

    set completeViolRepts [list]

    foreach violRept $violRepts {
	set curFile   [lindex $violRept 0]
	set curNViols [lindex $violRept 1]
	set curNC     [lindex $violRept 2]

	if {$useNonviolCount} {
	    set curViolScore \
		[expr 1- (($curNViols - $minNViols) / double($nViolRange))]
	} else {
	    set curViolScore \
		[expr ($curNViols - $minNViols) / double($nViolRange)]
	}

	set curNCScore [expr 1 - (($curNC - $minNC) / double($ncRange))]
	set curScore \
	    [expr ((1 - $completenessWeight) * $curViolScore) + \
		 ($completenessWeight * $curNCScore)]
	
	lappend completeViolRepts \
	    [list $curFile $curNViols $curNC \
		 $curViolScore $curNCScore $curScore]
    }

    set sortedViolRepts [lsort -real -index end $completeViolRepts]

    foreach vr $sortedViolRepts {
	puts $vr
    }

    set numToReturn    [expr ceil($frac * [llength $sortedViolRepts])]
    
    if {$numToReturn > [llength $sortedViolRepts]} {
	set numToReturn [llength $sortedViolRepts]
    }
    
    if {$numToReturn < 1} {
	set numToReturn 1
    }
    
    set retVal [list]
    
    for {set x 0} {$x < $numToReturn} {incr x} {
	set curViolRept [lindex $sortedViolRepts $x]
	lappend retVal [lindex $curViolRept 0]
    }
    
    return $retVal
}

#
# Returns a histogram of number of violations of a given set of peaks
# in a given set of <fname, atomPosArray> pairs.  Can optionally give 
# a detailed report, listing filenames and numbers of violations.
# For summarizing violations in all structs or converged structs
#

proc reportEnsemblePeakViolations args {

    set files      [requiredFlagVal $args -files]
    set peakList   [requiredFlagVal $args -peakList]
    set violCutoff [flagVal $args -violCutoff 1.0]
    set histoTitle [requiredFlagVal $args -histogramTitle]
    set detailed   [flagExists $args -detailed]

    set violList [list]
    foreach file $files {
	
	xplorSim setAtomPosArr [lindex $file 1]
	
	set curNViols [numViolatedPeaks -peakList $peakList -violCutoff $violCutoff]
	lappend violList $curNViols
    }

    set retVal [list]

    set remark1 [histogram -data $violList -title $histoTitle]
    lappend retVal $remark1

    if {$detailed} {
	set remark2 "Detailed violation report:"

	foreach curNViols $violList file $files {
	    appendLineToString remark2 [format "%s has %d violations" [lindex $file 0] $curNViols]
	}
	lappend retVal $remark2
    }

    return $retVal
}


#
# count the number of violated peaks in the given list of peaks.
# Uses the lowestViolation method, which ignores issues of 
# activation or peakAssignment choice.  It tries all peakAssignments
# and reports the minimum.  Unassigned peaks return 
# a bestViolation value of -1, which doesn't hurt this 
# calculation
#

proc numViolatedPeaks args {

    set peakList   [requiredFlagVal $args -peakList]
    set violCutoff [flagVal $args -violCutoff 0.5]

    set curNViols 0
    foreach curPeak $peakList {

	Peak -this $curPeak

	set curViol [$curPeak lowestViolation]
	if {$curViol > $violCutoff} {
	    incr curNViols
	}
	rename $curPeak ""
    }
    
    return $curNViols
}


proc numNonviolatedPeaks args {

    set peakList   [requiredFlagVal $args -peakList]
    set violCutoff [flagVal $args -violCutoff 0.5]

    set curNonViols 0
    foreach curPeak $peakList {

	Peak -this $curPeak

	if {[$curPeak isAssigned]} {
	    set curViol [$curPeak lowestViolation]
	    if {$curViol <= $violCutoff} {
		incr curNonViols
	    }
	}
	rename $curPeak ""
    }
    
    return $curNonViols
}


proc numActivePeaks args {
    
    set peakList [requiredFlagVal $args -peakList] 

    set curNActive 0
    foreach curPeak $peakList {
	Peak -this $curPeak
	if {[$curPeak isActive]} {
	    incr curNActive
	}
	rename $curPeak ""
    }

    return $curNActive
}


#
# Given a list of <filename, atomPosArray> pairs of converged structures 
# and a list of peaks, 
# use the coords to determine the likelihood of each peakAssignment
#
# Could be replaced with something that would calculate previous likelihood
# of any pair of ShiftAssignments, whether they're used in a PeakAssignment
# or not.
#

proc determinePeakAssignmentLikelihoodsFromConvergedStructs args {
    
    set files      [requiredFlagVal $args -files]
    set peakList   [requiredFlagVal $args -peakList]
    set violCutoff [flagVal $args -violCutoff 0.5]
    
    #
    # init array of times each PeakAssignment is violated
    #
    
    foreach curPeak $peakList {

	Peak -this $curPeak

	foreach curPA [$curPeak peakAssignments] {

	    PeakAssignment -this $curPA
	    set name [format "%s--%s" [$curPeak name] [$curPA name]]
	    set timesViolated($name) 0
	    rename $curPA ""
	}
	rename $curPeak ""
    }
    
    #
    # record how often each PeakAssignment is violated 
    # in the converged structures
    #

    set fileCount 0
    
    foreach file $files {
		
	xplorSim setAtomPosArr [lindex $file 1]

	foreach curPeak $peakList {

	    Peak -this $curPeak
	    set curUpBound  [$curPeak upBound]
	    set curLowBound [$curPeak lowBound]

	    foreach curPA [$curPeak peakAssignments] {

		PeakAssignment -this $curPA
		set name [format "%s--%s" [$curPeak name] [$curPA name]]
		set curViol [$curPA violation $curUpBound $curLowBound]
		if {$curViol > $violCutoff} {   
		    incr timesViolated($name)
		}
		rename $curPA ""
	    }
	    rename $curPeak ""
	}
    }
    
    #
    # convert timesViolated to likelihood for each peakAssignment
    #
    
    set numStructs [expr double([llength $files])]
    
    foreach curPeak $peakList {

	Peak -this $curPeak
	
	foreach curPA [$curPeak peakAssignments] {

	    PeakAssignment -this $curPA
	    set name [format "%s--%s" [$curPeak name] [$curPA name]]
	    set newL [expr ($numStructs - $timesViolated($name)) / $numStructs]
	    $curPA setPreviousLikelihood $newL
	    rename $curPA ""
	}
	rename $curPeak ""
    }

    return ""
}


proc reportNOEaccuracy args {

    set pot                  [requiredFlagVal $args -pot]
    set referenceStructFile  [requiredFlagVal $args -referenceStructureFilename]
    set violCutoff           [requiredFlagVal $args -violCutoff]
    set highLikelihoodCutoff [requiredFlagVal $args -highLikelihoodCutoff]
    set lowLikelihoodCutoff  [flagVal $args -lowLikelihoodCutoff [expr 1.0 - $highLikelihoodCutoff]]

    #
    # read in the reference coords
    #
	
    readPDB -fileName $referenceStructFile
	
    set retVal [list]

    #
    # select longrange peaks
    #

    set longRange   [selectNOEpeaks -from [$pot peaks] -proc {isLongRangePeak}]
    set nLR         [llength $longRange]

    #
    # Report fraction of high-likelihood longrange peaks that are good.
    #

    set  highLR     [selectNOEpeaks -from $longRange -proc [list likelihoodGreaterThanOrEqualTo $highLikelihoodCutoff]]
    set nHighLR     [llength $highLR]

    set nHighLRGood [llength [selectNOEpeaks -from $highLR -proc [list violLessThan $violCutoff]]]
    set nHighLRBad  [llength [selectNOEpeaks -from $highLR -proc [list violGreaterThan $violCutoff]]]

    if {$nHighLR > 0} {
	set fracHighLRGood [expr 100 * $nHighLRGood / double($nHighLR)]
	set fracHighLRBad  [expr 100 * $nHighLRBad  / double($nHighLR)]

	set line1 [format "%d of %d long range peaks have likelihood >= %f" \
		       $nHighLR $nLR $highLikelihoodCutoff]
	set line2 [format "%f %% of them agree with reference structure %s to within %f A" $fracHighLRGood $referenceStructFile $violCutoff]
	set line3 [format "%f %% of them do not agree with reference structure %s to within %f A" $fracHighLRBad $referenceStructFile $violCutoff]

	lappend retVal [format "%s\n%s\n%s" $line1 $line2 $line3]
    } else {
	lappend retVal "No long range peaks have likelihood > $highLikelihoodCutoff"
    }

    #
    # Report fraction of low-likelihood longrange peaks that are bad.
    #

    set  lowLR      [selectNOEpeaks -from $longRange -proc [list likelihoodLessThanOrEqualTo $lowLikelihoodCutoff]]
    set nLowLR      [llength $lowLR]

    set nLowLRGood  [llength [selectNOEpeaks -from $lowLR -proc [list violLessThan $violCutoff]]]
    set nLowLRBad   [llength [selectNOEpeaks -from $lowLR -proc [list violGreaterThan $violCutoff]]]

    if {$nLowLR > 0} {

	set fracLowLRGood [expr 100 * $nLowLRGood / double($nLowLR)]
	set fracLowLRBad  [expr 100 * $nLowLRBad  / double($nLowLR)]

	set line1 [format "%d of %d long range peaks have likelihood <= %f" \
		       $nLowLR $nLR $lowLikelihoodCutoff]
	set line2 [format "%f %% of them agree with reference structure %s to within %f A" $fracLowLRGood $referenceStructFile $violCutoff]
	set line3 [format "%f %% of them do not agree with reference structure %s to within %f A" $fracLowLRBad $referenceStructFile $violCutoff]
	
	lappend retVal [format "%s\n%s\n%s" $line1 $line2 $line3]
    } else {
	lappend retVal "No long range peaks have likelihood < $lowLikelihoodCutoff"
    }

    #
    # Report fraction of medium-likelihood longrange peaks that are good.
    #

    set  medLR      [selectNOEpeaks -from $longRange -proc [list likelihoodGreaterThan $lowLikelihoodCutoff]]
    set  medLR      [selectNOEpeaks -from $medLR     -proc [list likelihoodLessThan    $highLikelihoodCutoff]]
    set nMedLR      [llength $medLR]

    set nMedLRGood  [llength [selectNOEpeaks -from $medLR -proc [list violLessThan $violCutoff]]]
    set nMedLRBad   [llength [selectNOEpeaks -from $medLR -proc [list violGreaterThan $violCutoff]]]

    if {$nMedLR > 0} {
	set fracMedLRGood [expr 100 * $nMedLRGood / double($nMedLR)]
	set fracMedLRBad  [expr 100 * $nMedLRBad  / double($nMedLR)]

	set line1 [format "%d of %d long range peaks have likelihood between %f and  %f" \
		       $nMedLR $nLR $lowLikelihoodCutoff $highLikelihoodCutoff]
	set line2 [format "%f %% of them agree with reference structure %s to within %f A" $fracMedLRGood $referenceStructFile $violCutoff]
	set line3 [format "%f %% of them do not agree with reference structure %s to within %f A" $fracMedLRBad $referenceStructFile $violCutoff]

	lappend retVal [format "%s\n%s\n%s" $line1 $line2 $line3]
    } else {
	lappend retVal "No long range peaks have likelihood between $lowLikelihoodCutoff and $highLikelihoodCutoff"
    }

    #
    # Report fraction of good longrange peaks that are high-likelihood
    #

    set  goodLR     [selectNOEpeaks -from $longRange -proc [list violLessThan $violCutoff]]
    set nGoodLR     [llength $goodLR]

    set nGoodLRHigh [llength [selectNOEpeaks -from $goodLR -proc [list likelihoodGreaterThanOrEqualTo $highLikelihoodCutoff]]]
    set nGoodLRLow  [llength [selectNOEpeaks -from $goodLR -proc [list likelihoodLessThanOrEqualTo $lowLikelihoodCutoff]]]
    set nGoodLRMed  [expr $nGoodLR - ($nGoodLRHigh + $nGoodLRLow)]

    if {$nGoodLR > 0} {

	set fracGoodLRHigh [expr 100 * $nGoodLRHigh / double($nGoodLR)]
	set fracGoodLRLow  [expr 100 * $nGoodLRLow  / double($nGoodLR)]
	set fracGoodLRMed  [expr 100 * $nGoodLRMed  / double($nGoodLR)]

	set line1 [format "%d of %d long range peaks agree with reference structure %s to within %f A" \
		       $nGoodLR $nLR $referenceStructFile $violCutoff]
	set line2 [format "%f %% of them have likelihood >= %f" $fracGoodLRHigh $highLikelihoodCutoff]
	set line3 [format "%f %% of them have likelihood <= %f" $fracGoodLRLow  $lowLikelihoodCutoff]
	set line4 [format "%f %% of them have likelihood between %f and %f" $fracGoodLRMed  $lowLikelihoodCutoff $highLikelihoodCutoff]
	
	lappend retVal [format "%s\n%s\n%s\n%s" $line1 $line2 $line3 $line4]
    } else {
	lappend retVal [format "No long range peaks agree with reference structure %s to within %f A" \
			    $referenceStructFile $violCutoff]
    }

    #
    # Report fraction of bad longrange peaks that are low-likelihood
    #

    set badLR       [selectNOEpeaks -from $longRange -proc [list violGreaterThan $violCutoff]]
    set nBadLR      [llength $badLR]

    set nBadLRHigh  [llength [selectNOEpeaks -from $badLR -proc [list likelihoodGreaterThanOrEqualTo $highLikelihoodCutoff]]]
    set nBadLRLow   [llength [selectNOEpeaks -from $badLR -proc [list likelihoodLessThanOrEqualTo $lowLikelihoodCutoff]]]  
    set nBadLRMed   [expr $nBadLR - ($nBadLRHigh + $nBadLRLow)]

    if {$nBadLR > 0} {

	set fracBadLRHigh [expr 100 * $nBadLRHigh / double($nBadLR)]
	set fracBadLRLow  [expr 100 * $nBadLRLow  / double($nBadLR)]
	set fracBadLRMed  [expr 100 * $nBadLRMed  / double($nBadLR)]

	set line1 [format "%d of %d long range peaks do not agree with reference structure %s to within %f A" \
		       $nBadLR $nLR $referenceStructFile $violCutoff]
	set line2 [format "%f %% of them have likelihood <= %f" $fracBadLRLow $lowLikelihoodCutoff]
	set line3 [format "%f %% of them have likelihood >= %f" $fracBadLRHigh $highLikelihoodCutoff]
	set line4 [format "%f %% of them have likelihood between %f and %f" $fracBadLRMed  $lowLikelihoodCutoff $highLikelihoodCutoff]

	lappend retVal [format "%s\n%s\n%s\n%s" $line1 $line2 $line3 $line4]
    } else {
	lappend retVal [format "No long range peaks do not agree with reference structure %s to within %f A" \
			    $referenceStructFile $violCutoff]
    }


    #
    # simulate fbad-long w/ current likelihoods
    #

    lappend retVal [format "Frac bad long range forces = %f" [fracBadLRInfo -peakList [$pot peaks] -violCutoff $violCutoff -numIterations 10 -useLikelihoods]]
    
    #
    # select shortrange peaks
    #

    set shortRange   [selectNOEpeaks -from [$pot peaks] -proc {isShortRangePeak}]
    set nLR         [llength $shortRange]

    #
    # Report fraction of high-likelihood shortrange peaks that are good.
    #

    set  highLR     [selectNOEpeaks -from $shortRange -proc [list likelihoodGreaterThanOrEqualTo $highLikelihoodCutoff]]
    set nHighLR     [llength $highLR]

    set nHighLRGood [llength [selectNOEpeaks -from $highLR -proc [list violLessThan $violCutoff]]]
    set nHighLRBad  [llength [selectNOEpeaks -from $highLR -proc [list violGreaterThan $violCutoff]]]

    if {$nHighLR > 0} {
	set fracHighLRGood [expr 100 * $nHighLRGood / double($nHighLR)]
	set fracHighLRBad  [expr 100 * $nHighLRBad  / double($nHighLR)]

	set line1 [format "%d of %d short range peaks have likelihood >= %f" \
		       $nHighLR $nLR $highLikelihoodCutoff]
	set line2 [format "%f %% of them agree with reference structure %s to within %f A" $fracHighLRGood $referenceStructFile $violCutoff]
	set line3 [format "%f %% of them do not agree with reference structure %s to within %f A" $fracHighLRBad $referenceStructFile $violCutoff]

	lappend retVal [format "%s\n%s\n%s" $line1 $line2 $line3]
    } else {
	lappend retVal "No short range peaks have likelihood > $highLikelihoodCutoff"
    }

    #
    # Report fraction of low-likelihood shortrange peaks that are bad.
    #

    set  lowLR      [selectNOEpeaks -from $shortRange -proc [list likelihoodLessThanOrEqualTo $lowLikelihoodCutoff]]
    set nLowLR      [llength $lowLR]

    set nLowLRGood  [llength [selectNOEpeaks -from $lowLR -proc [list violLessThan $violCutoff]]]
    set nLowLRBad   [llength [selectNOEpeaks -from $lowLR -proc [list violGreaterThan $violCutoff]]]

    if {$nLowLR > 0} {

	set fracLowLRGood [expr 100 * $nLowLRGood / double($nLowLR)]
	set fracLowLRBad  [expr 100 * $nLowLRBad  / double($nLowLR)]

	set line1 [format "%d of %d short range peaks have likelihood <= %f" \
		       $nLowLR $nLR $lowLikelihoodCutoff]
	set line2 [format "%f %% of them agree with reference structure %s to within %f A" $fracLowLRGood $referenceStructFile $violCutoff]
	set line3 [format "%f %% of them do not agree with reference structure %s to within %f A" $fracLowLRBad $referenceStructFile $violCutoff]
	
	lappend retVal [format "%s\n%s\n%s" $line1 $line2 $line3]
    } else {
	lappend retVal "No short range peaks have likelihood < $lowLikelihoodCutoff"
    }

    #
    # Report fraction of medium-likelihood shortrange peaks that are good.
    #

    set  medLR      [selectNOEpeaks -from $shortRange -proc [list likelihoodGreaterThan $lowLikelihoodCutoff]]
    set  medLR      [selectNOEpeaks -from $medLR     -proc [list likelihoodLessThan    $highLikelihoodCutoff]]
    set nMedLR      [llength $medLR]

    set nMedLRGood  [llength [selectNOEpeaks -from $medLR -proc [list violLessThan $violCutoff]]]
    set nMedLRBad   [llength [selectNOEpeaks -from $medLR -proc [list violGreaterThan $violCutoff]]]

    if {$nMedLR > 0} {
	set fracMedLRGood [expr 100 * $nMedLRGood / double($nMedLR)]
	set fracMedLRBad  [expr 100 * $nMedLRBad  / double($nMedLR)]

	set line1 [format "%d of %d short range peaks have likelihood between %f and  %f" \
		       $nMedLR $nLR $lowLikelihoodCutoff $highLikelihoodCutoff]
	set line2 [format "%f %% of them agree with reference structure %s to within %f A" $fracMedLRGood $referenceStructFile $violCutoff]
	set line3 [format "%f %% of them do not agree with reference structure %s to within %f A" $fracMedLRBad $referenceStructFile $violCutoff]

	lappend retVal [format "%s\n%s\n%s" $line1 $line2 $line3]
    } else {
	lappend retVal "No short range peaks have likelihood between $lowLikelihoodCutoff and $highLikelihoodCutoff"
    }

    #
    # Report fraction of good shortrange peaks that are high-likelihood
    #

    set  goodLR     [selectNOEpeaks -from $shortRange -proc [list violLessThan $violCutoff]]
    set nGoodLR     [llength $goodLR]

    set nGoodLRHigh [llength [selectNOEpeaks -from $goodLR -proc [list likelihoodGreaterThanOrEqualTo $highLikelihoodCutoff]]]
    set nGoodLRLow  [llength [selectNOEpeaks -from $goodLR -proc [list likelihoodLessThanOrEqualTo $lowLikelihoodCutoff]]]
    set nGoodLRMed  [expr $nGoodLR - ($nGoodLRHigh + $nGoodLRLow)]

    if {$nGoodLR > 0} {

	set fracGoodLRHigh [expr 100 * $nGoodLRHigh / double($nGoodLR)]
	set fracGoodLRLow  [expr 100 * $nGoodLRLow  / double($nGoodLR)]
	set fracGoodLRMed  [expr 100 * $nGoodLRMed  / double($nGoodLR)]

	set line1 [format "%d of %d short range peaks agree with reference structure %s to within %f A" \
		       $nGoodLR $nLR $referenceStructFile $violCutoff]
	set line2 [format "%f %% of them have likelihood >= %f" $fracGoodLRHigh $highLikelihoodCutoff]
	set line3 [format "%f %% of them have likelihood <= %f" $fracGoodLRLow  $lowLikelihoodCutoff]
	set line4 [format "%f %% of them have likelihood between %f and %f" $fracGoodLRMed  $lowLikelihoodCutoff $highLikelihoodCutoff]
	
	lappend retVal [format "%s\n%s\n%s\n%s" $line1 $line2 $line3 $line4]
    } else {
	lappend retVal [format "No short range peaks agree with reference structure %s to within %f A" \
			    $referenceStructFile $violCutoff]
    }

    #
    # Report fraction of bad shortrange peaks that are low-likelihood
    #

    set badLR       [selectNOEpeaks -from $shortRange -proc [list violGreaterThan $violCutoff]]
    set nBadLR      [llength $badLR]

    set nBadLRHigh  [llength [selectNOEpeaks -from $badLR -proc [list likelihoodGreaterThanOrEqualTo $highLikelihoodCutoff]]]
    set nBadLRLow   [llength [selectNOEpeaks -from $badLR -proc [list likelihoodLessThanOrEqualTo $lowLikelihoodCutoff]]]  
    set nBadLRMed   [expr $nBadLR - ($nBadLRHigh + $nBadLRLow)]

    if {$nBadLR > 0} {

	set fracBadLRHigh [expr 100 * $nBadLRHigh / double($nBadLR)]
	set fracBadLRLow  [expr 100 * $nBadLRLow  / double($nBadLR)]
	set fracBadLRMed  [expr 100 * $nBadLRMed  / double($nBadLR)]

	set line1 [format "%d of %d short range peaks do not agree with reference structure %s to within %f A" \
		       $nBadLR $nLR $referenceStructFile $violCutoff]
	set line2 [format "%f %% of them have likelihood <= %f" $fracBadLRLow $lowLikelihoodCutoff]
	set line3 [format "%f %% of them have likelihood >= %f" $fracBadLRHigh $highLikelihoodCutoff]
	set line4 [format "%f %% of them have likelihood between %f and %f" $fracBadLRMed  $lowLikelihoodCutoff $highLikelihoodCutoff]

	lappend retVal [format "%s\n%s\n%s\n%s" $line1 $line2 $line3 $line4]
    } else {
	lappend retVal [format "No short range peaks do not agree with reference structure %s to within %f A" \
			    $referenceStructFile $violCutoff]
    }


    #
    # simulate fbad-short w/ current likelihoods
    #

    lappend retVal [format "Frac bad short range forces = %f" [fracBadLRInfo -peakList [$pot peaks] -violCutoff $violCutoff -numIterations 10 -useLikelihoods]]
    
    return $retVal
}


proc reportNOEprecision args {

    set pot                   [requiredFlagVal $args -pot]
    set heavyatomSelection    [requiredFlagVal $args -definedHeavyatoms]
    set highLikelihoodCutoff  [requiredFlagVal $args -highLikelihoodCutoff]
    set lowLikelihoodCutoff   [flagVal $args -lowLikelihoodCutoff \
				   [expr 1.0 - $highLikelihoodCutoff]]

    set retVal [list]

    #
    # extract the longrange likelihoods and plot their histogram
    #
    
    set longRange [selectNOEpeaks -from [$pot peaks] -proc {isLongRangePeak}]

    set temp [list]
    foreach curPeak $longRange {
	Peak -this $curPeak
	lappend temp [$curPeak previousLikelihood]
	rename $curPeak ""
    }

    lappend retVal [histogram -data $temp -min 0 -max 1 -title "Longrange peak likelihoods"]

    #
    # compute NOE precision (# high-likelihood peaks / residue)
    #

    set nHighLlongRange [llength [selectNOEpeaks -from $longRange -proc [list likelihoodGreaterThanOrEqualTo $highLikelihoodCutoff]]]
    set noePrec [expr $nHighLlongRange / double([numResiduesInSelection $heavyatomSelection])]

    #
    # compute NOE discrimination (fraction of LR peaks w/ high or low likelihood)
    #

    set nLowLlongRange [llength [selectNOEpeaks -from $longRange -proc [list likelihoodLessThanOrEqualTo $lowLikelihoodCutoff]]]
    if {[llength $longRange] > 0} {
	set noeDisc [expr ($nHighLlongRange + $nLowLlongRange) / double([llength $longRange])]
    } else {
	set noeDisc 0
    }
   
    set line1 [format "Defined region of structure is %s" [$heavyatomSelection string]]
    set line2 [format "Number of high-likelihood long range peaks/residue in defined region: %f" $noePrec]
    lappend retVal [format "%s\n%s" $line1 $line2]

    lappend retVal [format "Long-range NOE discrimination: %f %%" [expr 100.0 * $noeDisc]]

    return $retVal
}


# 
# Standard summary of the results of a Marvin pass.
# Given PDB filenames & original shiftAssign/peak filenames, 
# find converged structures, 
# generate new likelihoods,
# write new shiftAssign/peak files,
# and generate average structure
#

proc summarizeMarvinResults args {

    set fnames                    [requiredFlagVal $args -fileNames]
    set origPeakFilenames         [requiredFlagVal $args -origPeakFile]
    set completePeakFilenames     [flagVal $args -completePeakFile]
    set exceptionFilenames        [flagVal $args -explicitExceptionFile]
    set origSAfilenames           [requiredFlagVal $args -origShiftAssignmentFile]
    set newPeakFilenames          [requiredFlagVal $args -newPeakFile]
    set newSAfilenames            [requiredFlagVal $args -newShiftAssignmentFile]
    set avgPDBfile                [requiredFlagVal $args -avgPDBfile]
    set referenceStructFilename   [flagVal $args -referenceStructureFile ""]
    set violCutoff                [flagVal $args -violCutoff 0.5]
    set fracConverged             [flagVal $args -fracConverged 0.1]
    set highLikelihoodCutoff      [flagVal $args -highLikelihoodCutoff 0.9]
    set noeCompletenessCutoff     [flagVal $args -noeCompletenessCutoff 0.05]
    set scatterCutoff             [flagVal $args -scatterCutoff 0.05]
    set invBound                  [flagVal $args -inverseBound 4.0]
    set invMeth                   [flagVal $args -inverseMethylCorrection 0.0]
    set SAsForCompleteness        [flagVal $args -shiftAssignmentFilesToUseForNOECompleteness]
    set completenessWeight        [flagVal $args -completenessWeight 0.0]
    set completePeakHighLikelihoodCutoff [flagVal $args -completePeakHighLikelihoodCutoff 0.8]
    set completePeakViolCutoff           [flagVal $args -completePeakViolCutoff 5.0]

    # change default definedBackbone to include DNA/RNA atoms?  Which atoms?

    set backboneSelection  [flagVal $args -definedBackbone   [AtomSel -args "(name ca or name c or name n)"]]
    set heavyatomSelection [flagVal $args -definedHeavyatoms [AtomSel -args "(not name h*)"]]


    #
    # read the PDB files to be analyzed
    #

    updateUser [format "Reading all PDB files\n"]

    set allFiles [eval grabPDBFiles $args]

    #
    # read original shiftAssignment and peak files
    #

    set potCount 0
    set marvPots [list]
    set completenessPots [list]

    updateUser [format "Reading all shiftAssignment and peak files\n"]

    foreach origSAfilename $origSAfilenames origPeakFilename $origPeakFilenames exceptionFilename $exceptionFilenames {

	set curPot [create_MarvinNOEPotential [format "marv%d" [incr potCount]]]

	readShiftAssignments \
	    -fileName $origSAfilename \
	    -pot $curPot

	readMarvinPeaks \
	    -fileName $origPeakFilename \
	    -pot $curPot

	if {$exceptionFilename != ""} {
	    puts "reading exceptions"
	    readExplicitInverseExceptions \
		-fileName $exceptionFilename \
		-pot $curPot
	}


	lappend marvPots $curPot

	if {[lsearch -exact $SAsForCompleteness $origSAfilename] != -1} {
	    lappend completenessPots $curPot
	}
    }

    updateUser [format "\nSelecting converged structures on the basis of long range peaks from all spectra and completeness from defined spectra\n"]

    marvinPyth command {import pasd}
    set nMono [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.nMono} \
			   {} {out}] \
		      0] 1]
    set aveExp [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.aveExp} \
			   {} {out}] \
		      0] 1]


    # 
    # select long range peaks from the original peak file(s) for convergence evaluation.
    # Uses the default minimum primary sequence distance to define
    # "long range"
    #

    set allLongRange [list]

    foreach curPot $marvPots {
	$curPot updatePrimarySeqDists
	set curLongRange [selectNOEpeaks -from [$curPot peaks] -proc {isLongRangePeak}]
	set allLongRange [concat $allLongRange $curLongRange]
    }

    #
    # select the converged structures.
    #

    set goodFiles [chooseBestFraction -files $allFiles \
		       -peakList $allLongRange \
		       -violCutoff $violCutoff \
		       -completenessPots $completenessPots \
		       -inverseBound $invBound \
		       -inverseMethylCorrection $invMeth \
		       -completenessWeight $completenessWeight \
		       -frac $fracConverged]

    #
    # report statistics of peak violations used to select converged structs
    #

    set ensemblePeakRemarks [list]

    lappend ensemblePeakRemarks [format "Choosing best %d of %d total structs, based on num longrange peaks with violations < %f A" \
				     [llength $goodFiles] [llength $allFiles] $violCutoff]
    
    set ensemblePeakRemarks [concat $ensemblePeakRemarks [reportEnsemblePeakViolations -files $allFiles  -peakList $allLongRange -violCutoff $violCutoff -histogramTitle "Num violated longrange peaks in all structures"]]
    set ensemblePeakRemarks [concat $ensemblePeakRemarks [reportEnsemblePeakViolations -files $goodFiles -peakList $allLongRange -violCutoff $violCutoff -histogramTitle "Num violated longrange peaks in converged structures" -detailed ]]



    #
    # read in complete peak files if they're defined
    #

    updateUser "Adding high-likelihood peakAssignments from complete-peak-assignment files that were previously filtered out\n"

    set pasAddedFromCompletePeakTables [list]
	
    foreach curPot $marvPots origSAfilename $origSAfilenames completePeakFilename $completePeakFilenames {

	if {$completePeakFilename == ""} {
	    continue
	}

	#
	# record the PAs in the filtered pot to 
	# ensure that we don't add PAs that are equivalent to
	# existing ones
	#

	array unset equivPAs 

	foreach curPeak [$curPot peaks] {
	    Peak -this $curPeak
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA 
		set curIndex [join [list [$curPA fromAssignmentName] [$curPA toAssignmentName]]]
		set equivPAs($curIndex) 1
		rename $curPA ""
	    }
	    rename $curPeak ""
	}
	    
	#
	# create a new pot to hold the complete-PA data and read them in
	#
	
	set completePot [create_MarvinNOEPotential [format "marv%d" [incr potCount]]]
	
	readShiftAssignments \
	    -fileName $origSAfilename \
	    -pot $completePot
	
	readMarvinPeaks \
	    -fileName $completePeakFilename \
	    -pot $completePot

	# 
	# calculate likelihoods for all peakAssigns in the filtered-PA pot (in order to identify peaks that already have a high-l PA)
	#
	
	determinePeakAssignmentLikelihoodsFromConvergedStructs \
	    -files $goodFiles \
	    -peakList [$curPot peaks] \
	    -violCutoff $violCutoff 

	
	# 
	# calculate likelihoods for all peakAssigns in the complete-PA pot (in order to identify high-likelihood ones)
	#
	
	determinePeakAssignmentLikelihoodsFromConvergedStructs \
	    -files $goodFiles \
	    -peakList [$completePot peaks] \
	    -violCutoff $completePeakViolCutoff 
		
	foreach curCompletePeak [$completePot peaks] {
	    Peak -this $curCompletePeak

	    foreach curCompletePA [$curCompletePeak peakAssignments] {
		PeakAssignment -this $curCompletePA 

		if {[$curCompletePA previousLikelihood] >= $completePeakHighLikelihoodCutoff} {
		    
		    #
		    # look up the corresponding peak/peakAssignment in the filtered-PA pot
		    #
		    
		    if {! [$curPot hasPeakNamed [$curCompletePeak name]]} {
			error [format "Peak %s exists in complete-peak-assignment file, but not filtered-peak-assignment file" [$curCompletePeak name]]
		    }
		    
		    set curPeak [$curPot peakNamed [$curCompletePeak name]]
		    
		    #
		    # if the filtered-PA peak doesn't already have a high-likelihood PA,
		    #

		    if {[$curPeak previousLikelihood] < $highLikelihoodCutoff} {

			#
			# make sure the new PA from the complete table is not equivalent to an existing PA anywhere in the filtered table
			#

			set curIndex [join [list [$curCompletePA fromAssignmentName] [$curCompletePA toAssignmentName]]]
			if {! [info exists equivPAs($curIndex)]} {
			    

			    set newPA [PeakAssignment -args [$curCompletePA name]]
			    $newPA setNMono $nMono
			    $newPA setAveExp $aveExp

			    if {! [$curPot hasShiftAssignmentNamed [$curCompletePA fromAssignmentName]]} {
				error [format "Missing shiftAssignment %s" $fromAssignmentName]
			    }
			    if {! [$curPot hasShiftAssignmentNamed [$curCompletePA toAssignmentName]]} {
				error [format "Missing shiftAssignment %s" $toAssignmentName]
			    }
			    
			    set fromSA [$curPot shiftAssignmentNamed [$curCompletePA fromAssignmentName]]
			    set toSA   [$curPot shiftAssignmentNamed [$curCompletePA toAssignmentName]]
			    ShiftAssignment -this $fromSA
			    ShiftAssignment -this $toSA
			    $newPA setFromAssignment [$fromSA cget -this] 
			    $newPA setToAssignment   [$toSA cget -this]
			    $newPA setUpBoundCorrection  [$curCompletePA upBoundCorrection]
			    $newPA setLowBoundCorrection [$curCompletePA lowBoundCorrection]
			    $newPA setNote [$curCompletePA note]
			    $newPA appendToNote [format "Added from %s during summarization" $completePeakFilename]
			    
			    if {[$curCompletePA isGood]} {
				$newPA setGood
			    } else {
				$newPA setBad
			    }
			    
			    lappend pasAddedFromCompletePeakTables [list [$curPot instanceName] [$curPeak name] [$newPA name] [$newPA isGood] [$newPA isLongRange 6] [$curPeak isGood]]
			    
			    $curPeak addPeakAssignment [$newPA cget -this]
			    $newPA -disown
			    rename $newPA ""
			    rename $fromSA ""
			    rename $toSA ""
			}
		    }
		    rename $curPeak ""
		}
		rename $curCompletePA ""
	    }
	    rename $curCompletePeak ""
	}
	
	rename $completePot ""
    }

	    
    updateUser [format "Calculating new likelihoods for all peakAssignments and shiftAssignments for each spectrum\n"]

    foreach curPot $marvPots newPeakFilename $newPeakFilenames newSAfilename $newSAfilenames {

	set peakRemarks [list]

	#
	# report any peakAssigns that were added from the raw peak files
	#

	set curAddedPAs [list]
	foreach elem $pasAddedFromCompletePeakTables {
	    if {[lindex $elem 0] == [$curPot instanceName]} {
		lappend curAddedPAs $elem
	    }
	}

	if {[llength $curAddedPAs] > 0} {
	    set tempRem [format "Added %d high-likelihood peakAssignments from the unfiltered peak table:" [llength $curAddedPAs]]
	    set nGoodAddedPAs 0
	    set nPAsAddedToBadPeaks 0
	    set nLRAddedPAs 0
	    set nGoodLRAddedPAs 0
	    set nGoodPAsAddedToBadPeaks 0
	    foreach elem $curAddedPAs {
		if {[lindex $elem 3]} {
		    incr nGoodAddedPAs
		}
		if {[lindex $elem 4]} {
		    incr nLRAddedPAs
		}
		if {[lindex $elem 3] && [lindex $elem 4]} {
		    incr nGoodLRAddedPAs
		}
		if {! [lindex $elem 5]} {
		    incr nPAsAddedToBadPeaks
		    if {[lindex $elem 3]} {
			incr nGoodPAsAddedToBadPeaks
		    }
		}
	    }
	    appendLineToString tempRem [format "   %d of the added peak assignments are marked good" $nGoodAddedPAs]
	    appendLineToString tempRem [format "   %d of the added peak assignments are long range" $nLRAddedPAs]
	    appendLineToString tempRem [format "      Of the added long range peak assignments, %d are marked good" $nGoodLRAddedPAs]
	    appendLineToString tempRem [format "   %d peak assignments were added to peaks that had no marked-good peakAssignments available" $nPAsAddedToBadPeaks]
	    appendLineToString tempRem [format "      Of the peak assignments added to bad peaks, %d are marked good" $nGoodPAsAddedToBadPeaks]
	    lappend peakRemarks $tempRem
	}

	#
	# calculate new peakAssignment likelihoods
	#
	
	determinePeakAssignmentLikelihoodsFromConvergedStructs \
	    -files $goodFiles \
	    -peakList [$curPot peaks] \
	    -violCutoff $violCutoff 

	#
	# report statistics on new peakAssignment likelihoods
	#
	
	set peakRemarks [concat $ensemblePeakRemarks $peakRemarks]
	
	set peakRemarks [concat $peakRemarks [reportNOEprecision -pot $curPot -definedHeavyatoms $heavyatomSelection -highLikelihoodCutoff $highLikelihoodCutoff]]
	
	if {$referenceStructFilename != ""} {
	    set peakRemarks [concat $peakRemarks [reportNOEaccuracy -pot $curPot -referenceStructureFilename $referenceStructFilename -violCutoff $violCutoff -highLikelihoodCutoff $highLikelihoodCutoff]]
	}

	#
	# tack in more statistics
	#

	initialPeakAnalysis \
	    -pot $curPot \
	    -referenceStructureFile $referenceStructFilename \
	    -violCutoff $violCutoff \
	    -minLikelihood $highLikelihoodCutoff \
	    -description "pass summary" \
	    -remarksVariableName peakRemarks
	
	#
	# write out the new peak file
	#
	
	writeMarvinPeaks \
	    -fileName $newPeakFilename \
	    -peakList [$curPot peaks] \
	    -remarks $peakRemarks

	#
	# calc new ShiftAssignment likelihoods
	#
	
	set saRemarks   [list]
	lappend saRemarks [format "Calculating shiftAssignment likelihoods with noeCompleteness cutoff of %f" \
			       $noeCompletenessCutoff]
	lappend saRemarks [format "And proton peak position scatter cutoff of %f ppm" $scatterCutoff]
	lappend saRemarks [format "And NOE violation cutoff of %f A" $violCutoff]
	lappend saRemarks [format "And inverse bound of %f A and inverse methyl correction of %f A" \
			       $invBound $invMeth]
    

	#
	# calculate new shiftAssignment likelihoods
	#

	determineSALikelihoodsFromConvergedStructs \
	    -files $goodFiles \
	    -pot $curPot \
	    -violCutoff $violCutoff \
	    -noeCompletenessCutoff $noeCompletenessCutoff \
	    -scatterCutoff $scatterCutoff \
	    -inverseBound $invBound \
	    -inverseMethylCorrection $invMeth 

	lappend saRemarks [reportShiftAssignLikelihoods -pot $curPot]

	writeShiftAssignments \
	    -fileName $newSAfilename \
	    -shiftAssignments [$curPot shiftAssignments] \
	    -remarks $saRemarks

    }

    #
    # calc minimized mean struct
    #

    updateUser [format "Calculating mean structure\n"]

    XplorCommand "set message off echo off print off end" 

    set pdbRemarks [list]
    
    lappend pdbRemarks [format "Note:  Backbone atoms in defined structure region are\n (%s)" \
			    [$backboneSelection string]]
    lappend pdbRemarks [format "Note:  Heavy atoms in defined structure region are\n (%s)" \
			    [$heavyatomSelection string]]

    meanStruct -files $goodFiles -selection $backboneSelection
    
    set curRem "minimized average PDB file created from these files:"
    foreach curFile $goodFiles {
	set curRem [format "%s\n   %s" $curRem [lindex $curFile 0]]
    }
    
    lappend pdbRemarks $curRem
    
    #
    # check coordinate accuracy, if there's a reference structure
    #
    
    if {$referenceStructFilename != ""} {
	
	set heavyDefinedRMSD [compareToReference -referenceFileName $referenceStructFilename \
				  -selection $heavyatomSelection]
	
	set bbnDefinedRMSD [compareToReference -referenceFileName $referenceStructFilename \
				-selection $backboneSelection]
	
	set line1 [format "Accuracy of defined structure region of current coords vs %s :" \
		       $referenceStructFilename]
	
	set line2 [format "   %f (backbone)   %f (heavy atoms)"\
		       $bbnDefinedRMSD $heavyDefinedRMSD]
	
	lappend pdbRemarks [format "%s\n%s" $line1 $line2]
    }
    

    # check structure precision
    
    lappend pdbRemarks [structPrecision -files $goodFiles -verbose \
			    -backboneSelection $backboneSelection \
			    -heavyatomSelection $heavyatomSelection]
    
    # write mean coords
    
    writePDB -fileName $avgPDBfile -remarks $pdbRemarks
    
    # clean up

    foreach curPot $marvPots {
	rename $curPot ""
    }
    
    return ""
}



proc fracBadLRInfo args {

    set peakList         [requiredFlagVal $args -peakList]
    set violCutoff       [flagVal $args -violCutoff 0.5]
    set useGoodFlag      [flagExists $args -useGoodBadFlag]
    set longRangeCutoff  [flagVal $args -longRangeCutoff 6]
    set minLikelihood    [flagVal $args -minPreviousLikelihood -1]
    set maxLikelihood    [flagVal $args -maxPreviousLikelihood 2]
    set minFiltersFailed [flagVal $args -minFiltersFailed -1]
    set maxFiltersFailed [flagVal $args -maxFiltersFailed 9999]
    set nIters           [flagVal $args -numIterations 1]
    set useLikelihoods   [flagExists $args -useLikelihoods]
    
    #
    # determines the fraction of the longrange forces that are good, based either on current coords 
    # or on the good/bad flag of each PeakAssignment (if called with -useGoodBadFlag)
    #
    # by default, all peakAssignments are presumed to be active.  But if -useLikelihoods is set, 
    # then a random set of active peakAssignments is generated using the previousLikelihoods, and 
    # fbadlong is calculated using those active peakAssignments.  Setting -numIterations changes the
    # number of such random activation sets that are generated, and their mean fbadlong is returned.
    #
    # If min/maxPreviousLikelihood is set, then activation is limited to peakAssignments with 
    # previousLikelihoods < min and/or > max.  
    #
    # Note that peakAssignments with no value for the previousLikelihood are always
    # be presumed to be active, so that for scripts that don't calculate previousLikelihoods, all 
    # peakAssignments will be assumed to be active
    #

    if {$nIters < 1} {
	error "fracBadLRInfo called with numIterations < 1"
    }

    set totFbadLong 0
    for {set count 0} {$count < $nIters} {incr count} {

	set totGoodLRforce 0
	set totLRforce 0

	foreach curPeak $peakList {

	    Peak -this $curPeak

	    set upBound  [$curPeak upBound]
	    set lowBound [$curPeak lowBound]
	    
	    set nCurActivePAs 0
	    set nCurActiveLRPAs 0
	    set nCurActiveGoodLRPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		
		PeakAssignment -this $curPA 
		
		set curPAactive 0

		if {([$curPA numFiltersFailed] >= $minFiltersFailed) && 
		    ([$curPA numFiltersFailed] <= $maxFiltersFailed)} {

		    if { (! [$curPA hasPreviousLikelihood]) ||
			  (([$curPA previousLikelihood] >= $minLikelihood) &&
			   ([$curPA previousLikelihood] <= $maxLikelihood) &&
			   (! $useLikelihoods || ([$curPA previousLikelihood] > [uniformRandom]))) } {
			       
			set curPAactive 1
		    }
		}
		
		if {$curPAactive} {
		    incr nCurActivePAs
		    
		    if {[$curPA isLongRange $longRangeCutoff]} {
			incr nCurActiveLRPAs
			
			if {$useGoodFlag} {
			    if {[$curPA isGood]} {
				incr nCurActiveGoodLRPAs
			    }
			} else {
			    if {[$curPA violation $upBound $lowBound] <= $violCutoff} {
				incr nCurActiveGoodLRPAs
			    }
			}
		    }
		}
		
		rename $curPA ""
	    }
	    
	    if {$nCurActivePAs > 0} {
		set forcePerAssign [expr 1.0 / double($nCurActivePAs)]
		set totLRforce     [expr $totLRforce     + ($nCurActiveLRPAs     * $forcePerAssign)]
		set totGoodLRforce [expr $totGoodLRforce + ($nCurActiveGoodLRPAs * $forcePerAssign)]
	    }
	    
	    rename $curPeak ""
	}
	
	if {$totLRforce == 0} {
	    set curFbadLong 1
	} else {
	    set curFbadLong [expr 1 - ($totGoodLRforce / (double($totLRforce)))]
	}

	set totFbadLong [expr $totFbadLong + $curFbadLong]
    }

    return [expr $totFbadLong / double($nIters)]
}


proc markGoodPeakAssignments args {

    set violCutoff [requiredFlagVal $args -violCutoff 0.5]
    set peakList   [requiredFlagVal $args -peakList]
    set remVar     [flagVal $args -remarksVariableName ""]

    #
    # set the official "good/bad" flag for each peakAssignment
    #

    set nBadPAs 0
    set nGoodPAs 0
    set nBadPeaks 0
    set nGoodPeaks 0
    set nUnassignedPeaks 0

    foreach p $peakList {
	    
	Peak -this $p	    		
	    
	set curPeakGood 0
	
	foreach pa [$p peakAssignments] {		
		
	    PeakAssignment -this $pa
	    set v [$pa violation [$p upBound] [$p lowBound]]

	    if {($v != -1) && ($v < $violCutoff)} {
		$pa setGood
		set curPeakGood 1
		incr nGoodPAs
	    } else {
		$pa setBad
		incr nBadPAs
	    }
		
	    rename $pa ""
	}

	if {! [$p isAssigned]} {
	    incr nUnassignedPeaks
	} elseif {$curPeakGood} {
	    incr nGoodPeaks
	} else {
	    incr nBadPeaks
	}

	rename $p ""
    }
    
    if {$remVar != ""} {
	upvar $remVar r
	set curRem [format "MarkGoodPeakAssignments:  using distance violation cutoff of %f A" $violCutoff]
	appendLineToString curRem [format "  %d peakAssignments were marked good, %d were marked bad" $nGoodPAs $nBadPAs]
	appendLineToString curRem [format "  %d peaks had a good peakAssignment, %d peaks did not" $nGoodPeaks $nBadPeaks]
	appendLineToString curRem [format "  %d peaks had no peakAssignments at all" $nUnassignedPeaks]
	lappend r $curRem
    }
}
    





proc initialPeakAnalysis args {
    
    set pot                 [flagVal $args -pot ""]
    set allPeaks            [flagVal $args -peakList [list]]
    set referenceStructFile [flagVal $args -referenceStructureFile ""]
    set violCutoff          [flagVal $args -violCutoff 0.5]
    set longRangeCutoff     [flagVal $args -longRangeCutoff 5]
    set remVar              [flagVal $args -remarksVariableName ""]
    set minLikelihoodCutoff [flagVal $args -minLikelihood -1]
    set maxLikelihoodCutoff [flagVal $args -maxLikelihood 2]
    set description         [flagVal $args -description "initialMatch"]
    set isVerbose           [flagExists $args -verbose]

    set remList [list]

    #
    # Note the previous likelihood range we're using
    #

    lappend remList [format "InitialPeakAnalysis at step %s, only including peakAssignments with previousLikelihoods in range %f .. %f" $description $minLikelihoodCutoff $maxLikelihoodCutoff]

    if {$pot != ""} {
	set allPeaks [$pot peaks]
    } 

    #
    # Gather up the peakAssignments with previousLikelihoods in the prescribed
    # range
    #

    set peakData [list]

    foreach curPeak $allPeaks {
	Peak -this $curPeak 
	
	set curPAdata [list]
	set curPSDs   [list]

	foreach curPA [$curPeak peakAssignments] {
	    PeakAssignment -this $curPA

	    $curPA updatePrimarySeqDist 

	    if {([$curPA previousLikelihood] >= $minLikelihoodCutoff) &&
		([$curPA previousLikelihood] <= $maxLikelihoodCutoff)} {

		lappend curPAdata [list [$curPA name] [$curPA cget -this] [$curPA primarySeqDist] [$curPA previousLikelihood]]
		lappend curPSDs [$curPA primarySeqDist]
	    }    
	    rename $curPA ""
	}

	if {[llength $curPSDs] == 0} {
	    set curPSD -1
	} else {
	    set curPSD [min $curPSDs]
	}
	
	lappend peakData [list [$curPeak name] [$curPeak cget -this] $curPSD $curPAdata]
	rename $curPeak ""
    }

    #
    # mark the primary sequence distance class of each selected
    # peak and peakAssignment
    #

    if {$isVerbose} {

	foreach peakDatum $peakData {
	    
	    set curPeak    [lindex $peakDatum 1]
	    set curPeakPSD [lindex $peakDatum 2]
	    set curPAData  [lindex $peakDatum end]
	    
	    Peak -this $curPeak
	    
	    foreach paDatum $curPAData {
		
		set curPA [lindex $paDatum 1]
		PeakAssignment -this $curPA
		
		if {[$curPA isIntraresidue]} {
		    $curPA appendToNote "intraresidue"
		    
		} elseif {[$curPA isSequential]} {
		    $curPA appendToNote "sequential"
		    
		} elseif {[$curPA isShortRange $longRangeCutoff]} {
		    $curPA appendToNote "short range"
		    
		} else {
		    $curPA appendToNote "long range"
		}
		
		rename $curPA ""
	    }
	    
	    if {[llength $curPAData] == 0} {
		$curPeak appendToNote [format "No peakAssignments in previousLikelihood range %f  .. %f" $minLikelihoodCutoff $maxLikelihoodCutoff]
	    } else {
		
		if {$curPeakPSD == 0} {
		    $curPeak appendToNote [format "intraresidue, with previousLikelihood range %f  .. %f" $minLikelihoodCutoff $maxLikelihoodCutoff]
		} elseif {$curPeakPSD == 1} {
		    $curPeak appendToNote [format "sequential, with previousLikelihood range %f  .. %f" $minLikelihoodCutoff $maxLikelihoodCutoff]
		} elseif {$curPeakPSD < $longRangeCutoff} {
		    $curPeak appendToNote [format "short range, with previousLikelihood range %f  .. %f" $minLikelihoodCutoff $maxLikelihoodCutoff]
		} else {
		    $curPeak appendToNote [format "long range, with previousLikelihood range %f  .. %f" $minLikelihoodCutoff $maxLikelihoodCutoff]
		}
	    }
	    
	    rename $curPeak ""
	}
    }
	
    #
    # note how many peaks are unassigned
    #

    set unassigned [list]
    foreach peakDatum $peakData {
	if {[llength [lindex $peakDatum end]] == 0} {
	    lappend unassigned $peakDatum
	}
    }
    
    lappend remList [format "Number of unassigned peaks: %d (%f%% of total)" \
			 [llength $unassigned] [expr 100 * [llength $unassigned] / double([llength $peakData])]]

    if { [llength $unassigned] == [llength $peakData] } {
	puts "initialPeakAnalysis: Error. No peaks were assigned!"
    }
    
		     
    #
    # divide peaks into primary sequence distance ranges to get statistics
    #

    set longRange  [list]
    set shortRange [list]
    set sequential [list]
    set intrares   [list]

    foreach peakDatum $peakData {
	set curPSD [lindex $peakDatum 2]
	if {$curPSD == -1} {
	    continue
	}
	
	if {$curPSD == 0} {
	    lappend intrares $peakDatum
	} elseif {$curPSD == 1} {
	    lappend sequential $peakDatum
	} elseif {$curPSD < $longRangeCutoff} {
	    lappend shortRange $peakDatum
	} else {
	    lappend longRange $peakDatum
	}
    }

        
    lappend remList [format "Number of intraresidue peaks: %d" [llength $intrares]]
    lappend remList [format "Number of sequential peaks: %d"   [llength $sequential]]
    lappend remList [format "Number of short range peaks: %d"  [llength $shortRange]]
    lappend remList [format "Number of long range peaks: %d"   [llength $longRange]]
    
    #
    # mark each peak with its degeneracy
    # and gather statistics
    #
        
    set allDegens [list]
    set nonZeroDegens [list]
    set numSingleAssignPeaks 0

    foreach peakDatum $peakData {
	set curDegen [llength [lindex $peakDatum end]]
	set curPeak  [lindex $peakDatum 1]

	Peak -this $curPeak

	$curPeak appendToNote [format "degeneracy %d, with previousLikelihood range %f  .. %f" $curDegen $minLikelihoodCutoff $maxLikelihoodCutoff]

	lappend allDegens $curDegen

	if {$curDegen > 0} {
	    lappend nonZeroDegens $curDegen
	} 

	if {$curDegen == 1} {
	    incr numSingleAssignPeaks
	}

	rename $curPeak ""
    }
    
    set m   [mean $nonZeroDegens]
    set sd  [standardDeviation $nonZeroDegens]
    set tot [sum $allDegens]
    
    lappend remList [histogram -data $allDegens -min 0 -nBins [max $allDegens] -title "Distribution of degeneracy"]

    lappend remList [format "Total number of peak assignments: %d" $tot]
    lappend remList [format "Number of peakAssignments per assigned peak: %f +/- %f" $m $sd]
    
    lappend remList [format "Number of single-assignment peaks: %d" $numSingleAssignPeaks]
    
    #
    # if we have a reference structure, generate stats based on it
    #
    
    if {$referenceStructFile != ""} {
	
	lappend remList "reference structure is $referenceStructFile"
	
	readPDB -fileName $referenceStructFile	

	# set the official "good/bad" flag for each peakAssignment

	markGoodPeakAssignments -violCutoff $violCutoff -peakList $allPeaks \
	    -remarksVariableName remList

	#
	# mark each peak and peakAssignment as good or bad, with degree of violation
	#

	set newPeakData [list]
	    
	foreach peakDatum $peakData {

	    set curPeak    [lindex $peakDatum 1]
	    set curPeakPSD [lindex $peakDatum 2]
	    set curPAdata  [lindex $peakDatum end]

	    Peak -this $curPeak

	    set curPeakViols [list]
	    set newPAData    [list]

	    foreach paDatum $curPAdata {
		set curPA [lindex $paDatum 1]
		PeakAssignment -this $curPA 

		set v [$curPA violation [$curPeak upBound] [$curPeak lowBound]]
		
		lappend curPeakViols $v
		lappend newPAData \
		    [list [$curPA name] [$curPA cget -this] \
			 [$curPA primarySeqDist] \
			 [$curPA previousLikelihood] $v]
		
		if {$isVerbose} {
		    if {$v > $violCutoff} {
			$curPA appendToNote [format "Bad.  Violation in reference struct is %f" $v]
		    } elseif {$v != -1} {
			$curPA appendToNote "Good."
		    }
		}

		rename $curPA ""
	    }

	    if {[llength $curPeakViols] == 0} {
		set curPeakMinViol -1
	    } else {
		set curPeakMinViol [min $curPeakViols]
	    }

	    lappend newPeakData [list [$curPeak name] [$curPeak cget -this] $curPeakPSD $curPeakMinViol $newPAData]

	    if {$isVerbose} {
		if {$curPeakMinViol > $violCutoff} {
		    $curPeak appendToNote \
			[format "Bad.  Lowest violation among peakAssignments with previousLikelihoods in range %f .. %f is %f" \
			     $minLikelihoodCutoff $maxLikelihoodCutoff \
			     $curPeakMinViol]
		} else {
		    $curPeak appendToNote "Good."
		}
	    }

	    rename $curPeak ""
	}

	#
	# analyze good/bad peaks as a function of primary sequence distance range
	#

	set glr [list]
	set gsr [list]
	set gsq [list]
	set gir [list]

	foreach peakDatum $newPeakData {
	    set curPSD  [lindex $peakDatum 2] 
	    if {$curPSD == -1} {
		continue
	    }

	    set curViol [lindex $peakDatum 3]

	    if {$curViol < $violCutoff} {
		if {$curPSD == 0} {
		    lappend gir $peakDatum
		} elseif {$curPSD == 1} {
		    lappend gsq $peakDatum
		} elseif {$curPSD < $longRangeCutoff} {
		    lappend gsr $peakDatum
		} else {
		    lappend glr $peakDatum
		}
	    }
	}
	   	
	if {[llength $longRange] > 0} {
	    set goodlrpct [expr 100.0 * ([llength $glr] / double([llength $longRange]))]
	} else {
	    set goodlrpct 0.0
	}
	
	if {[llength $shortRange] > 0} {
	    set goodsrpct [expr 100.0 * ([llength $gsr] / double([llength $shortRange]))]
	} else {
	    set goodsrpct 0.0
	}
	
	if {[llength $sequential] > 0} {
	    set goodsqpct [expr 100.0 * ([llength $gsq] / double([llength $sequential]))]
	} else {
	    set goodsqpct 0.0
	}
	
	if {[llength $intrares] > 0} {
	    set goodirpct [expr 100.0 * ([llength $gir] / double([llength $intrares]))]
	} else {
	    set goodirpct 0.0
	}
	
	set overall [concat $intrares $sequential $shortRange $longRange]
	set goodOverall [concat $gir $gsq $gsr $glr]

	if {[llength $overall] > 0} {
	    set goodOverallPct [expr 100.0 * ([llength $goodOverall] / double([llength $overall]))]
	} else {
	    set goodOverallPct 0.0
	}
		
	set badInfo [fracBadLRInfo -peakList $allPeaks \
			 -minPreviousLikelihood $minLikelihoodCutoff \
			 -maxPreviousLikelihood $maxLikelihoodCutoff \
			 -violCutoff $violCutoff \
			 -longRangeCutoff $longRangeCutoff]
	
	lappend remList [format "Number of good intraresidue peaks: %d of %d (%f%% good)" \
			     [llength $gir] [llength $intrares] $goodirpct]
	lappend remList [format "Number of good sequential peaks: %d of %d (%f%% good)" \
			     [llength $gsq] [llength $sequential] $goodsqpct]
	lappend remList [format "Number of good short range peaks: %d of %d (%f%% good)" \
			     [llength $gsr] [llength $shortRange] $goodsrpct]
	lappend remList [format "Number of good long range peaks: %d of %d (%f%% good)" \
			     [llength $glr] [llength $longRange] $goodlrpct]

	lappend remList [format "Overall, number of good peaks:  %d of %d peaks with assignments (%f%% good)" \
			     [llength $goodOverall] [llength $overall] $goodOverallPct]

	lappend remList [format "Fraction of starting long-range forces that are bad (among peakAssignments with previousLikelihoods in range %f  %f ): %f" \
			     $minLikelihoodCutoff $maxLikelihoodCutoff $badInfo]
	
    }
    
    if {$remVar != ""} {
	upvar $remVar r
	foreach elem $remList {
	    lappend r $elem
	}
    }
}


proc removeIntraresiduePeaks args {

    set pot      [requiredFlagVal $args -pot]
    set peakList [flagVal $args -peakList [$pot peaks]]

    set nRemoved 0
    
    foreach p $peakList {
	
	Peak -this $p
	
	if {[$p isIntraresidue]} {
	    $pot removePeakNamed [$p name]
	    $p -acquire
	    incr nRemoved
	}

	rename $p ""
    }
    
    puts stdout "Removed $nRemoved intraresidue peaks\n"
}

proc removeLowLikelihoodPeakAssigns args {

    set pot                 [flagVal $args -pot ""]
    set allPeaks            [flagVal $args -peakList [list]]
    set likelihoodCutoff    [flagVal $args -likelihoodCutoff 0.9]
    
    
    set numPeakAssignsRemoved 0
    set numPeaksUnassignedAtStart 0
    set numPeaksUnassignedAtEnd 0

    if {$pot != ""} {
	set allPeaks [$pot peaks]
    } 

    foreach curPeak $allPeaks {

	Peak -this $curPeak

	if {! [$curPeak isAssigned]} {
	    incr numPeaksUnassignedAtStart
	}

	foreach curPA [$curPeak peakAssignments] {

	    PeakAssignment -this $curPA
	    
	    if {[$curPA previousLikelihood] < $likelihoodCutoff} {
		
		incr numPeakAssignsRemoved

		$curPeak appendToNote [format "Removed peak assignment %s with likelihood %f" \
					   [$curPA name] [$curPA previousLikelihood]]
				
		$curPeak removePeakAssignmentNamed [$curPA name]
	    }
	    
	    rename $curPA ""
	}

	if {! [$curPeak isAssigned]} {
	    incr numPeaksUnassignedAtEnd
	}

	rename $curPeak ""
    }

    updateUser [format "Removed %d peak assignments (unassigning %d whole peaks) because their likelihoods were < %f." \
		    $numPeakAssignsRemoved \
		    [expr $numPeaksUnassignedAtEnd - $numPeaksUnassignedAtStart] \
		    $likelihoodCutoff]
}


proc activateNonviolatedPeakAssigns args {

    #
    # Given a distance violation cutoff, 
    # activate all peakAssigns that are not currently violated
    # and deactivate all peakAssigns that are
    #
    
    set pot                 [flagVal $args -pot ""]
    set peakList            [flagVal $args -peakList [list]]
    set violCutoff          [flagVal $args -violCutoff 0.5]
    
    if {[llength $peakList] == 0} {
	set peakList [$pot peaks]
    }

    foreach curPeak $peakList {
	
	Peak -this $curPeak
	set upBound  [$curPeak upBound]
	set lowBound [$curPeak lowBound]
	
	foreach curPA [$curPeak peakAssignments] {

	    PeakAssignment -this $curPA
	    set curViol [$curPA violation $upBound $lowBound]
	    
	    if {$curViol < $violCutoff} {
		$curPA activate
	    } else {
		$curPA inactivate
	    }
	    
	    rename $curPA ""
	}
	
	rename $curPeak ""
    }
}

proc numPAsInFailedFilterRange args {

    set pot              [requiredFlagVal $args -pot]
    set minFiltersFailed [requiredFlagVal $args -minFiltersFailed]
    set maxFiltersFailed [requiredFlagVal $args -maxFiltersFailed]
    set lrOnly           [flagExists $args -longRangeOnly]

    set retVal 0

    foreach curPeak [$pot peaks] {
	Peak -this $curPeak
	foreach curPA [$curPeak peakAssignments] {
	    PeakAssignment -this $curPA

	    if {(! $lrOnly) || [$curPA isLongRange 6]} {
		
		if {([$curPA numFiltersFailed] >= $minFiltersFailed) &&
		    ([$curPA numFiltersFailed] <= $maxFiltersFailed)} {

		    incr retVal
		}
	    }
	    
	    rename $curPA ""
	}
	rename $curPeak ""
    }

    return $retVal
}

proc numPAsInLikelihoodRange args {

    set pot           [requiredFlagVal $args -pot]
    set minLikelihood [requiredFlagVal $args -minLikelihood]
    set maxLikelihood [requiredFlagVal $args -maxLikelihood]
    set lrOnly        [flagExists $args -longRangeOnly]

    set retVal 0

    foreach curPeak [$pot peaks] {
	Peak -this $curPeak
	foreach curPA [$curPeak peakAssignments] {
	    PeakAssignment -this $curPA

	    if {(! $lrOnly) || [$curPA isLongRange 6]} {
		
		if {([$curPA previousLikelihood] >= $minLikelihood) &&
		    ([$curPA previousLikelihood] <= $maxLikelihood)} {

		    incr retVal
		}
	    }
	    
	    rename $curPA ""
	}
	rename $curPeak ""
    }

    return $retVal
}


proc residueByResiduePeakReport args {

    set potList          [requiredFlagVal $args -pot]
    set likelihoodCutoff [flagVal $args -minLikelihood 0.9]
    set longRangeCutoff  [flagVal $args -longRangeCutoff 6]
    set verb             [flagExists $args [list -verb -verbose]]

    #
    # for each peakAssignment with high-enough likelihood, 
    # record the residues it's involved in 
    #

    foreach pot $potList {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak

	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA previousLikelihood] < $likelihoodCutoff} {
		    rename $curPA ""
		    continue
		}
		
		if {[$curPA isIntraresidue]} {
		    set curPAClass "intra"
		} elseif {[$curPA isSequential]} {
		    set curPAClass "seq"
		} elseif {[$curPA isShortRange $longRangeCutoff]} {
		    set curPAClass "SR"
		} else {
		    set curPAClass "LR"
		}
		
		set fromSA [$curPA fromAssignment]
		set toSA   [$curPA toAssignment]
		
		ShiftAssignment -this $fromSA
		ShiftAssignment -this $toSA
		
		set fromSel [$fromSA protonSelection]
		set toSel   [$toSA   protonSelection]
		
		AtomSel -this $fromSel
		AtomSel -this $toSel
		
		set fromResidues [residuesInSelection $fromSel]
		set toResidues   [residuesInSelection $toSel]
		
		#
		# if the lists of from or to residues are only one element long, 
		# make a copy of it for printing purposes that's not a list (to avoid extraneous braces in output)
		#
		
		if {[llength $fromResidues] == 1} {
		    set cleanFromResidues [lindex $fromResidues 0]
		} else {
		    set cleanFromResidues $fromResidues
		}
		
		if {[llength $toResidues] == 1} {
		    set cleanToResidues [lindex $toResidues 0]
		} else {
		    set cleanToResidues $toResidues
		}
		
		foreach res $fromResidues {
		    lappend residueData($res) [list $curPAClass [join [list [$curPeak name] [$curPA name]]] [join [list [$fromSel string] [$toSel string]]] $cleanToResidues]
		}
		
		foreach res $toResidues {
		    lappend residueData($res) [list $curPAClass [join [list [$curPeak name] [$curPA name]]] [join [list [$fromSel string] [$toSel string]]] $cleanFromResidues]
		}
		
		rename $fromSel ""
		rename $toSel ""
		rename $fromSA ""
		rename $toSA ""
		rename $curPA ""
	    }
	    rename $curPeak ""
	}
    }

    set potNames ""
    foreach pot $potList {
	set potNames [format "%s %s" $potNames [$pot instanceName]]
    }

    set retVal [format "Residue-by-residue peak report for potential(s) %s\n" $potNames]

    appendLineToString retVal "Only peakAssignments with likelihood > $likelihoodCutoff are included\n"

    if {$verb} {
	appendLineToString retVal "For each residue, its specific NOE connectivities are printed with the name of the connected residue, the peak's name, the peakAssignment's name, and the peakAssignment's selections\n"
    }

    foreach curRes [lsort -dictionary [array names residueData]] {

	set curIntra [list]
	set curSeq   [list]
	set curSR    [list]
	set curLR    [list]

	foreach peakRecord $residueData($curRes) {

	    switch [lindex $peakRecord 0] {

		"intra" {lappend curIntra $peakRecord}
		"seq"   {lappend curSeq   $peakRecord}
		"SR"    {lappend curSR    $peakRecord}
		"LR"    {lappend curLR    $peakRecord}
	    }
	}

	# 
	# eliminate entries with duplicate peak/peakAssign data (usually only intraresidue)
	#

	set curIntra [lsort -unique -dictionary -index 1 $curIntra]
	set curSeq   [lsort -unique -dictionary -index 1 $curSeq]
	set curSR    [lsort -unique -dictionary -index 1 $curSR]
	set curLR    [lsort -unique -dictionary -index 1 $curLR]
    
	appendLineToString retVal [format "Residue %s :  %d intraresidue,  %d sequential, %d short range, %d long range" \
				       $curRes [llength $curIntra] [llength $curSeq] [llength $curSR] [llength $curLR]]

	if {$verb} {

	    #
	    # sort entries by other-residue-name for printing
	    #

	    set curSeq   [lsort -dictionary -index 3 $curSeq]
	    set curSR    [lsort -dictionary -index 3 $curSR]
	    set curLR    [lsort -dictionary -index 3 $curLR]
	    
	    appendLineToString retVal "\n   intraresidue connectivities: "
	    if {[llength $curIntra] == 0} {
		appendLineToString retVal "      none"
	    }

	    foreach peakRecord $curIntra {
		appendLineToString retVal [format "      %12s   %-20s   %s" [lindex $peakRecord 3] [lindex $peakRecord 1] [lindex $peakRecord 2]]
	    }

	    appendLineToString retVal "\n   sequential connectivities: "
	    if {[llength $curSeq] == 0} {
		appendLineToString retVal "      none"
	    }

	    foreach peakRecord $curSeq {
		appendLineToString retVal [format "      %12s   %-20s   %s" [lindex $peakRecord 3] [lindex $peakRecord 1] [lindex $peakRecord 2]]
	    }

	    appendLineToString retVal "\n   short range connectivities: "
	    if {[llength $curSR] == 0} {
		appendLineToString retVal "      none"
	    }

	    foreach peakRecord $curSR {
		appendLineToString retVal [format "      %12s   %-20s   %s" [lindex $peakRecord 3] [lindex $peakRecord 1] [lindex $peakRecord 2]]
	    }

	    appendLineToString retVal "\n   long range connectivities: "
	    if {[llength $curLR] == 0} {
		appendLineToString retVal "      none"
	    }

	    foreach peakRecord $curLR {
		appendLineToString retVal [format "      %12s   %-20s   %s" [lindex $peakRecord 3] [lindex $peakRecord 1] [lindex $peakRecord 2]]
	    }

	    appendLineToString retVal [format "\n"]
	}
    }

    return $retVal
}


#
# a quick report on how things are going during the process of filtering
#

proc currentFilterSummary {pot} {

    set maxFiltersFailed 0

    set gsr [list]
    set glr [list]
    set bsr [list]
    set blr [list]
    foreach curPeak [$pot peaks] {
	Peak -this $curPeak
	foreach curPA [$curPeak peakAssignments] {
	    PeakAssignment -this $curPA

	    set maxFiltersFailed [max [list $maxFiltersFailed [$curPA numFiltersFailed]]]

	    if {[$curPA isGood]} {
		if {[$curPA isLongRange 6]} {
		    lappend glr [$curPA numFiltersFailed]
		} else {
		    lappend gsr [$curPA numFiltersFailed]
		}
	    } else {
		if {[$curPA isLongRange 6]} {
		    lappend blr [$curPA numFiltersFailed]
		} else {
		    lappend bsr [$curPA numFiltersFailed]
		}
	    }
	    rename $curPA ""
	}
	rename $curPeak ""
    }
    
    puts [histogram -data [concat $gsr $glr] -min 0 -max $maxFiltersFailed -integerBins -title "Overall good PAs' numFailedFilters"]
    puts [histogram -data $gsr -min 0 -max $maxFiltersFailed -integerBins -title "Short range good PAs' numFailedFilters"]
    puts [histogram -data $glr -min 0 -max $maxFiltersFailed -integerBins -title "Long range good PAs' numFailedFilters"]
    
    puts [histogram -data [concat $bsr $blr] -min 0 -max $maxFiltersFailed -integerBins -title "Overall bad PAs' numFailedFilters"]
    puts [histogram -data $bsr -min 0 -max $maxFiltersFailed -integerBins -title "Short range bad PAs' numFailedFilters"]
    puts [histogram -data $blr -min 0 -max $maxFiltersFailed -integerBins -title "Long range bad PAs' numFailedFilters"]
    
    for {set nFF 0} {$nFF <= $maxFiltersFailed} {incr nFF} { 
	puts [format "Num tot PAs w/ %d filters failed is %d" $nFF [numPAsInFailedFilterRange -pot $pot -minFiltersFailed $nFF -maxFiltersFailed $nFF]]
	puts [format "Num  LR PAs w/ %d filters failed is %d" $nFF  [numPAsInFailedFilterRange -pot $pot -minFiltersFailed $nFF -maxFiltersFailed $nFF -longRangeOnly]]
	puts [format "fbad-long w/ just those PAs active: %f"  [fracBadLRInfo -peakList [$pot peaks] \
								    -minFiltersFailed $nFF \
								    -maxFiltersFailed $nFF \
								    -useGoodBadFlag]]
    }
}





proc newSummarizeMarvinResults args {

    set fnames                           [requiredFlagVal $args -fileNames]
    set potList                          [requiredFlagVal $args [list -pot -potList]]
    set newPeakFilenames                 [flagVal $args [list -newPeakFile -newPeakFiles]]
    set newSAfilenames                   [flagVal $args [list -newShiftAssignmentFile -newShiftAssignmentFiles]]
    set avgPDBfile                       [flagVal $args -avgPDBfile]
    set referenceStructFilename          [flagVal $args -referenceStructureFile ""]
    set violCutoff                       [flagVal $args -violCutoff 0.5]
    set fracConverged                    [flagVal $args -fracConverged 0.1]
    set highLikelihoodCutoff             [flagVal $args -highLikelihoodCutoff 0.9]
    set noeCompletenessCutoff            [flagVal $args -noeCompletenessCutoff 0.0]
    set scatterCutoff                    [flagVal $args -scatterCutoff 999.9]
    set completenessWeight               [flagVal $args -completenessWeight 0.0]
    set viewScriptFileName               [flagVal $args -viewScriptFileName]
    set viewScriptPsfFileName            [flagVal $args -viewScriptPsfFileName]
    set useSingleAssignmentBehavior      [flagExists $args -useSingleAssignmentBehavior]
    set numTalos                         [flagVal $args -numTalosRestraints 0]

    # change default definedBackbone to include DNA/RNA atoms?  Which atoms?

    set backboneSelection  [flagVal $args -definedBackbone   [AtomSel -args "(name ca or name c or name n)"]]
    set heavyatomSelection [flagVal $args -definedHeavyatoms [AtomSel -args "(not name h*)"]]


    #
    # read the PDB files to be analyzed
    #

    updateUser [format "Reading all PDB files\n"]

    set allFiles [eval grabPDBFiles $args]

    #
    # count the of non-violated peaks in each structure
    #

    set allPeaks [list]
    foreach curPot $potList {
	set allPeaks [concat $allPeaks [$curPot peaks]]
    }

    #
    # select the converged structures.
    #

    set goodFiles [chooseBestFraction -files $allFiles \
		       -peakList $allPeaks \
		       -violCutoff $violCutoff \
		       -completenessPots $potList \
		       -inverseBound 4.0 \
		       -inverseMethylCorrection 0.0 \
		       -completenessWeight $completenessWeight \
		       -useNonviolatedCount \
		       -numTalosRestraints $numTalos \
		       -frac $fracConverged]

    #
    # create symbolic links to them, for convenience
    #

    set temp -1
    foreach file $goodFiles {
	set oldname [lindex $file 0]
	set p  [string last . $oldname]
	set newname [string replace $oldname $p end [format "_converged_%d.pdb" [incr temp]]] 
	if {! [file exists $newname]} {
	    exec ln -s $oldname $newname
	}
    }

    #
    # record their num nonviol peaks and NOE completenesses
    #

    set numNonViolPeaks [list]
    
    foreach file $allFiles {
	
	xplorSim setAtomPosArr [lindex $file 1]
	
	lappend numNonViolPeaks [numNonviolatedPeaks -peakList $allPeaks -violCutoff $violCutoff]
    }

    set goodFilesNumNonViolPeaks [list]

    foreach file $goodFiles {

    	xplorSim setAtomPosArr [lindex $file 1]
	
	lappend goodFilesNumNonViolPeaks [numNonviolatedPeaks -peakList $allPeaks -violCutoff $violCutoff]
    }


    #
    # create the view script for vmd-xplor
    #
    
    if {$viewScriptFileName != ""} {

	set temp [list]
	lappend temp $avgPDBfile
	
	if {$referenceStructFilename != ""} {
	    lappend temp $referenceStructFilename
	}

	foreach elem $goodFiles {
	    lappend temp [lindex $elem 0]
	}

	createViewScript \
	    -scriptFileName $viewScriptFileName \
	    -psfFileName $viewScriptPsfFileName \
	    -pdbFileNames $temp \
	    -fitToFileName $avgPDBfile \
	    -selection $backboneSelection
    }

    #
    # report statistics of peak violations used to select converged structs
    #

    set ensemblePeakRemarks [list]

    lappend ensemblePeakRemarks [format "Chose best %d of %d total structs, based on num nonviolated active longrange peaks and overall NOE completeness (completeness weight = %f) " \
				     [llength $goodFiles] [llength $allFiles] $completenessWeight]
    
    lappend ensemblePeakRemarks [histogram -data $numNonViolPeaks     -title "Num nonviolated peaks in all spectra in all structures"]    
    lappend ensemblePeakRemarks [histogram -data $goodFilesNumNonViolPeaks -title "Num nonviolated peaks in all spectra in converged structures"]

    set curRem "Detailed violation and completeness report: "
    foreach f $goodFiles np $goodFilesNumNonViolPeaks {
	appendLineToString curRem [format "%s has %d nonviolated peaks" \
				       [lindex $f 0] $np]
    }

	    
    updateUser [format "Calculating new likelihoods for all peakAssignments and shiftAssignments for each spectrum\n"]

    foreach curPot $potList newPeakFilename $newPeakFilenames newSAfilename $newSAfilenames {

	set peakRemarks [list]

	#
	# calculate new peakAssignment likelihoods
	#
	
	determinePeakAssignmentLikelihoodsFromConvergedStructs \
	    -files $goodFiles \
	    -peakList [$curPot peaks] \
	    -violCutoff $violCutoff 

	#
	# report statistics on new peakAssignment likelihoods
	#
	
	set peakRemarks [concat $ensemblePeakRemarks $peakRemarks]
	
	set peakRemarks [concat $peakRemarks [reportNOEprecision -pot $curPot -definedHeavyatoms $heavyatomSelection -highLikelihoodCutoff $highLikelihoodCutoff]]
	
	if {$referenceStructFilename != ""} {
	    set peakRemarks [concat $peakRemarks [reportNOEaccuracy -pot $curPot -referenceStructureFilename $referenceStructFilename -violCutoff $violCutoff -highLikelihoodCutoff $highLikelihoodCutoff]]
	}

	#
	# tack in more statistics
	#

	initialPeakAnalysis \
	    -pot $curPot \
	    -referenceStructureFile $referenceStructFilename \
	    -violCutoff $violCutoff \
	    -minLikelihood $highLikelihoodCutoff \
	    -description "pass summary" \
	    -remarksVariableName peakRemarks
	
	#
	# write out the new peak file
	#
	
	if {$newPeakFilename != ""} {
	    writeMarvinPeaks \
		-fileName $newPeakFilename \
		-peakList [$curPot peaks] \
		-remarks $peakRemarks
	}

	#
	# calc new ShiftAssignment likelihoods
	#
	
	set invBound [$curPot inverseBound]
	set invMeth  [$curPot inverseMethylCorrection]

	set saRemarks   [list]
	lappend saRemarks [format "Calculating shiftAssignment likelihoods with noeCompleteness cutoff of %f" \
			       $noeCompletenessCutoff]
	lappend saRemarks [format "And proton peak position scatter cutoff of %f ppm" $scatterCutoff]
	lappend saRemarks [format "And NOE violation cutoff of %f A" $violCutoff]
	lappend saRemarks [format "And inverse bound of %f A and inverse methyl correction of %f A" \
			       $invBound $invMeth]
    

	#
	# calculate new shiftAssignment likelihoods
	#

	determineSALikelihoodsFromConvergedStructs \
	    -files $goodFiles \
	    -pot $curPot \
	    -violCutoff $violCutoff \
	    -noeCompletenessCutoff $noeCompletenessCutoff \
	    -scatterCutoff $scatterCutoff \
	    -inverseBound $invBound \
	    -inverseMethylCorrection $invMeth 

	lappend saRemarks [reportShiftAssignLikelihoods -pot $curPot]

	if {$newSAfilename != ""} {
	    writeShiftAssignments \
		-fileName $newSAfilename \
		-shiftAssignments [$curPot shiftAssignments] \
		-remarks $saRemarks
	}

    }

    #
    # calc minimized mean struct
    #

    updateUser [format "Calculating mean structure\n"]

    XplorCommand "set message off echo off print off end" 

    set pdbRemarks [list]
    
    lappend pdbRemarks [format "Note:  Backbone atoms in defined structure region are\n (%s)" \
			    [$backboneSelection string]]
    lappend pdbRemarks [format "Note:  Heavy atoms in defined structure region are\n (%s)" \
			    [$heavyatomSelection string]]

    meanStruct -files $goodFiles -selection $backboneSelection
    
    set curRem "minimized average PDB file created from these files:"
    foreach curFile $goodFiles {
	set curRem [format "%s\n   %s" $curRem [lindex $curFile 0]]
    }
    
    lappend pdbRemarks $curRem
    
    #
    # check coordinate accuracy, if there's a reference structure
    #
    
    if {$referenceStructFilename != ""} {
	
	set heavyDefinedRMSD [compareToReference -referenceFileName $referenceStructFilename \
				  -selection $heavyatomSelection]
	
	set bbnDefinedRMSD [compareToReference -referenceFileName $referenceStructFilename \
				-selection $backboneSelection]
	
	set line1 [format "Accuracy of defined structure region of current coords vs %s :" \
		       $referenceStructFilename]
	
	set line2 [format "   %f (backbone)   %f (heavy atoms)"\
		       $bbnDefinedRMSD $heavyDefinedRMSD]
	
	lappend pdbRemarks [format "%s\n%s" $line1 $line2]
    }
    

    # check structure precision
    
    lappend pdbRemarks [structPrecision -files $goodFiles -verbose \
			    -backboneSelection $backboneSelection \
			    -heavyatomSelection $heavyatomSelection]
    
    # write mean coords
    
    if {$avgPDBfile != ""} {
	writePDB -fileName $avgPDBfile -remarks $pdbRemarks
    } else {
	puts $pdbRemarks
    }
    
    return ""
}


proc completenessSummarizeMarvinResults args {

    set fnames                           [requiredFlagVal $args -fileNames]
    set potList                          [requiredFlagVal $args [list -pot -potList]]
    set newPeakFilenames                 [flagVal $args [list -newPeakFile -newPeakFiles]]
    set newSAfilenames                   [flagVal $args [list -newShiftAssignmentFile -newShiftAssignmentFiles]]
    set violCutoff                       [flagVal $args -violCutoff 0.5]
    set fracConverged                    [flagVal $args -fracConverged 0.1]


    set allFiles [eval grabPDBFiles $args]

    foreach file $allFiles {
	set coordIndex([lindex $file 0]) [lindex $file 1]
    }

    #
    # for each file, eval its NOE completeness, and record the individual SA counts
    #

    foreach file $allFiles {

	set curFname [lindex $file 0]
	xplorSim setAtomPosArr [lindex $file 1]
	
	updateUser [format "Checking completeness of file %s \r" $curFname]
	
	foreach pot $potList {
	    
	    # 
	    # activate all PAs that are violated by < violCutoff
	    #

	    $pot activateAllAssigns

	    foreach curPeak [$pot peaks] {
		Peak -this $curPeak
		foreach curPA [$curPeak peakAssignments] {
		    PeakAssignment -this $curPA
		    set curViol [$curPA violation [$curPeak upBound] [$curPeak lowBound]]
		    if {$curViol > $violCutoff} {
			$curPA inactivate
		    }
		    rename $curPA ""
		}
		rename $curPeak ""
	    }

	    $pot updateNoeCompleteness
	    
	    foreach curSA [$pot shiftAssignments] {
		ShiftAssignment -this $curSA
		
		lappend nAcc([$curSA name]) [list $curFname [$curSA numAccountedNeighbors]]
		lappend nUnacc([$curSA name]) [list $curFname [$curSA numUnaccountedNeighbors]]
		
		rename $curSA ""
	    }
	}
    }
    

    #
    # for each PA, find the best nConverged structures as measured by local completeness around its SAs
    # and record their distances
    #

    set nConverged [expr int(ceil($fracConverged * [llength $allFiles]) - 1)]

    foreach pot $potList {
	
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak 
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		updateUser [format "Checking %s %s %s \r" [$pot instanceName] [$curPeak name] [$curPA name]]
		
		set curPAresults [list]
		
		set curFrom [$curPA fromAssignmentName]
		set curTo   [$curPA toAssignmentName]
		
		set fromAccounted $nAcc($curFrom)
		set fromUnaccounted $nUnacc($curFrom)
		set toAccounted $nAcc($curTo)
		set toUnaccounted $nUnacc($curTo)

		#
		# these lists are already in order by filename, so I can just blow through them
		#

		set localCompletenesses [list]
		
		foreach fa $fromAccounted fu $fromUnaccounted ta $toAccounted tu $toUnaccounted {
		    
		    set curFname [lindex $fa 0]
		    set curLocalAccounted [expr [lindex $fa 1] + [lindex $ta 1]]
		    set curLocalUnaccounted [expr [lindex $fu 1] + [lindex $tu 1]]
		    set curLocalCompleteness [expr $curLocalAccounted / double($curLocalAccounted + $curLocalUnaccounted)]
		    lappend localCompletenesses [list $curFname $curLocalCompleteness]
		}
		
		set localCompletenesses [lsort -index 1 -real -decreasing $localCompletenesses]
		set convergedScore [lindex [lindex $localCompletenesses $nConverged] 1]
		
		set convergedStructs [list]
		foreach lc $localCompletenesses {
		    if {[lindex $lc 1] >= $convergedScore} {
			lappend convergedStructs [lindex $lc 0]
		    }
		}
		
		set nOK 0
		foreach fname $convergedStructs {
		    xplorSim setAtomPosArr $coordIndex($fname)
		    set curV [$curPA violation [$curPeak upBound] [$curPeak lowBound]]
		    if {$curV <= $violCutoff} {
			incr nOK
		    }
		}
	    
		$curPA setPreviousLikelihood [expr $nOK / double([llength $convergedStructs])]

		rename $curPA ""
	    }
	    rename $curPeak ""
	}
    }

    
    foreach pot $potList peakName $newPeakFilenames saName $newSAfilenames {
	
	if {($peakName == "") || ($saName == "")} {
	    continue
	}

	writeMarvinPeaks \
	    -fileName $peakName \
	    -peakList [$pot peaks] 
	
	writeShiftAssignments \
	    -fileName $saName \
	    -shiftAssignments [$pot shiftAssignments] 
    }
    
    return ""
}
