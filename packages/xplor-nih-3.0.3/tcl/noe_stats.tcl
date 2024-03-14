package provide noestats 0.1
package require marvin

namespace eval NoeStats {

namespace export newStats

#
# new, activation-aware NOE statistics procs
# To be called by noeStats, below
#

proc activeFracBadLRInfo args {

    set potList [requiredFlagVal $args [list -pot -potList]]

    set totLRActive 0
    set totBadLRActive 0

    foreach pot $potList {
	
	#
	# find the active LR and bLR PAs, and return appropriate fractions
	#
	    	
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set curPeakNumActivePAs 0
	    set curPeakNumActiveLRPAs 0
	    set curPeakNumActiveBadLRPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA isActive]} {
		    
		    incr curPeakNumActivePAs
		    
		    if {[$curPA isLongRange [$pot longRangePrimarySeqCutoff]]} {
			
			incr curPeakNumActiveLRPAs
			
			if {! [$curPA isGood]} {
			    
			    incr curPeakNumActiveBadLRPAs
			}
		    }
		}
		
		rename $curPA ""
	    }
	    
	    if {$curPeakNumActiveLRPAs > 0} {
		    
		set totLRActive [expr $totLRActive + ($curPeakNumActiveLRPAs / double($curPeakNumActivePAs))]
		set totBadLRActive [expr $totBadLRActive + ($curPeakNumActiveBadLRPAs / double($curPeakNumActivePAs))]
	    }
	    
	    rename $curPeak ""
	}
    }
	
    if {$totLRActive > 0} {
	set retVal [expr $totBadLRActive / double($totLRActive)]
    } else {
	set retVal 0
    }

    return $retVal
}


proc activeNumGoodLRPeaks args {

    set potList [requiredFlagVal $args [list -pot -potList]]

    set totNumGLRPeaks 0
    
    foreach pot $potList {
	    
	#
	# find the peaks that only have active LR PAs, and include an active good LR PA, and return their count
	#
	# This definition agrees with older routines that count good long range peaks w/o looking at activation 
	# when the potential uses split assignment behavior.  When using singleAssignmentBehavior, this definition 
	# reduces to gLR peak == peak with ONLY gLR PA active.
	#
	
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set hasActiveGLR 0
	    set hasActiveSR 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		    
		if {[$curPA isActive]} {
		    
		    if {! [$curPA isLongRange [$pot longRangePrimarySeqCutoff]]} {
			set hasActiveSR 1
		    } 
		    
		    if {[$curPA isLongRange [$pot longRangePrimarySeqCutoff]] &&
			[$curPA isGood]} {
			
			set hasActiveGLR 1
		    } 
		}
		
		rename $curPA ""
	    }
	    
	    if {$hasActiveGLR && (! $hasActiveSR)} {
		
		incr totNumGLRPeaks 
	    }
	    
	    rename $curPeak ""
	}
    }

    return $totNumGLRPeaks
}



proc activeNumLRPeaks args {

    set potList [requiredFlagVal $args [list -pot -potList]]


    set totNumLRPeaks 0
    
    foreach pot $potList {
	    
	#
	# find the peaks that only have active LR PAs, and return their count
	#
	    
	
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set hasActiveLR 0
	    set hasActiveSR 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA isActive]} {
		    
		    if {[$curPA isLongRange [$pot longRangePrimarySeqCutoff]]} {
			
			set hasActiveLR 1
		    } else {
			set hasActiveSR 1
		    }
		}
		
		rename $curPA ""
	    }
	    
	    if {$hasActiveLR && (! $hasActiveSR)} {
		
		incr totNumLRPeaks 
	    }
	    
	    rename $curPeak ""
	}
    }

    return $totNumLRPeaks
}


proc activeNumNonviolatedLRPeaks args {

    set potList [requiredFlagVal $args [list -pot -potList]]
    set cutoff  [flagVal $args -cutoff 0.5]

    set totNumNonviolLRPeaks 0
    
    foreach pot $potList {
	    
	#
	# find the peaks that only have active LR PAs, and have at least one active PA with viol <= cutoff, and return their count
	#
	    
	
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set hasActiveSR 0
	    set hasNonviolLR 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA isActive]} {
		    
		    if {[$curPA isLongRange [$pot longRangePrimarySeqCutoff]]} {
			
			set curViol [$curPA violation [$curPeak upBound] [$curPeak lowBound]]
			if {$curViol <= $cutoff} {
			    set hasNonviolLR 1
			}

		    } else {
			set hasActiveSR 1
		    }
		}
		
		rename $curPA ""
	    }
	    
	    if {(! $hasActiveSR) && $hasNonviolLR} {
		
		incr totNumNonviolLRPeaks 
	    }
	    
	    rename $curPeak ""
	}
    }

    return $totNumNonviolLRPeaks
}



proc activeNumPeaks args {

    set potList [requiredFlagVal $args [list -pot -potList]]


    set totNumActivePeaks 0

    foreach pot $potList {
	    
	#
	# find the peaks that have an active PA, and return their count
	#
	
	
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set hasActivePA 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA isActive]} {
		    
		    set hasActivePA 1
		}
		
		rename $curPA ""
	    }
	    
	    if {$hasActivePA} {
		
		incr totNumActivePeaks 
	    }
	    
	    rename $curPeak ""
	}
    }

    return $totNumActivePeaks
}


proc activeNumNonviolatedPeaks args {

    set potList [requiredFlagVal $args [list -pot -potList]]
    set cutoff  [flagVal $args -cutoff 0.5]


    set totNumNonviolPeaks 0

    foreach pot $potList {
	    
	#
	# find the peaks that have an active PA with violation <= cutoff, and return their count
	#
	
	
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set hasNonviolPA 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA isActive]} {
		    
		    set curViol [$curPA violation [$curPeak upBound] [$curPeak lowBound]]

		    if {$curViol <= $cutoff} {
			set hasNonviolPA 1
		    }
		}
		
		rename $curPA ""
	    }
	    
	    if {$hasNonviolPA} {
		
		incr totNumNonviolPeaks 
	    }
	    
	    rename $curPeak ""
	}
    }

    return $totNumNonviolPeaks
}




proc activeFracGoodPeaks args {

    set potList [requiredFlagVal $args [list -pot -potList]]

    set totActiveGoodPeaks 0
    set totActivePeaks 0

    foreach pot $potList {
	    
	#
	# find the peaks that have active PAs, and return the fraction of them that have good active PAs
	#
	# This definition agrees with older reporting procs that don't use activation when pot uses
	# splitAssignmentBehavior.  When pot uses singleAssignmentBehavior, this definition reduces to 
	# good active peak == peak with only good active PA
	#
	
	
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set hasActivePAs 0
	    set hasActiveGoodPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA isActive]} {
		    
		    set hasActivePAs 1
		    
		    if {[$curPA isGood]} {
			
			set hasActiveGoodPAs 1
		    } 
		}
		
		rename $curPA ""
	    }
	    
	    if {$hasActivePAs} {
		incr totActivePeaks
		if {$hasActiveGoodPAs} {
		    incr totActiveGoodPeaks
		}
	    }
	    
	    rename $curPeak ""
	}
    }
	
    if {$totActivePeaks > 0} {
	set retVal [expr $totActiveGoodPeaks / double($totActivePeaks)]
    } else {
	set retVal 0
    }

    return $retVal
}



proc activeFracGoodLRPeaks args {

    set potList [requiredFlagVal $args [list -pot -potList]]

    set totActiveGLRPeaks 0
    set totActiveLRPeaks 0

    foreach pot $potList {
	
	#
	# find the peaks that only have LR PAs active, and return the fraction of them that have a good active PA
	#
	# This definition agrees with non-activation-aware reporting procs for splitAssignmentBehavior, 
	# and reduces to gLR peak == peak with only gLR PAs active for singleAssignmentBehavior
	#
	
	
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set hasActiveSRPAs 0
	    set hasActiveLRPAs 0
	    set hasActiveGLRPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA isActive]} {
		    
		    if {! [$curPA isLongRange [$pot longRangePrimarySeqCutoff]]} {
			
			set hasActiveSRPAs 1
		    } 
		    
		    if {[$curPA isLongRange [$pot longRangePrimarySeqCutoff]]} {
			
			set hasActiveLRPAs 1
			
			if {[$curPA isGood]} {
			    
			    set hasActiveGLRPAs 1 
			} 
		    }
		}
		
		rename $curPA ""
	    }
	    
	    if {$hasActiveLRPAs && (! $hasActiveSRPAs)} {
		incr totActiveLRPeaks
		if {$hasActiveGLRPAs} {
		    incr totActiveGLRPeaks
		}
	    }
	    
	    rename $curPeak ""
	}
    }
	
    if {$totActiveLRPeaks > 0} {
	set retVal [expr $totActiveGLRPeaks / double($totActiveLRPeaks)]
    } else {
	set retVal 0
    }
    
    return $retVal
}


proc activePACompleteness args {

    set potList [requiredFlagVal $args [list -pot -potList]]

    #
    # Determine fraction of close SA pairs that correspond to an active PA
    #
    # Depends on current coordinates!
    #
    # Just uses the pot's code, to ensure that everything is consistent
    #

    set totClose 0
    set totClosePAAccounted 0

    foreach pot $potList {

	$pot updateNoeCompleteness
	
	set totClose [expr $totClose + [$pot numCloseShiftAssignmentPairs]]
	set totClosePAAccounted [expr $totClosePAAccounted + [$pot numClosePAAccountedShiftAssignmentPairs]]

    }

    if {$totClose > 0} {
	set curCompleteness [expr $totClosePAAccounted / double($totClose)]
    } else {
	set curCompleteness 1
    }
	
    return $curCompleteness
}


proc activeOverallCompleteness args {

    set potList [requiredFlagVal $args [list -pot -potList]]

    #
    # Determine fraction of close SA pairs that correspond to an active PA or an exception
    #
    # Depends on current coordinates!
    #
    # Just uses the pot's code, to ensure that everything is consistent
    #

    set totClose 0
    set totCloseAccounted 0

    foreach pot $potList {

	$pot updateNoeCompleteness
	
	set totClose [expr $totClose + [$pot numCloseShiftAssignmentPairs]]
	set totCloseAccounted [expr $totCloseAccounted + [$pot numCloseAccountedShiftAssignmentPairs]]

    }

    if {$totClose > 0} {
	set curCompleteness [expr $totCloseAccounted / double($totClose)]
    } else {
	set curCompleteness 1
    }
	
    return $curCompleteness
}


proc activeFracInverseRestraints args {

    set potList [requiredFlagVal $args [list -pot -potList]]

    #
    # Determine fraction of SA pairs that do not correspond to an active PA or an exception
    #
    
    foreach pot $potList {
	foreach curFrom [$pot fromShiftAssignments] {
	    ShiftAssignment -this $curFrom
	    
	    foreach curTo [$pot toShiftAssignments] {
		ShiftAssignment -this $curTo
		
		set curPairName [join [list [$curFrom name] [$curTo name]]]
		set isAccounted($curPairName) 0
		set isException($curPairName) 0
		
		rename $curTo ""
	    }
	    rename $curFrom ""
	}
	
	
	foreach {curFromName curToName} [$pot explicitInverseExceptions] {
	    set curPairName [join [list $curFromName $curToName]]
	    set isException($curPairName) 1
	}
    }
    

    foreach curPairName [array names isAccounted] {
	set isAccounted($curPairName) 0
    }

    foreach pot $potList {	    
	
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
		
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA isActive]} {
		    
		    set curPairName [join [list [$curPA fromAssignmentName] [$curPA toAssignmentName]]]
		    set isAccounted($curPairName) 1
		}
		rename $curPA ""
	    }
	    rename $curPeak ""
	}
    }

    set totActivePairs 0
    set totPairs 0

    foreach curPairName [array names isAccounted] {

	incr totPairs 
	    
	if {! ($isAccounted($curPairName) || $isException($curPairName))} {
	    incr totActivePairs
	}
    }

    if {$totPairs > 0} {
	set retVal [expr $totActivePairs / double($totPairs)]
    } else {
	set retVal -1
    }

    return $retVal
}



proc newDegeneracy args {

    set potList [requiredFlagVal $args [list -pot -potList]]

    #
    # Return the average number of PAs with non-zero likelihoods per assigned peak.
    # Insensitive to activation state.
    #
	
    set degens [list]
    
    foreach pot $potList {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	
	    set curDegen 0
	
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
	    
		if {[$curPA previousLikelihood] > 0} {
		    incr curDegen
		}
	    
		rename $curPA ""
	    }

	    if {$curDegen > 0} {
		lappend degens $curDegen
	    }

	    rename $curPeak ""
	}
    }
    
    if {[llength $degens] > 0} {
	set retVal [list [mean $degens] [standardDeviation $degens]]
    } else {
	set retVal [list -1 -1]
    }

    return $retVal
}




proc newFracBadLowLikelihoodPeaks args {

    set potList [requiredFlagVal $args [list -pot -potList]]
    set cutoff  [flagVal $args -cutoff 0.1]
        
    #
    # find the assigned peaks that have max prev likelihood <= cutoff  and return the fraction of them that have only bad PAs
    #
    # Insensitive to activation state
    #
	
    set totLowLikelihoodPeaks 0
    set totBadLowLikelihoodPeaks 0
	
    foreach pot $potList {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set nPAs 0
	    set nLowLikelihoodPAs 0
	    set nBadLowLikelihoodPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		incr nPAs
		
		if {[$curPA previousLikelihood] <= $cutoff} {
		    incr nLowLikelihoodPAs
		    
		    if {! [$curPA isGood]} {
			incr nBadLowLikelihoodPAs
		    }
		}
		
		rename $curPA ""
	    }
	    
	    if {($nLowLikelihoodPAs == $nPAs) && ($nPAs > 0)} {
		incr totLowLikelihoodPeaks 

		if {$nBadLowLikelihoodPAs == $nLowLikelihoodPAs} {
		    incr totBadLowLikelihoodPeaks
		}
	    }
	    
	    rename $curPeak ""
	}
    }
	
    if {$totLowLikelihoodPeaks > 0} {
	set retVal [expr $totBadLowLikelihoodPeaks / double($totLowLikelihoodPeaks)]
    } else {
	set retVal -1
    }
    
    return $retVal
}



proc newFracLowLikelihoodPeaksBad args {

    set potList [requiredFlagVal $args [list -pot -potList]]
    set cutoff  [flagVal $args -cutoff 0.1]
        
    #
    # find the assigned peaks that have only bad PAs, and return the fraction of them with max prev likelihood <= cutoff
    #
    # Insensitive to activation state
    #
	
    set totBadPeaks 0
    set totBadLowLikelihoodPeaks 0
	
    foreach pot $potList {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set nPAs 0
	    set nBadPAs 0
	    set nBadLowLikelihoodPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		incr nPAs

		if {! [$curPA isGood]} {

		    incr nBadPAs

		    if {[$curPA previousLikelihood] <= $cutoff} {
			
			incr nBadLowLikelihoodPAs
		    }
		}

		rename $curPA ""
	    }

	    if {($nPAs > 0) && ($nBadPAs == $nPAs)} {
		incr totBadPeaks 
		
		if {$nBadLowLikelihoodPAs == $nBadPAs} {
		    incr totBadLowLikelihoodPeaks
		}
	    }
	
	    rename $curPeak ""
	}
    }
	
    if {$totBadPeaks > 0} {
	set retVal [expr $totBadLowLikelihoodPeaks / double($totBadPeaks)]
    } else {
	set retVal -1
    }

    return $retVal
}



proc newFracGoodHighLikelihoodPeaks args {

    set potList    [requiredFlagVal $args [list -pot -potList]]
    set highCutoff [flagVal $args -highCutoff 0.9]
    set lowCutoff  [flagVal $args -lowCutoff  0.1]

    #
    # find the peaks w/ high-likelihood PAs and no intermediate-likelihood PAs.  
    # Return fraction of them where all high-likelihood PAs are good.
    #
    # Insensitive to activation state
    #
	
    set totHighLikelihoodPeaks 0
    set totGoodHighLikelihoodPeaks 0
	
    foreach pot $potList {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set nHighLikelihoodPAs 0
	    set nGoodHighLikelihoodPAs 0
	    set nIntermediateLikelihoodPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA previousLikelihood] >= $highCutoff} {
		    incr nHighLikelihoodPAs
		    if {[$curPA isGood]} {
			incr nGoodHighLikelihoodPAs
		    }
		} 
		
		if {([$curPA previousLikelihood] > $lowCutoff) &&
		    ([$curPA previousLikelihood] < $highCutoff)} {
		    
		    incr nIntermediateLikelihoodPAs
		}
		
		rename $curPA ""
	    }
	
	    if {($nHighLikelihoodPAs > 0) && ($nIntermediateLikelihoodPAs == 0)} {
		incr totHighLikelihoodPeaks 

		if {$nGoodHighLikelihoodPAs == $nHighLikelihoodPAs} {
		    incr totGoodHighLikelihoodPeaks
		}
	    }
	
	    rename $curPeak ""
	}
    }
    
    if {$totHighLikelihoodPeaks > 0} {
	set retVal [expr $totGoodHighLikelihoodPeaks / double($totHighLikelihoodPeaks)]
    } else {
	set retVal -1
    }

    return $retVal
}



proc newFracHighLikelihoodPeaksGood args {

    set potList    [requiredFlagVal $args [list -pot -potList]]
    set highCutoff [flagVal $args -highCutoff 0.9]
    set lowCutoff  [flagVal $args -lowCutoff  0.1]

    #
    # find the peaks that have good PAs.  Return fraction of them that have no intermediate-likelihood 
    # PAs, 1+ high-likelihood PAs,  and only good high-likelihood PAs
    #
    # Insensitive to activation state
    #
	
    set totGoodPeaks 0
    set totGoodHighLikelihoodPeaks 0
	
    foreach pot $potList {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set nGoodPAs 0
	    set nHighLikelihoodPAs 0
	    set nGoodHighLikelihoodPAs 0
	    set nIntermediateLikelihoodPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA isGood]} {
		    
		    incr nGoodPAs
		    
		    if {[$curPA previousLikelihood] >= $highCutoff} {
			
			incr nGoodHighLikelihoodPAs
		    }
		}
		
		if {[$curPA previousLikelihood] >= $highCutoff} {
		    incr nHighLikelihoodPAs
		} 
		
		if {([$curPA previousLikelihood] > $lowCutoff) &&
		    ([$curPA previousLikelihood] < $highCutoff)} {
		    
		    incr nIntermediateLikelihoodPAs
		}
		
		rename $curPA ""
	    }
	
	    if {$nGoodPAs > 0} {
		incr totGoodPeaks 

		if {($nHighLikelihoodPAs > 0) && 
		    ($nGoodHighLikelihoodPAs == $nHighLikelihoodPAs) &&
		    ($nIntermediateLikelihoodPAs == 0) } {
		    
		    incr totGoodHighLikelihoodPeaks
		}
	    }
	
	    rename $curPeak ""
	}
    }
    
    if {$totGoodPeaks > 0} {
	set retVal [expr $totGoodHighLikelihoodPeaks / double($totGoodPeaks)]
    } else {
	set retVal -1
    }

    return $retVal
}


proc numHighLikelihoodPeaks args {

    set potList    [requiredFlagVal $args [list -pot -potList]]
    set highCutoff [flagVal $args -highCutoff 0.9]
    set lowCutoff  [flagVal $args -lowCutoff  0.1]

    #
    # find the peaks w/ high-likelihood PAs and no intermediate-likelihood PAs.  
    # Return their count.
    #
    # Insensitive to activation state
    #
	
    set totHighLikelihoodPeaks 0
	
    foreach pot $potList {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set nHighLikelihoodPAs 0
	    set nProblemPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA previousLikelihood] >= $highCutoff} {
		    incr nHighLikelihoodPAs
		} 
		
		if {([$curPA previousLikelihood] > $lowCutoff) &&
		    ([$curPA previousLikelihood] < $highCutoff)} {
		    
		    incr nProblemPAs
		}
		
		rename $curPA ""
	    }
	
	    if {($nHighLikelihoodPAs > 0) && ($nProblemPAs == 0)} {
		incr totHighLikelihoodPeaks 
	    }
	
	    rename $curPeak ""
	}
    }
    
    return $totHighLikelihoodPeaks
}


proc numGoodHighLikelihoodPeaks args {

    set potList    [requiredFlagVal $args [list -pot -potList]]
    set highCutoff [flagVal $args -highCutoff 0.9]
    set lowCutoff  [flagVal $args -lowCutoff  0.1]

    #
    # find the peaks w/ only good high-likelihood PAs and no intermediate-likelihood PAs.  
    # Return their count.
    #
    # Insensitive to activation state
    #
	
    set totHighLikelihoodPeaks 0
	
    foreach pot $potList {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set nGoodHighLikelihoodPAs 0
	    set nProblemPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA previousLikelihood] >= $highCutoff} {
		    if {[$curPA isGood]} {
			incr nGoodHighLikelihoodPAs
		    } else {
			incr nProblemPAs
		    }
		} 
		
		if {([$curPA previousLikelihood] > $lowCutoff) &&
		    ([$curPA previousLikelihood] < $highCutoff)} {
		    
		    incr nProblemPAs
		}
		
		rename $curPA ""
	    }
	
	    if {($nGoodHighLikelihoodPAs > 0) && ($nProblemPAs == 0)} {

		incr totHighLikelihoodPeaks 
	    }
	
	    rename $curPeak ""
	}
    }
    
    return $totHighLikelihoodPeaks
}


proc numHighLikelihoodLRPeaks args {

    set potList    [requiredFlagVal $args [list -pot -potList]]
    set highCutoff [flagVal $args -highCutoff 0.9]
    set lowCutoff  [flagVal $args -lowCutoff  0.1]

    #
    # find the peaks where all PAs are either high-likelihood LR, or low likelihood.
    # Return their count.
    #
    # Insensitive to activation state
    #
	
    set totHighLikelihoodLRPeaks 0
	
    foreach pot $potList {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set nHighLikelihoodLRPAs 0
	    set nProblemPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA previousLikelihood] >= $highCutoff} {

		    if {[$curPA isLongRange [$pot longRangePrimarySeqCutoff]]} {
			incr nHighLikelihoodLRPAs
		    } else {
			incr nProblemPAs
		    }
		}

		if {([$curPA previousLikelihood] > $lowCutoff) &&
		    ([$curPA previousLikelihood] < $highCutoff)} {
		    
		    incr nProblemPAs
		}
		
		rename $curPA ""
	    }
	
	    if {($nHighLikelihoodLRPAs > 0) && ($nProblemPAs == 0)} {
		incr totHighLikelihoodLRPeaks 
	    }
	
	    rename $curPeak ""
	}
    }
    
    return $totHighLikelihoodLRPeaks
}


proc numGoodHighLikelihoodLRPeaks args {

    set potList    [requiredFlagVal $args [list -pot -potList]]
    set highCutoff [flagVal $args -highCutoff 0.9]
    set lowCutoff  [flagVal $args -lowCutoff  0.1]

    #
    # find the peaks where all PAs are either high-likelihood LR, or low likelihood.
    # Return their count.
    #
    # Insensitive to activation state
    #
	
    set totHighLikelihoodLRPeaks 0
	
    foreach pot $potList {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set nGoodHighLikelihoodLRPAs 0
	    set nProblemPAs 0
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		if {[$curPA previousLikelihood] >= $highCutoff} {

		    if {[$curPA isGood] && \
			    [$curPA isLongRange \
				 [$pot longRangePrimarySeqCutoff]]} {
			incr nGoodHighLikelihoodLRPAs
		    } else {
			incr nProblemPAs
		    }
		}

		if {([$curPA previousLikelihood] > $lowCutoff) &&
		    ([$curPA previousLikelihood] < $highCutoff)} {
		    
		    incr nProblemPAs
		}
		
		rename $curPA ""
	    }
	
	    if {($nGoodHighLikelihoodLRPAs > 0) && ($nProblemPAs == 0)} {
		incr totHighLikelihoodLRPeaks 
	    }
	
	    rename $curPeak ""
	}
    }
    
    return $totHighLikelihoodLRPeaks
}


   

proc newDiscrimination args {

    set potList    [requiredFlagVal $args [list -pot -potList]]
    set highCutoff [flagVal $args -highCutoff 0.9]
    set lowCutoff  [flagVal $args -lowCutoff 0.1]

    
    #
    # return the fraction of all assigned peaks that do not have have PAs with likelihoods > lowCutoff and < highCutoff
    #
    # Insensitive to activation state
    #

    set totAssignedPeaks 0
    set totExtremeLikelihoodPeaks 0
	
    foreach pot $potList {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    set nPAs 0
	    set nIntermediateLikelihoodPAs 0

	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		incr nPAs
		if {([$curPA previousLikelihood] > $lowCutoff) &&
		    ([$curPA previousLikelihood] < $highCutoff)} {
		    
		    incr nIntermediateLikelihoodPAs
		}
		rename $curPA ""
	    }

	    if {$nPAs > 0} {

		incr totAssignedPeaks 

		if {$nIntermediateLikelihoodPAs == 0} {
		    incr totExtremeLikelihoodPeaks
		}
	    }

	    rename $curPeak ""
	}
    }

    if {$totAssignedPeaks > 0} {
	set retVal [expr $totExtremeLikelihoodPeaks / double($totAssignedPeaks)]
    } else {
	set retVal -1
    }
    
    return $retVal
}


proc newFracGoodExceptions args {
    
    set potList [requiredFlagVal $args [list -pot -potList]]

    #
    # Return fraction of exceptions that correspond to close SA pairs
    #
    # Depends on current coordinates!
    #

    set totGoodExceptions 0
    set totExceptions 0

    foreach pot $potList {
	set cutoff [$pot inverseBound]
	foreach {curFromName curToName} [$pot explicitInverseExceptions] {
	    set curFrom [$pot shiftAssignmentNamed $curFromName]
	    set curTo   [$pot shiftAssignmentNamed $curToName]
	    ShiftAssignment -this $curFrom
	    
	    if {[$curFrom distanceToShiftAssignment $curTo] <= $cutoff} {
		incr totGoodExceptions
	    }
	    
	    incr totExceptions
	    rename $curFrom ""
	}
    }

    if {$totExceptions > 0} {
	set retVal [expr $totGoodExceptions / double($totExceptions)]
    } else {
	set retVal -1
    }
    
    return $retVal
}


proc numSAsWithNoHighLikelihoodPAs args {

    set potList    [requiredFlagVal $args [list -pot -potList]]
    set highCutoff [flagVal $args -highCutoff 0.9]

    array unset nHighLikelihoodPAs

    foreach pot $potList {

	foreach curSA [$pot shiftAssignments] {
	    ShiftAssignment -this $curSA
	    set nHighLikelihoodPAs([$curSA name]) 0
	    rename $curSA ""
	}

	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		if {[$curPA previousLikelihood] >= $highCutoff} {
		    incr nHighLikelihoodPAs([$curPA fromAssignmentName])
		    incr nHighLikelihoodPAs([$curPA toAssignmentName])
		}
		rename $curPA ""
	    }
	    rename $curPeak ""
	}
    }

    set retVal 0

    foreach saName [array names nHighLikelihoodPAs] {
	if {$nHighLikelihoodPAs($saName) == 0} {
	    incr retVal
	}
    }

    return $retVal
}




proc newStats args {

    set potList [requiredFlagVal $args [list -pot -potList]]
    set label   [flagVal $args -label " "]
    set nTimes  [flagVal $args -nTimes 20]


    set retVal [format "\nNew NOE stats report of pot(s)"]
    foreach pot $potList {

	if {[$pot usesSingleAssignmentBehavior]} {
	    set behav "single"
	} else {
	    set behav "split"
	}

	set retVal [format "%s %s using %s assignment behavior" $retVal [$pot instanceName] $behav]
    }

    set retVal [format "%s at stage %s : " $retVal $label]

    #
    # since updateActivation takes a while, I'm moving the loop to here
    #

    set allDat [list]

    for {set count 0} {$count < $nTimes} {incr count} {

	foreach pot $potList {

	    #
	    # update activation, using previous likelihoods only
	    #
	    
	    $pot setPreviousLikelihoodWeight 1.0
	    $pot setDistanceViolationWeight  0.0
	    $pot setNoeCompletenessWeight    0.0
	    $pot setScatterWeight            0.0
	    $pot disallowShiftAssignmentInactivation 
	    $pot setInverseBound 4.0
	    
	    $pot cacheNeighbors
	    
	    $pot updateActivation
	}

	set curDat [list]

	updateUser [format "Calculating activation-sensitive NOE statistics, cycle %d of %d \r" $count $nTimes]

	lappend curDat [activeFracBadLRInfo -potList $potList]
	lappend curDat [activeNumGoodLRPeaks -potList $potList]
	lappend curDat [activeNumLRPeaks -potList $potList]
	lappend curDat [activeNumPeaks -potList $potList]
	lappend curDat [activeFracGoodPeaks -potList $potList]
	lappend curDat [activeFracGoodLRPeaks -potList $potList]
	lappend curDat [activePACompleteness -potList $potList]
	lappend curDat [activeOverallCompleteness -potList $potList]
	lappend curDat [activeFracInverseRestraints -potList $potList]	

	lappend allDat $curDat
    }
    
    set curDat [list]
    foreach elem $allDat {
	lappend curDat [lindex $elem 0]
    }
    appendLineToString retVal [format "Frac bad/LR info = %f +/- %f" [mean $curDat] [standardDeviation $curDat]]


    set curDat [list]
    foreach elem $allDat {
	lappend curDat [lindex $elem 1]
    }
    appendLineToString retVal [format "Num good LR peaks = %f +/- %f" [mean $curDat] [standardDeviation $curDat]]

    set curDat [list]
    foreach elem $allDat {
	lappend curDat [lindex $elem 2]
    }
    appendLineToString retVal [format "Num LR peaks = %f +/- %f" [mean $curDat] [standardDeviation $curDat]]

    set curDat [list]
    foreach elem $allDat {
	lappend curDat [lindex $elem 3]
    }
    appendLineToString retVal [format "Num active peaks = %f +/- %f" [mean $curDat] [standardDeviation $curDat]]

    set curDat [list]
    foreach elem $allDat {
	lappend curDat [lindex $elem 4]
    }
    appendLineToString retVal [format "Frac good/active peaks = %f +/- %f" [mean $curDat] [standardDeviation $curDat]]

    set curDat [list]
    foreach elem $allDat {
	lappend curDat [lindex $elem 5]
    }
    appendLineToString retVal [format "Frac good/LR peaks = %f +/- %f" [mean $curDat] [standardDeviation $curDat]]

    set curDat [list]
    foreach elem $allDat {
	lappend curDat [lindex $elem 6]
    }
    appendLineToString retVal [format "PA completeness =  %f +/- %f" [mean $curDat] [standardDeviation $curDat]]

    set curDat [list]
    foreach elem $allDat {
	lappend curDat [lindex $elem 7]
    }
    appendLineToString retVal [format "Overall completeness =  %f +/- %f" [mean $curDat] [standardDeviation $curDat]]

    set curDat [list]
    foreach elem $allDat {
	lappend curDat [lindex $elem 8]
    }
    appendLineToString retVal [format "Frac inverse restraints active =  %f +/- %f" [mean $curDat] [standardDeviation $curDat]]

    #
    # Now report statistics that are insensitive to activation state
    #

    set temp [newDegeneracy -potList $potList]
    appendLineToString retVal [format "NOE degeneracy = %f +/- %f" [lindex $temp 0] [lindex $temp 1]]
    appendLineToString retVal [format "Frac bad/low likelihood peaks = %f" [newFracBadLowLikelihoodPeaks -potList $potList]]
    appendLineToString retVal [format "Frac low likelihood/ bad peaks = %f" [newFracLowLikelihoodPeaksBad -potList $potList]]
    appendLineToString retVal [format "Frac good/high likelihood peaks = %f" [newFracGoodHighLikelihoodPeaks -potList $potList]]
    appendLineToString retVal [format "Frac high likelihood/ good peaks = %f" [newFracHighLikelihoodPeaksGood -potList $potList]]
    appendLineToString retVal [format "NOE discrimination = %f" [newDiscrimination -potList $potList]]
    appendLineToString retVal [format "Frac good/exceptions = %f" [newFracGoodExceptions -potList $potList]]
    appendLineToString retVal [format "Num high likelihood peaks = %d" [numHighLikelihoodPeaks -potList $potList]]
    appendLineToString retVal [format "Num good high likelihood peaks = %d" [numGoodHighLikelihoodPeaks -potList $potList]]
    appendLineToString retVal [format "Num high likelihood LR peaks = %d" [numHighLikelihoodLRPeaks -potList $potList]]
    appendLineToString retVal [format "Num good high likelihood LR peaks = %d" [numGoodHighLikelihoodLRPeaks -potList $potList]]
    appendLineToString retVal [format "Num SAs with no high likelihood PAs = %d" [numSAsWithNoHighLikelihoodPAs -potList $potList]]

    return $retVal
}








}

namespace import NoeStats::*
