package provide marvin 1.0

package require pasdpot
package require marvinassignment


proc symmetryFilter args {

    #
    # Given a set of shiftAssignments and peaks/peakAssignments,
    # predict (using the shiftAssignments) whether it's possible to see 
    # a peak with a symmetric peakAssignment.  If it's possible, check each 
    # peak for a match.  Then evaluate the actual matching peak(s)
    # unfolded peak position relative to the original peak's to ensure
    # that they're within the given shift tolerance of each other.  If 
    # a match is possible, but no peak has such a symmetric peakAssignment 
    # within the shift tolerance, check to see if that peak's upper distance
    # bound is within the enforcement range.  If so, mark that peakAssignment 
    # as failing the filter. 
    #
    # A special mode for testing purposes:  
    #
    # -noLink ensures that linked peakAssignments (for symmetry pairs) are not made.
    #

    set pot            [requiredFlagVal $args -pot]
    set shiftAssigns   [flagVal $args -shiftAssignments [$pot shiftAssignments]]
    set peaks          [flagVal $args -peaks [$pot peaks]]
    set tolerance      [flagVal $args -tolerance 0.03]
    set maxUpperBound  [flagVal $args -maxUpperBound 99999.9]
    set remVar         [flagVal $args -remarksVariableName ""]
    set noLinkMode     [flagExists $args -noLink]

    #
    # for each shiftAssignment, record its to-from partner name
    #

    foreach curAssign $shiftAssigns {

	ShiftAssignment -this $curAssign
	set curName [$curAssign name]
	
	if {[$curAssign hasToFromPartnerName]} {
	    set toFromPartner($curName) [$curAssign toFromPartnerName]
	} else {
	    set toFromPartner($curName) "NONE"
	}

	rename $curAssign ""
    }
    
    #
    # for each Peak/PeakAssignment
    # record Peak, PeakAssign, from, and to shiftAssign names in a list,
    # with an array index based on the names of its from and to shiftAssigns
    #

    set peakAssData [list]

    set count 0

    foreach p $peaks {

	updateUser [format "Recording data for peak %d of %d \r" \
			[incr count] [llength $peaks]]

	Peak -this $p

	foreach pa [$p peakAssignments] {

	    PeakAssignment -this $pa

	    if {[$pa hasUnfoldedFromProtonPeakPosition]} {
		set fromProtonPP [$pa unfoldedFromProtonPeakPosition]
	    } else {
		set fromProtonPP "NONE"
	    }
	    
	    if {[$pa hasUnfoldedToProtonPeakPosition]} {
		set toProtonPP   [$pa unfoldedToProtonPeakPosition]
	    } else {
		set toProtonPP "NONE"
	    }

	    lappend peakAssData [list [$p name] [$pa name] [$pa fromAssignmentName] [$pa toAssignmentName] \
				    $fromProtonPP $toProtonPP]
	    
	    lappend peakAssIndex([join [list [$pa fromAssignmentName] [$pa toAssignmentName]]]) \
		[expr [llength $peakAssData] - 1]
	    
	    rename $pa ""
	}

	rename $p ""
    }

    #
    # For each PeakAssignment, see if a symmetry match is possible.  
    # If so, see if one exists.  
    # If the expected symmetry match is missing, fail the PeakAssignment.
    #

    set count 0
    set numNoPossSymmPeakAssigns 0
    set numPeakAssignsKeptForHighUpperBound 0

    foreach curData $peakAssData {

	updateUser [format "Examining peak assignment %d of %d for symmetry \r" \
			[incr count] [llength $peakAssData]]

	set curPeakName       [lindex $curData 0]
	set curPAName         [lindex $curData 1]
	set curFromAssignName [lindex $curData 2]
	set curToAssignName   [lindex $curData 3]
	set curFromProtonPP   [lindex $curData 4]
	set curToProtonPP     [lindex $curData 5]

	set partnerFromAssignName $toFromPartner($curToAssignName)
	set partnerToAssignName   $toFromPartner($curFromAssignName)

	set curPeak [$pot peakNamed $curPeakName]
	Peak -this $curPeak

	set curPA [$curPeak peakAssignmentNamed $curPAName]
	PeakAssignment -this $curPA

	# If a symmetry match is impossible, leave the PeakAssignment in place

	if {($partnerFromAssignName == "NONE") || ($partnerToAssignName == "NONE")} {

	    $curPA appendToNote "No possible symmetry" 
	    incr numNoPossSymmPeakAssigns 

	    rename $curPA ""
	    rename $curPeak ""
	    continue
	}

	# Look up matching peak assignments based on the partner shiftAssign 

	set shiftAssignMatchList [list]

	if {[info exists peakAssIndex([join [list $partnerFromAssignName $partnerToAssignName]])]} {
	    set shiftAssignMatchList $peakAssIndex([join [list $partnerFromAssignName $partnerToAssignName]])
	}

	#
	# make sure a peak can't be its own symmetry partner
	#

	set cleanShiftAssignMatchList [list]

	foreach i $shiftAssignMatchList {
	    
	    set matchData           [lindex $peakAssData $i]
	    set matchPeakName       [lindex $matchData 0]

	    if {$matchPeakName != $curPeakName} {
		lappend cleanShiftAssignMatchList $i
	    }
	}

	#
	# check each shift assign match's unfolded peak positions to ensure they're within tolerance
	#

	set peakPosnMatchList [list]

	foreach i $cleanShiftAssignMatchList {

	    set matchData           [lindex $peakAssData $i]
	    set matchPeakName       [lindex $matchData 0]
	    set matchPAName         [lindex $matchData 1]
	    set matchFromAssignName [lindex $matchData 2]
	    set matchToAssignName   [lindex $matchData 3]
	    set matchFromProtonPP   [lindex $matchData 4]
	    set matchToProtonPP     [lindex $matchData 5]

	    # don't enforce peak position matching if folded peak positions aren't known

	    if {($curFromProtonPP == "NONE") ||
		($curToProtonPP == "NONE") ||
		($matchFromProtonPP == "NONE") ||
		($matchToProtonPP == "NONE")} {
		
		lappend peakPosnMatchList $i

	    } else {

		set delta [expr sqrt(pow($curFromProtonPP - $matchToProtonPP, 2) + \
					 pow($curToProtonPP - $matchFromProtonPP, 2))]

		if {$delta <= $tolerance} {
		    lappend peakPosnMatchList $i
		}
	    }
	}

	if {[llength $peakPosnMatchList] == 0} {
	    
	    # If no expected symmetry match is found, either mark the peakAssignment as 
	    # having no symmetry partner, or mark it for deletion, depending on the peak's 
	    # upper bound (since NOE completeness drops with upper bound)

	    if {[$curPeak upBound] <= $maxUpperBound} {
		
		lappend PAsToFail([$curPeak name]) [$curPA name]
	    
	    } else {

		$curPA appendToNote "No symmetry partner found, but assignment kept for high upper bound" 
		incr numPeakAssignsKeptForHighUpperBound
	    }

	} else {

	    #
	    # mark the assignment with the name(s) of all matching symmetry partners
	    #

	    foreach i $peakPosnMatchList {

		set matchData           [lindex $peakAssData $i]
		set matchPeakName       [lindex $matchData 0]
		set matchPAName         [lindex $matchData 1]
		set matchFromAssignName [lindex $matchData 2]
		set matchToAssignName   [lindex $matchData 3]

		if {! $noLinkMode} {
		    $curPA addLinkedPeakAssignmentName $matchPeakName $matchPAName
		}
		$curPA appendToNote "Symmetry partner is peak $matchPeakName assignment $matchPAName" 
	    }
	}
	
	rename $curPA ""
	rename $curPeak ""
    }
	
    #
    # Now record statistics on failed peakAssigns
    #

    set numPeakAssignsFailed 0
    set numGoodPeakAssignsFailed 0
    set numPeaksFailed 0
    set numGoodPeaksFailed 0

    foreach peakName [array names PAsToFail] {
	
	set curPeak [$pot peakNamed $peakName]
	Peak -this $curPeak 
	set curPeakGood [$curPeak isGood]
	
	if {[llength $PAsToFail($peakName)] == [$curPeak numPeakAssignments]} {
	    incr numPeaksFailed 
	    if {$curPeakGood} {
		incr numGoodPeaksFailed
	    }
	}
	
	foreach paName $PAsToFail($peakName) {
		
	    set curPA [$curPeak peakAssignmentNamed $paName]
	    PeakAssignment -this $curPA 

	    $curPA incrementNumFiltersFailed 
	    incr numPeakAssignsFailed
	    if {[$curPA isGood]} {
		incr numGoodPeakAssignsFailed
	    }
		
	    rename $curPA ""
	}
	
	rename $curPeak ""
    }

    #
    # summarize results
    #

    if {$remVar != ""} {

	if {$maxUpperBound < 99999.9} {
	    set curRpt "SymmetryFilter, using shift tolerance of $tolerance and maxUpperBound of $maxUpperBound :"
	} else {
	    set curRpt "SymmetryFilter, using shift tolerance of $tolerance :"
	}
    
	appendLineToString curRpt "$numNoPossSymmPeakAssigns peak assignments have no possible symmetry partner"
	
	if {$numPeakAssignsKeptForHighUpperBound > 0} {
	    appendLineToString curRpt "$numPeakAssignsKeptForHighUpperBound asymmetric peak assignments were kept because they had  high upper bounds"
	}

	appendLineToString curRpt "$numPeakAssignsFailed peak assignments failed this filter for missing expected symmetry"
	    
	if {$numGoodPeakAssignsFailed > 0} {
	    appendLineToString curRpt "   Of them, $numGoodPeakAssignsFailed PeakAssignments were marked as good"
	}

	appendLineToString curRpt "$numPeaksFailed peaks had all of their peakAssignments fail this filter for missing expected symmetry"
	    
	if {$numGoodPeaksFailed > 0} {
	    appendLineToString curRpt "   Of them, $numGoodPeaksFailed Peaks were marked as good"
	}

	upvar $remVar tempRem
	lappend tempRem $curRpt
    }
    return ""
}


proc primarySequenceDistanceFilter args {

    set pot            [requiredFlagVal $args -pot]
    set intra          [flagExists $args -intraresidue]
    set bbnseq         [flagExists $args -backbone-sequential]
    set seq            [flagExists $args -sequential]
    set remVar         [flagVal $args -remarksVariableName ""]

    #
    # the PSD filter will increment the number of failed filters for  any peakAssignment
    # of a peak if that peak also has an intraresidue PA with equal or fewer failed filters 
    #
    
    # IMPROVEMENT for later:  Check for equivalent intra/seq PAs.  If a particular intra/seq PA 
    # can be assigned to >1 peak (and one hasn't been filtered by the stripe filter), don't apply 
    # the PSD filter to either peak

    set PAsToFail [list]

    #
    # Filter peaks with intraresidue peakAssignment(s)
    #

    if {$intra || $bbnseq || $seq} {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak

	    #
	    # find the lowest number of failed filters of this peak's intraresidue PAs
	    # Note that peaks w/ no intrares PAs will get lowestFailedFilters of 999
	    #

	    set lowestFailedFilters 999
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		if {[$curPA isIntraresidue]} {
		    set lowestFailedFilters [min [list $lowestFailedFilters [$curPA numFiltersFailed]]]
		}
		rename $curPA ""
	    }
		
	    #
	    # record peakAssignments with PSDs > intraresidue and
	    # numFailedFilters >= lowestFailedFilters
	    #
		
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA

		if {(! [$curPA isIntraresidue]) && 
		    ([$curPA numFiltersFailed] >= $lowestFailedFilters)} {
			
		    lappend PAsToFail [$curPeak name]
		    lappend PAsToFail [$curPA name]
		}	
		rename $curPA ""
	    }
	    rename $curPeak ""
	}
	
    }


    #
    # Filter peaks with bbn-seq peakAssignment(s)
    #

    if {$bbnseq} {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak

	    #
	    # find the lowest number of failed filters of this peak's intraresidue PAs
	    # Note that peaks w/ no bbn-seq PAs will get lowestFailedFilters of 999
	    #

	    set lowestFailedFilters 999
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		if {[$curPA isBackboneSequential]} {
		    set lowestFailedFilters [min [list $lowestFailedFilters [$curPA numFiltersFailed]]]
		}
		rename $curPA ""
	    }
		
	    #
	    # record peakAssignments with PSDs > intraresidue and
	    # numFailedFilters >= lowestFailedFilters
	    #
		
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA

		if {(! [$curPA isBackboneSequential]) && 
		    ([$curPA numFiltersFailed] >= $lowestFailedFilters)} {
			
		    lappend PAsToFail [$curPeak name]
		    lappend PAsToFail [$curPA name]
		}	
		rename $curPA ""
	    }
	    rename $curPeak ""
	}
	
    }


    #
    # Filter peaks with sequential peakAssignment(s)
    #

    if {$seq} {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak

	    #
	    # find the lowest numFailedFilters of this peak's sequential PAs
	    # Note that peaks w/ no sequential PAs will get lowestFailedFilters of 999
	    #
	    # Should this value be calculated using sequential PAs only, or w/ seq + intra PAs
	    # together?
	    #

	    set lowestFailedFilters 999
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		if {[$curPA isSequential]} {
		    set lowestFailedFilters [min [list $lowestFailedFilters [$curPA numFiltersFailed]]]
		}
		rename $curPA ""
	    }
		
	    #
	    # record peakAssignments with PSDs > sequential and 
	    # numFailedFilters >= lowestFailedFilters 
	    #
		
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA

		if {(! ([$curPA isSequential] || [$curPA isIntraresidue])) && 
		    ([$curPA numFiltersFailed] >= $lowestFailedFilters)} {

		    lappend PAsToFail [$curPeak name]
		    lappend PAsToFail [$curPA name]
		}
		rename $curPA ""
	    }
	    rename $curPeak ""
	}
    }

    set nPAsFailed 0
    set nGoodPAsFailed 0

    set peaksWithPAsFailed [list]

    foreach {peakName paName} $PAsToFail {

	set curPeak [$pot peakNamed $peakName]
	Peak -this $curPeak
	set curPA [$curPeak peakAssignmentNamed $paName]
	PeakAssignment -this $curPA

	if {[$curPA isGood]} {
	    incr nGoodPAsFailed 
	}
	
	incr nPAsFailed
	lappend peaksWithPAsFailed [$curPeak name]

	$curPA incrementNumFiltersFailed

	rename $curPA ""
	rename $curPeak ""
    }

    set nPeaksWithPAsFailed [llength [lsort -unique $peaksWithPAsFailed]]

    

    if {$remVar != ""} {
	upvar $remVar tempRem

	set rpt "PrimarySequenceDistance filter: "
	if {$intra} {
	    appendLineToString rpt "   Applying filter to peaks with intra-residue peakAssignments"
	}

	if {$bbnseq} {
	    appendLineToString rpt "   Applying filter to peaks with backbone-sequential peakAssignments"
	}

	if {$seq} {
	    appendLineToString rpt "   Applying filter to peaks with sequential peakAssignments"
	}

	set nPeaksFiltered 0
	set nPeaksAssigned 0

	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    
	    if {[$curPeak isAssigned]} {
		incr nPeaksAssigned
	    }

	    if {($intra && [$curPeak isIntraresidue]) ||
		($bbnseq  && [$curPeak isBackboneSequential]) ||
		($seq && [$curPeak isSequential])} {

		incr nPeaksFiltered
	    }
	    rename $curPeak ""
	}

	appendLineToString rpt [format "   This filter was applied to %d peaks (%f%% of total assigned peaks)" $nPeaksFiltered [expr 100 * $nPeaksFiltered / double($nPeaksAssigned)]]
	appendLineToString rpt [format "   %d peakAssignments on %d peaks (%f%% of filtered peaks) failed this filter." $nPAsFailed $nPeaksWithPAsFailed \
				    [expr 100 * $nPeaksWithPAsFailed / double($nPeaksFiltered)]]

	if {$nGoodPAsFailed > 0} {
	    appendLineToString rpt [format "   %d of those peakAssignments were marked as good" $nGoodPAsFailed]
	}

	lappend tempRem $rpt
    }
    return ""
}
