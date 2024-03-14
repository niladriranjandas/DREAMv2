#
# procedures to produce a randomized NOE list
#
# JJK 4/24/02
#

package provide noerandomizer 4.0
package require marvin

#
# Given a list of shiftAssignments, bounds, and minimal/maximal acceptable violations,
# return a pair of randomly chosen shiftAssigns that have an acceptable violation
#

proc generateRandomShiftAssignmentPair args {

    set saList        [requiredFlagVal $args -shiftAssignments]
    set upBound       [flagVal $args -upBound 5.0]
    set minAcceptViol [flagVal $args -minAcceptableViol 0.5]
    set maxAcceptViol [flagVal $args -maxAcceptableViol 99999999.9]

    set minUpBound [expr $upBound + $minAcceptViol]
    set maxUpBound [expr $upBound + $maxAcceptViol]

    set fromSAs [list]
    set toSAs   [list]

    foreach elem $saList {
	ShiftAssignment -this $elem
	if {[$elem isFrom]} {
	    lappend fromSAs [$elem cget -this]
	} else {
	    lappend toSAs [$elem cget -this]
	}
	rename $elem ""
    }


    set randomFrom [chooseRandomItem $fromSAs]
    set randomTo   [chooseRandomItem $toSAs]

    ShiftAssignment -this $randomFrom
    ShiftAssignment -this $randomTo

    set curDist [$randomFrom distanceToShiftAssignment [$randomTo cget -this]]

    #
    # make sure its violation is between min & max acceptable viols
    #

    while {! (($curDist > $minUpBound) &&
	      ($curDist < $maxUpBound))}  {

	rename $randomFrom ""
	rename $randomTo ""

	set randomFrom [chooseRandomItem $fromSAs]
	set randomTo   [chooseRandomItem $toSAs]

	ShiftAssignment -this $randomFrom
	ShiftAssignment -this $randomTo

	set curDist [$randomFrom distanceToShiftAssignment [$randomTo cget -this]]
    }

    return [list [$randomFrom cget -this] [$randomTo cget -this]]
}



#
# Generate a bunch of single-assignment random peaks
# and add them to a pot
#

proc addBadPeaks args {

    set outputPot     [requiredFlagVal $args -outputPot]
    set saList        [requiredFlagVal $args -shiftAssignments]
    set numToAdd      [flagVal $args -numToAdd 500]
    set namePrefix    [flagVal $args -namePrefix "random"]
    set upBound       [flagVal $args -upBound 5.0]
    set lowBound      [flagVal $args -lowBound 1.8]
    set minAcceptViol [flagVal $args -minAcceptableViol 0.5]
    set maxAcceptViol [flagVal $args -maxAcceptableViol 99999999.9]

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
    # for each new peak to generate,
    #
    
    for {set count 1} {$count <= $numToAdd} {incr count} {
	
	#
	# generate a random peakAssignment
	#
		
	set SAs [generateRandomShiftAssignmentPair \
		     -shiftAssignments $saList \
		     -upBound $upBound \
		     -minAcceptableViol $minAcceptViol \
		     -maxAcceptableViol $maxAcceptViol]

	set tempPA [PeakAssignment -args [format "%s_%d_0" $namePrefix $count]]	  
	$tempPA setNMono $nMono
	$tempPA setAveExp $aveExp
	
	$tempPA setFromAssignment [lindex $SAs 0]
	$tempPA setToAssignment   [lindex $SAs 1]
	$tempPA appendToNote "Randomly generated bad peakAssignment"

	#
	# generate a random peak
	# and attach the random peakAssignment
	#
	
	set tempPeak [Peak -args [format "%s_%d" $namePrefix $count]]
	$tempPeak setUpBound $upBound
	$tempPeak setLowBound $lowBound
	$tempPeak appendToNote "Randomly generated bad peak"
	
	$tempPeak addPeakAssignment [$tempPA cget -this]
	$outputPot addPeak [$tempPeak cget -this]

	$tempPA -disown
	rename $tempPA ""
	    
	$tempPeak -disown
	rename $tempPeak ""
    }
}



# 
# Given a pot, add a given number of random bad peakAssignments
# to randomly chosen elements of its peaks list
#

proc addBadPeakAssigns args {

    set outputPot     [requiredFlagVal $args -outputPot]
    set saList        [requiredFlagVal $args -shiftAssignments]
    set numToAdd      [flagVal $args -numToAdd 500]
    set minAcceptViol [flagVal $args -minAcceptableViol 0.5]
    set maxAcceptViol [flagVal $args -maxAcceptableViol 99999999.9]

    set peakList [$outputPot peaks]

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
    # for each new peakAssignment to add,
    #
    
    for {set count 1} {$count <= $numToAdd} {incr count} {

	#
	# grab a peak at random
	#

	set curPeak [chooseRandomItem $peakList]

	Peak -this $curPeak

	set curName [$curPeak name]
	set nAssigns [$curPeak numPeakAssignments]

	set upBound [$curPeak upBound]
	set lowBound [$curPeak lowBound]

	#
	# generate a random peakAssignment
	#
		
	set SAs [generateRandomShiftAssignmentPair \
		     -shiftAssignments $saList \
		     -upBound $upBound \
		     -minAcceptableViol $minAcceptViol \
		     -maxAcceptableViol $maxAcceptViol]
	
	set tempPA [PeakAssignment -args [format "%s_%d" $curName $nAssigns]]	  
	$tempPA setNMono $nMono
	$tempPA setAveExp $aveExp
	
	$tempPA setFromAssignment [lindex $SAs 0]
	$tempPA setToAssignment   [lindex $SAs 1]
	$tempPA appendToNote "Randomly generated bad peakAssignment"

	#
	# attach the new peakAssignment to the chosen peak
	#

	$curPeak addPeakAssignment [$tempPA cget -this]
	    
	$tempPA -disown
	rename $tempPA ""
	rename $curPeak ""
    }
}



proc removeRandomPeaks args {

    set pot    [requiredFlagVal $args -pot]
    set frac   [flagVal $args -frac]
    set num    [flagVal $args -num]
    set remVar [flagVal $args -remarksVariableName ""]

    if {! ([flagExists $args -frac] || [flagExists $args -num])} {
	error "Neither -frac nor -num flags defined"
    }

    if {[flagExists $args -frac]} {
	set num [expr round($frac * [$pot numPeaks])]
	if {$num > [$pot numPeaks]} {
	    set $num [$pot numPeaks]
	}
    }

    set peakList [shuffledList [$pot peaks]]
    
    for {set count 0} {$count < $num} {incr count} {
	set curPeak [lindex $peakList $count]
	Peak -this $curPeak
	$pot removePeakNamed [$curPeak name]
	$curPeak -acquire
	rename $curPeak ""
    }

    if {$remVar != ""} {
	upvar $remVar tempRem
	
	set rpt [format "RemoveRandomPeaks: %d peaks were removed at random" $num]
	lappend tempRem $rpt
    }
    return ""
}


proc swapShifts args {
    
    set shiftList [requiredFlagVal $args -shiftList]
    set residA    [requiredFlagVal $args -resNumA]
    set residB    [requiredFlagVal $args -resNumB]
    set remVar    [flagVal $args -remarksVariableName ""]

    set nSwapped 0

    set retVal [list]

    foreach elem $shiftList {

	set curShift [lindex $elem 0]
	set curSel   [lindex $elem 1]
	    
	if {[string match "*resid $residA *" $curSel]} {
	    
	    regsub "resid $residA" $curSel "resid $residB" temp
	    set curSel $temp
	    incr nSwapped
	    
	} elseif {[string match "*resid $residB *" $curSel]} {
	    
	    regsub "resid $residB" $curSel "resid $residA" temp
	    set curSel $temp
	    incr nSwapped
	}

	lappend retVal [list $curShift $curSel]
    }


    if {$remVar != ""} {
	upvar $remVar tempRem
	
	set rpt [format "SwapShiftAssigns: %d shift table entries were swapped between residues %s and %s" $nSwapped $residA $residB]
	lappend tempRem $rpt
    }

    return $retVal
}
