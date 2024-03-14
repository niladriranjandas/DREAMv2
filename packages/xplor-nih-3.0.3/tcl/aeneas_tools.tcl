#
# TCL procedures for AENEAS:  Automated Experimental NoE Assignments based on Shifts
#
# (Props to Gabriel & Claudia Cornilescu for the name)
#
# A simple system for generating initial restraints for Marvin, based on 
# an input chemical shift table and an input NOE peak list (with shift locations 
# and intensities).
# 
# JJK 2/6/04
#

package provide aeneas 0.1
package require marvin
package require newstripe
package require primseqexcept

#
# Given an unfolded peak position from unfoldPeakPosnToMatch, a raw shiftAssignment, and a tolerance (all in PPM),
# return whether the unfolded position is within tolerance of the shiftAssignment.
#
# Note that this assumes that tolerance < spectralWidth
#

proc unfoldedPeakPosnMatches {unfoldedPeakPosn targetPosn tolerance} {

    return [expr abs([lindex $unfoldedPeakPosn 0] - $targetPosn) < $tolerance]
}


#
# Given a raw peak position, a raw shift assignment, and the top/bot edges of the spectrum (all in PPM), 
# a folded/aliased flag, and a sign changes flag,
#
# Return the un-folded/aliased peak positions closest to the shiftAssignment and their expected signs
# as [unfoldedPeakPos, expectedSign] pairs
#

proc unfoldPeakPosnToMatch {peakPosn targetPosn isFolded top bot signChangesUponFolding} {

    set specWid  [expr double($top - $bot)]
    set specWid2 [expr 2 * $specWid]

    if {! $isFolded} {

	#
	# find closest positive aliased position
	#

	set expectedSign      1
	set startPosn         $peakPosn
	set nAliasingsToApply [expr round(($targetPosn - $startPosn) / $specWid2)]
	set unaliasedPeakPosn [expr $startPosn + ($nAliasingsToApply * $specWid2)]
	
	set retVal1 [list $unaliasedPeakPosn $expectedSign]

	#
	# find closest negative aliased position
	#

	if {$signChangesUponFolding} {
	    set expectedSign -1
	} else {
	    set expectedSign 1
	}

	set startPosn         [expr $peakPosn - $specWid]
	set nAliasingsToApply [expr round(($targetPosn - $startPosn) / $specWid2)]
	set unaliasedPeakPosn [expr $startPosn + ($nAliasingsToApply * $specWid2)]

	set retVal2 [list $unaliasedPeakPosn $expectedSign]

    } else {

	set retVal [list]

	#
	# find closest positive folded position
	#

	
	set expectedSign     1
	set startPosn        $peakPosn
	set nFoldsToApply    [expr round(($targetPosn - $startPosn) / $specWid2)]
	set unfoldedPeakPosn [expr $startPosn + ($nFoldsToApply * $specWid2)]

	set retVal1 [list $unfoldedPeakPosn $expectedSign]

	#
	# try matching to negative folded position
	#

	if {$signChangesUponFolding} {
	    set expectedSign -1
	} else {
	    set expectedSign 1
	}

	set startPosn        [expr $peakPosn - 2 * double($peakPosn - $bot)]
	set nFoldsToApply    [expr round(($targetPosn - $startPosn) / $specWid2)]
	set unfoldedPeakPosn [expr $startPosn + ($nFoldsToApply * $specWid2)]

	set retVal2 [list $unfoldedPeakPosn $expectedSign]
    }

    return [list $retVal1 $retVal2]
}
     

#
# this is now only used to determine folded positions of shiftAssignments for 
# the diagonal check in creating explicit inverse exceptions
#

proc foldedOrAliasedPosition {rawShift top bot signChangeFlag isFolded} {

    set expectedSign 1

    set foldedShift $rawShift

    # this loops to allow aliased positions that are far outside the spectral
    # width to be calculated correctly

    set nFolds 0

    if {$isFolded} {

	while {($foldedShift > $top) || ($foldedShift < $bot)} {
	    incr nFolds
	    set expectedSign [expr $expectedSign * -1]
	    if {$foldedShift > $top} {
		set foldedShift [expr $top - ($foldedShift - $top)]
	    } elseif {$foldedShift < $bot} {
		set foldedShift [expr $bot + ($bot - $foldedShift)]
	    } 
	}
    
    } else {
	
	while {($foldedShift > $top) || ($foldedShift < $bot)} {
	    incr nFolds
	    set expectedSign [expr $expectedSign * -1]
	    if {$foldedShift > $top} {
		set foldedShift [expr $foldedShift - ($top - $bot)]
	    } elseif {$foldedShift < $bot} {
		set foldedShift [expr $foldedShift + ($top - $bot)]
	    } 
	}
    }

	
    if {! $signChangeFlag} {
	set $expectedSign 1
    }

    return [list $foldedShift $expectedSign]
}


proc recordUnfoldedPositions args {

    set pot                   [flagVal $args -pot]
    set saList                [flagVal $args -shiftAssignments ]

    if {($pot == "") && ($saList == "")} {
	error "Neither -pot nor -shiftAssignments defined in call to recordUnfoldedPositions"
    }

    if {$saList == ""} {
	set saList [$pot shiftAssignments]
    }

    set fromProtonRange       [flagVal $args -fromProtonSpectralRange [list 0 1]]
    set fromHeavyRange        [flagVal $args -fromHeavyatomSpectralRange [list 0 1]]
    set toProtonRange         [flagVal $args -toProtonSpectralRange [list 0 1]]
    set toHeavyRange          [flagVal $args -toHeavyatomSpectralRange [list 0 1]]

    set fromProtonSignChanges [flagVal $args -fromProtonSignChanges 0]
    set fromHeavySignChanges  [flagVal $args -fromHeavyatomSignChanges 0]
    set toProtonSignChanges   [flagVal $args -toProtonSignChanges 0]
    set toHeavySignChanges    [flagVal $args -toHeavyatomSignChanges 0]

    set foldedAlongFromProton [flagVal $args -shiftsFoldedAlongFromProtonDimension 0]
    set foldedAlongFromHeavy  [flagVal $args -shiftsFoldedAlongFromHeavyatomDimension 0]
    set foldedAlongToProton   [flagVal $args -shiftsFoldedAlongToProtonDimension 0]
    set foldedAlongToHeavy    [flagVal $args -shiftsFoldedAlongToHeavyatomDimension 0]

    set fromProtonTop [max $fromProtonRange]
    set fromProtonBot [min $fromProtonRange]
    set fromHeavyTop  [max $fromHeavyRange]
    set fromHeavyBot  [min $fromHeavyRange]
    set toProtonTop   [max $toProtonRange]
    set toProtonBot   [min $toProtonRange]
    set toHeavyTop    [max $toHeavyRange]
    set toHeavyBot    [min $toHeavyRange]

    set count 0

    foreach curSA $saList {

	ShiftAssignment -this $curSA 

	updateUser [format "Recording unfolded positions of %s (%d of %d) \r" [$curSA name] [incr count] [llength $saList]]

	set curSign 1

	if {[$curSA isFrom]} {

	    if {[$curSA hasProtonShift]} {
		set foldedPos  [foldedOrAliasedPosition [$curSA protonShift] $fromProtonTop $fromProtonBot $fromProtonSignChanges $foldedAlongFromProton]
		$curSA setFoldedProtonShift [lindex $foldedPos 0]
		set curSign [expr $curSign * [lindex $foldedPos 1]]
	    }

	    if {[$curSA hasHeavyatomShift]} {
		set foldedPos [foldedOrAliasedPosition [$curSA heavyatomShift] $fromHeavyTop $fromHeavyBot $fromHeavySignChanges $foldedAlongFromHeavy]
		$curSA setFoldedHeavyatomShift [lindex $foldedPos 0]
		set curSign [expr $curSign * [lindex $foldedPos 1]]
	    }

	} else {

	    if {[$curSA hasProtonShift]} {
		set foldedPos  [foldedOrAliasedPosition [$curSA protonShift] $toProtonTop $toProtonBot $toProtonSignChanges $foldedAlongToProton]
		$curSA setFoldedProtonShift [lindex $foldedPos 0]
		set curSign [expr $curSign * [lindex $foldedPos 1]]
	    }

	    if {[$curSA hasHeavyatomShift]} {
		set foldedPos [foldedOrAliasedPosition [$curSA heavyatomShift] $toHeavyTop $toHeavyBot $toHeavySignChanges $foldedAlongToHeavy]
		$curSA setFoldedHeavyatomShift [lindex $foldedPos 0]
		set curSign [expr $curSign * [lindex $foldedPos 1]]
	    }
	}

	$curSA setFoldedSign $curSign
	rename $curSA ""
    }
}

    
#
# sanity check that each peak's position is within the range
# specified by the top and bottom
#

proc peakOutsideSpectrum {peakPosn top bot} {

    if {($peakPosn > $top) || ($peakPosn < $bot)} {
	return 1
    } else {
	return 0
    }
}


#
# a way of caching chemical type selectivity for match*d, below.
#
 
proc shiftEntriesThatMatchChemicalType {shiftList selTargetString} {

    set retVal [list]

    set selTarget [AtomSel -args $selTargetString]

    foreach shiftEntry $shiftList {

	set curShift    [lindex $shiftEntry 0]        
	set cleanEntry [list $curShift]
	
	foreach sel [lrange $shiftEntry 1 end] {
	    if {[$sel intersects [$selTarget cget -this]]} {
		lappend cleanEntry $sel
	    }
	}
	
	if {[llength $cleanEntry] > 1} {
	    lappend retVal $cleanEntry
	}
    }
    
    rename $selTarget ""
    
    return $retVal
}


#
# A tighter tolerance check is possible once all proposed assignment
# shifts are known.  Uses an ellipsoid that just fits inside the 
# rectange defined by the tolerances.  Its volume is slightly more
# than half the volume of the rectangle.
#
# Modified to handle cases where 1 or more shifts are given the value "NONE".
# Just ignores that dimension, which doesn't break the math at all
#

proc ellipseCheck {peakShifts assignShifts tolerances} {

    set tot 0

    foreach peakShift $peakShifts assignShift $assignShifts tol $tolerances {

	if {($peakShift != "NONE") && ($assignShift != "NONE")} {

	    set tot [expr $tot - pow(($peakShift - $assignShift) / $tol, 2)]
	}
    }

    if {$tot > -1.0} {
	return 1
    } else {
       	return 0
    }
}

#
# A better check for use in simpleMatch--
# there's no reason errors in one atom's shift table entry
# should correlate with another atom's
#
# Modified to handle cases where 1 or more shifts are given the value "NONE".
# Just ignores that dimension, which doesn't break the logic at all
#

proc rectangleCheck {peakShifts assignShifts tolerances} {

    set ok 1

    foreach peakShift $peakShifts assignShift $assignShifts tol $tolerances {

	if {($peakShift != "NONE") && ($assignShift != "NONE")} {

	    if {[expr abs($peakShift - $assignShift)] > $tol} {
		set ok 0
		break
	    }
	}
    }

    return $ok
}


# 
# My version of STAPP
# 


proc match2d args {

    set pot                   [requiredFlagVal $args -pot]
    set fromProtonTol         [requiredFlagVal $args -fromProtonTolerancePPM]
    set fromProtonRange       [requiredFlagVal $args -fromProtonSpectralRangePPM]
    set toProtonTol           [requiredFlagVal $args -toProtonTolerancePPM]
    set toProtonRange         [requiredFlagVal $args -toProtonSpectralRangePPM]
    set peakList              [flagVal $args -peakList [$pot peaks]]
    set SAList                [flagVal $args -shiftAssignList [$pot shiftAssignments]]
    set remVar                [flagVal $args -remarksVariableName ""]
    set minSALikelihood       [flagVal $args -minShiftAssignmentLikelihood -1]
    set verboseOutsideSpectrum [flagExists $args -verboseOutsideSpectrum]

    set fromProtonSignChanges [flagExists $args [list -signChangesUponFoldingFromProtonDimension \
						     -signChangesUponAliasingFromProtonDimension]]

    set toProtonSignChanges   [flagExists $args [list -signChangesUponFoldingToProtonDimension \
						     -signChangesUponAliasingToProtonDimension]]

    set useIndividualTols     [flagExists $args -useIndividualTolerances]

    set foldedAlongFromProton [flagExists $args -shiftsFoldedAlongFromProtonDimension]
    set foldedAlongToProton   [flagExists $args -shiftsFoldedAlongToProtonDimension]

    marvinPyth command {import pasd}
    set defAveExp [lindex \
		    [lindex \
			 [marvinPyth command {out=pasd.aveExp} \
			      {} {out}] \
			 0] 1]
    set aveExp                [flagVal $args -aveExp $defAveExp]

    set unfoldedSign 0

    if {[flagExists $args -nonFoldedPeaksArePositive]} {
	set unfoldedSign 1
    }

    if {[flagExists $args -nonFoldedPeaksAreNegative]} {
	set unfoldedSign -1
    }

    if {$unfoldedSign == 0} {
	error "Missing either -nonFoldedPeaksArePositive or -nonFoldedPeaksAreNegative flag"
    }

    set ignoreSign [flagExists $args -ignoreSign]

    #
    # extract spectral top/bottom from range flags
    #

    set fromProtonTop [max $fromProtonRange]
    set fromProtonBot [min $fromProtonRange]
    set toProtonTop   [max $toProtonRange]
    set toProtonBot   [min $toProtonRange]

    #
    # assemble lists of from and to shiftAssigns that can be used--ie., 
    # those that have all needed chemical shifts, and appropriate previousLikelihoods
    #

    set fromSAs [list]
    set toSAs   [list]

    foreach curSA $SAList {
	ShiftAssignment -this $curSA
	if {[$curSA isFrom] && [$curSA hasProtonShift] && ([$curSA previousLikelihood] >= $minSALikelihood)} {
	    lappend fromSAs [$curSA cget -this]
	} elseif {[$curSA isTo] && [$curSA hasProtonShift] && ([$curSA previousLikelihood] >= $minSALikelihood)} {
	    lappend toSAs [$curSA cget -this]
	}
	rename $curSA ""
    }

    #
    # Need to record unfolded positions of each shiftAssign for use in 
    # explicit inverse exceptions.  Do it here because I have the spectral 
    # ranges and so forth
    #
    
    recordUnfoldedPositions \
	-shiftAssignments                        [concat $fromSAs $toSAs] \
	-fromProtonSpectralRange                 $fromProtonRange \
	-toProtonSpectralRange                   $toProtonRange \
	-fromProtonSignChanges                   $fromProtonSignChanges \
	-toProtonSignChanges                     $toProtonSignChanges \
	-shiftsFoldedAlongFromProtonDimension    $foldedAlongFromProton \
	-shiftsFoldedAlongToProtonDimension      $foldedAlongToProton 

    #
    # remove any linked peakAssignments that exist
    #

    $pot resetLinkedPeakAssignmentNames

    #
    # For each peak, try unfolding its position to best match all 
    # from and to shiftAssignments.  
    #
    # Generate all possible pairings of from and to ShiftAssignments 
    # whose shifts are within tolerance of the peak's unfolded position. 
    #
    # Make sure the peak's sign matches the expected value (which depends 
    # on how many times it was unfolded to match the shiftAssignments).  
    #
    # Record those pairings of ShiftAssignments as assignments
    #


    #
    # could get a big speedup if I 
    # 1.  sort from and to shift lists by folded proton shift
    # 2.  write a proc that, given the sorted shift lists and 
    #     max/min values of interest (from peak posn +/- tol),
    #     returns list of shift entries in range.
    # 3.  run foreaches over these relevant sublists
    #

    set numPeaks [llength $peakList]
    set curPeakNum 0
    set numUnassignablePeaks 0
    set numAssignedPeaks 0
    set totNumPeakAssignments 0
    set peaksOutsideSpectrum [list]

    marvinPyth command {import pasd}
    set nMono [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.nMono} \
			   {} {out}] \
		      0] 1]
    
    foreach curPeak $peakList {

	Peak -this $curPeak

	updateUser [format "Checking for matches against peak %s (%d of %d) \r" \
			[$curPeak name] [incr curPeakNum] $numPeaks]

	#
	# remove any existing shift and peakAssignments from the current peak
	#

	$curPeak removeAllPeakAssignments
	$curPeak removeAllShiftAssignments

	#
	# check that this peak's position is inside the spectrum
	#

	if {[peakOutsideSpectrum [$curPeak fromProtonShift] $fromProtonTop $fromProtonBot]} {

	    set err [format "from proton location at %f is outside of spectrum (%f .. %f)" \
			 [$curPeak fromProtonShift] $fromProtonBot $fromProtonTop ]


	    $curPeak appendToNote $err
	    
	    updateUser [format "Peak %s %s \n" [$curPeak name] $err]

	    lappend peaksOutsideSpectrum [format "Peak %s %s" [$curPeak name] $err]
	}

	if {[peakOutsideSpectrum [$curPeak toProtonShift] $toProtonTop $toProtonBot]} {

	    set err [format "to proton location at %f is outside of spectrum (%f .. %f)" \
			 [$curPeak toProtonShift] $toProtonBot $toProtonTop ]

	    $curPeak appendToNote $err
	    
	    updateUser [format "Peak %s %s \n" [$curPeak name] $err]

	    lappend peaksOutsideSpectrum [format "Peak %s %s" [$curPeak name] $err]
	}

	#
	# find from shiftAssigns that match the restraint's from dimensions.
	# Record each one's name, expected sign, and unfolded peak positions
	#

	set matchingFromShiftAssigns [list]

	foreach fromShiftAssign $fromSAs {	    
	    ShiftAssignment -this $fromShiftAssign

	    if {$useIndividualTols && [$fromShiftAssign hasProtonTolerance]} {
		set fpt [$fromShiftAssign protonTolerance]
	    } else {
		set fpt $fromProtonTol
	    }

	    if {! [$fromShiftAssign hasProtonShift]} {
		rename $fromShiftAssign ""
		continue
	    }

	    set unfoldedFromProtonPPs [unfoldPeakPosnToMatch [$curPeak fromProtonShift]    [$fromShiftAssign protonShift]    $foldedAlongFromProton $fromProtonTop $fromProtonBot $fromProtonSignChanges] 

	    foreach protonPP $unfoldedFromProtonPPs {
		if {[unfoldedPeakPosnMatches $protonPP [$fromShiftAssign protonShift] $fpt]} {
			    
		    if {! [$curPeak hasShiftAssignmentNamed [$fromShiftAssign name]]} {
			$curPeak addShiftAssignment [$fromShiftAssign cget -this]
		    }

		    set curMatchSign [lindex $protonPP 1]
		    lappend matchingFromShiftAssigns [list [$fromShiftAssign name] $curMatchSign [lindex $protonPP 0]]
		}
	    }
	    
	    rename $fromShiftAssign ""
	}
		
	#
	# find to shiftAssigns that match the restraint's to dimension
	# Record each one's name, expected sign, and unfolded peak position
	#

	set matchingToShiftAssigns [list]

	foreach toShiftAssign $toSAs {	    
	    ShiftAssignment -this $toShiftAssign
	    
	    if {$useIndividualTols && [$toShiftAssign hasProtonTolerance]} {
		set tpt [$toShiftAssign protonTolerance]
	    } else {
		set tpt $toProtonTol
	    }

	    if {! [$toShiftAssign hasProtonShift]} {
		rename $toShiftAssign ""
		continue
	    }

	    set unfoldedToProtonPPs [unfoldPeakPosnToMatch [$curPeak toProtonShift] [$toShiftAssign protonShift] $foldedAlongToProton $toProtonTop $toProtonBot $toProtonSignChanges]

	    foreach protonPP $unfoldedToProtonPPs {
		if {[unfoldedPeakPosnMatches $protonPP [$toShiftAssign protonShift] $tpt]} {

		    if {! [$curPeak hasShiftAssignmentNamed [$toShiftAssign name]]} {
			$curPeak addShiftAssignment [$toShiftAssign cget -this]
		    }

		    set curMatchSign [lindex $protonPP 1]
		    lappend matchingToShiftAssigns [list [$toShiftAssign name] $curMatchSign [lindex $protonPP 0]]
		}
	    }

	    rename $toShiftAssign ""
	}

	#
	# now generate all possible pairings of from and to shiftAssigns that match the 
	# current peak, check that their expected signs match the peak's, 
	# and generate corresponding PeakAssignments
	#

	foreach fromMatch $matchingFromShiftAssigns {
	    set fromShiftAssign      [$pot shiftAssignmentNamed [lindex $fromMatch 0]]
	    ShiftAssignment -this $fromShiftAssign

	    set fromSign             [lindex $fromMatch 1]
	    set unfoldedFromProtonPP [lindex $fromMatch 2]

	    foreach toMatch $matchingToShiftAssigns {
		set toShiftAssign      [$pot shiftAssignmentNamed [lindex $toMatch 0]]
		ShiftAssignment -this $toShiftAssign

		set toSign             [lindex $toMatch 1]
		set unfoldedToProtonPP [lindex $toMatch 2]

		set expectedSign [expr $unfoldedSign * $fromSign * $toSign]

		if {(($expectedSign == [sign [$curPeak intensity]]) || $ignoreSign)} {
			
		    #
		    # record the PeakAssignment in the Peak
		    #
		    
		    incr totNumPeakAssignments
	    
		    set curPAName [format "%s_%d" [$curPeak name] [$curPeak numPeakAssignments]]
		    set curPA [PeakAssignment -args $curPAName]
		    $curPA setNMono $nMono
		    $curPA setAveExp $aveExp
	    
		    $curPeak addPeakAssignment [$curPA cget -this]
	
		    $curPA setFromAssignment [$fromShiftAssign cget -this]
		    $curPA setToAssignment   [$toShiftAssign   cget -this]
		    
		    $curPA setUnfoldedFromProtonPeakPosition    $unfoldedFromProtonPP
		    $curPA setUnfoldedToProtonPeakPosition      $unfoldedToProtonPP

		    $curPA resetNumFiltersFailed 
		    $curPA setPreviousLikelihood 1.0
		    
		    $curPA -disown
		    rename $curPA ""
		}
		rename $toShiftAssign ""
	    }
	    rename $fromShiftAssign ""
	}

	#
	# if I found no peakAssignments, complain
	#
	
	if {! [$curPeak isAssigned]} {
	    $curPeak appendToNote "No 2d assignments found."	    
	    incr numUnassignablePeaks
	} else {
	    incr numAssignedPeaks
	}

	rename $curPeak ""
    }

    if {$remVar != ""} {
	upvar $remVar tempRem

	set newRem [format "Matching %d 2-D peaks with %d from and %d to shiftAssignments" \
			[llength $peakList] [llength $fromSAs] [llength $toSAs]]

	if {$useIndividualTols} {
	    appendLineToString newRem "Using each shiftAssign's individual tolerances.  Defaults are: "
	}

	appendLineToString newRem  "from proton tolerance:     $fromProtonTol ppm"
	appendLineToString newRem  "to proton tolerance:       $toProtonTol ppm"

	appendLineToString newRem  [format "%d peaks were assigned with a total of %d assignments (mean degeneracy: %f)" \
					$numAssignedPeaks $totNumPeakAssignments [expr $totNumPeakAssignments / max(1.,double($numAssignedPeaks))]]
	appendLineToString newRem  [format " %d peaks (%f%% of the total) were found to be without assignments." \
					$numUnassignablePeaks [expr 100.0 * $numUnassignablePeaks / max(1.,double($numUnassignablePeaks + $numAssignedPeaks))]]
	
	lappend tempRem $newRem

	set newRem [format "%d peaks were located outside of the spectrum" [llength $peaksOutsideSpectrum]]
	if {$verboseOutsideSpectrum} {
	    foreach outsideError $peaksOutsideSpectrum {
		appendLineToString newRem $outsideError
	    }
	}

	lappend tempRem $newRem
    }

    return ""
}





proc match3d args {

    set pot                   [requiredFlagVal $args -pot]
    set fromProtonTol         [requiredFlagVal $args -fromProtonTolerancePPM]
    set fromProtonRange       [requiredFlagVal $args -fromProtonSpectralRangePPM]
    set fromHeavyTol          [requiredFlagVal $args -fromHeavyatomTolerancePPM]
    set fromHeavyRange        [requiredFlagVal $args -fromHeavyatomSpectralRangePPM]
    set toProtonTol           [requiredFlagVal $args -toProtonTolerancePPM]
    set toProtonRange         [requiredFlagVal $args -toProtonSpectralRangePPM]
    set peakList              [flagVal $args -peakList [$pot peaks]]
    set SAList                [flagVal $args -shiftAssignList [$pot shiftAssignments]]
    set remVar                [flagVal $args -remarksVariableName ""]
    set minSALikelihood       [flagVal $args -minShiftAssignmentLikelihood -1]
    set verboseOutsideSpectrum [flagExists $args -verboseOutsideSpectrum]

    set fromProtonSignChanges [flagExists $args [list -signChangesUponFoldingFromProtonDimension \
						     -signChangesUponAliasingFromProtonDimension]]

    set fromHeavySignChanges  [flagExists $args [list -signChangesUponFoldingFromHeavyatomDimension \
						     -signChangesUponAliasingFromHeavyatomDimension]]

    set toProtonSignChanges   [flagExists $args [list -signChangesUponFoldingToProtonDimension \
						     -signChangesUponAliasingToProtonDimension]]

    set useIndividualTols     [flagExists $args -useIndividualTolerances]

    set foldedAlongFromProton [flagExists $args -shiftsFoldedAlongFromProtonDimension]
    set foldedAlongFromHeavy  [flagExists $args -shiftsFoldedAlongFromHeavyatomDimension]
    set foldedAlongToProton   [flagExists $args -shiftsFoldedAlongToProtonDimension]

    set unfoldedSign 0

    if {[flagExists $args -nonFoldedPeaksArePositive]} {
	set unfoldedSign 1
    }

    if {[flagExists $args -nonFoldedPeaksAreNegative]} {
	set unfoldedSign -1
    }

    if {$unfoldedSign == 0} {
	error "Missing either -nonFoldedPeaksArePositive or -nonFoldedPeaksAreNegative flag"
    }

    set ignoreSign [flagExists $args -ignoreSign]

    marvinPyth command {import pasd}
    set defAveExp [lindex \
		    [lindex \
			 [marvinPyth command {out=pasd.aveExp} \
			      {} {out}] \
			 0] 1]
    set aveExp                [flagVal $args -aveExp $defAveExp]

    #
    # extract spectral top/bottom from range flags
    #

    set fromProtonTop [max $fromProtonRange]
    set fromProtonBot [min $fromProtonRange]
    set fromHeavyTop  [max $fromHeavyRange]
    set fromHeavyBot  [min $fromHeavyRange]
    set toProtonTop   [max $toProtonRange]
    set toProtonBot   [min $toProtonRange]

    #
    # assemble lists of from and to shiftAssigns that can be used--ie., 
    # those that have all needed chemical shifts, and appropriate previousLikelihoods
    #

    set fromSAs [list]
    set toSAs   [list]

    foreach curSA $SAList {
	ShiftAssignment -this $curSA
	if {[$curSA isFrom] && [$curSA hasProtonShift] && [$curSA hasHeavyatomShift] && ([$curSA previousLikelihood] >= $minSALikelihood)} {
	    lappend fromSAs [$curSA cget -this]
	} elseif {[$curSA isTo] && [$curSA hasProtonShift] && ([$curSA previousLikelihood] >= $minSALikelihood)} {
	    lappend toSAs [$curSA cget -this]
	}
	rename $curSA ""
    }

    #
    # Need to record unfolded positions of each shiftAssign for use in 
    # explicit inverse exceptions.  Do it here because I have the spectral 
    # ranges and so forth
    #
    
    recordUnfoldedPositions \
	-shiftAssignments                        [concat $fromSAs $toSAs] \
	-fromProtonSpectralRange                 $fromProtonRange \
	-fromHeavyatomSpectralRange              $fromHeavyRange \
	-toProtonSpectralRange                   $toProtonRange \
	-fromProtonSignChanges                   $fromProtonSignChanges \
	-fromHeavyatomSignChanges                $fromHeavySignChanges \
	-toProtonSignChanges                     $toProtonSignChanges \
	-shiftsFoldedAlongFromProtonDimension    $foldedAlongFromProton \
	-shiftsFoldedAlongFromHeavyatomDimension $foldedAlongFromHeavy \
	-shiftsFoldedAlongToProtonDimension      $foldedAlongToProton 

    #
    # remove any linked peakAssignments that exist
    #

    $pot resetLinkedPeakAssignmentNames

    #
    # For each peak, try unfolding its position to best match all 
    # from and to shiftAssignments.  
    #
    # Generate all possible pairings of from and to ShiftAssignments 
    # whose shifts are within tolerance of the peak's unfolded position. 
    #
    # Make sure the peak's sign matches the expected value (which depends 
    # on how many times it was unfolded to match the shiftAssignments).  
    #
    # Record those pairings of ShiftAssignments as assignments
    #


    #
    # could get a big speedup if I 
    # 1.  sort from and to shift lists by folded proton shift
    # 2.  write a proc that, given the sorted shift lists and 
    #     max/min values of interest (from peak posn +/- tol),
    #     returns list of shift entries in range.
    # 3.  run foreaches over these relevant sublists
    #

    set numPeaks [llength $peakList]
    set curPeakNum 0
    set numUnassignablePeaks 0
    set numAssignedPeaks 0
    set totNumPeakAssignments 0
    set peaksOutsideSpectrum [list]

    marvinPyth command {import pasd}
    set nMono [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.nMono} \
			   {} {out}] \
		      0] 1]

    foreach curPeak $peakList {

	Peak -this $curPeak

	updateUser [format "Checking for matches against peak %s (%d of %d) \r" \
			[$curPeak name] [incr curPeakNum] $numPeaks]

	#
	# remove any existing shift and peakAssignments from the current peak
	#

	$curPeak removeAllPeakAssignments
	$curPeak removeAllShiftAssignments

	#
	# check that this peak's position is inside the spectrum
	#


	if {[peakOutsideSpectrum [$curPeak fromProtonShift] $fromProtonTop $fromProtonBot]} {

	    set err [format "from proton location at %f is outside of spectrum (%f .. %f)" \
			 [$curPeak fromProtonShift] $fromProtonTop $fromProtonBot]


	    $curPeak appendToNote $err
	    
	    updateUser [format "Peak %s %s \n" [$curPeak name] $err]

	    lappend peaksOutsideSpectrum [format "Peak %s %s" [$curPeak name] $err]
	}

	if {[peakOutsideSpectrum [$curPeak fromHeavyatomShift] $fromHeavyTop $fromHeavyBot]} {

	    set err [format "from heavyatom location at %f is outside of spectrum (%f .. %f)" \
			 [$curPeak fromHeavyatomShift] $fromHeavyTop $fromHeavyBot]

	    $curPeak appendToNote $err
	    
	    updateUser [format "Peak %s %s \n" [$curPeak name] $err]
	    
	    lappend peaksOutsideSpectrum [format "Peak %s %s" [$curPeak name] $err]
	}


	if {[peakOutsideSpectrum [$curPeak toProtonShift] $toProtonTop $toProtonBot]} {

	    set err [format "to proton location at %f is outside of spectrum (%f .. %f)" \
			 [$curPeak toProtonShift] $toProtonTop $toProtonBot]

	    $curPeak appendToNote $err
	    
	    updateUser [format "Peak %s %s \n" [$curPeak name] $err]

	    lappend peaksOutsideSpectrum [format "Peak %s %s" [$curPeak name] $err]
	}

	#
	# find from shiftAssigns that match the restraint's from dimensions.
	# Record each one's name, expected sign, and unfolded peak positions
	#

	set matchingFromShiftAssigns [list]

	foreach fromShiftAssign $fromSAs {	    
	    ShiftAssignment -this $fromShiftAssign

	    if {$useIndividualTols && [$fromShiftAssign hasProtonTolerance]} {
		set fpt [$fromShiftAssign protonTolerance]
	    } else {
		set fpt $fromProtonTol
	    }

	    if {$useIndividualTols && [$fromShiftAssign hasHeavyatomTolerance]} {
		set fht [$fromShiftAssign heavyatomTolerance]
	    } else {
		set fht $fromHeavyTol
	    }

	    if {(! [$fromShiftAssign hasProtonShift]) || (! [$fromShiftAssign hasHeavyatomShift])} {
		rename $fromShiftAssign ""
		continue
	    }

	    set unfoldedFromProtonPPs [unfoldPeakPosnToMatch [$curPeak fromProtonShift]    [$fromShiftAssign protonShift]    $foldedAlongFromProton $fromProtonTop $fromProtonBot $fromProtonSignChanges] 
	    set unfoldedFromHeavyPPs  [unfoldPeakPosnToMatch [$curPeak fromHeavyatomShift] [$fromShiftAssign heavyatomShift] $foldedAlongFromHeavy  $fromHeavyTop  $fromHeavyBot  $fromHeavySignChanges] 

	    foreach protonPP $unfoldedFromProtonPPs {
		if {[unfoldedPeakPosnMatches $protonPP [$fromShiftAssign protonShift] $fpt]} {
		    foreach heavyPP $unfoldedFromHeavyPPs {
			if {[unfoldedPeakPosnMatches $heavyPP [$fromShiftAssign heavyatomShift] $fht]} {
			    
			    if {! [$curPeak hasShiftAssignmentNamed [$fromShiftAssign name]]} {
				$curPeak addShiftAssignment [$fromShiftAssign cget -this]
			    }

			    set curMatchSign [expr [lindex $protonPP 1] * [lindex $heavyPP 1]]
			    lappend matchingFromShiftAssigns [list [$fromShiftAssign name] $curMatchSign [lindex $protonPP 0] [lindex $heavyPP 0]]
			}
		    }
		}
	    }
	    
	    rename $fromShiftAssign ""
	}

	if { [llength $matchingFromShiftAssigns] == 0 } {
	    updateUser [format "%s No matching fromShiftAssignments" \
			    [$curPeak name]]
			
	}
	    

	#
	# find to shiftAssigns that match the restraint's to dimension
	# Record each one's name, expected sign, and unfolded peak position
	#

	set matchingToShiftAssigns [list]

	foreach toShiftAssign $toSAs {    
	    ShiftAssignment -this $toShiftAssign
	    
	    if {$useIndividualTols && [$toShiftAssign hasProtonTolerance]} {
		set tpt [$toShiftAssign protonTolerance]
	    } else {
		set tpt $toProtonTol
	    }

	    if {! [$toShiftAssign hasProtonShift]} {
		rename $toShiftAssign ""
		continue
	    }

	    set unfoldedToProtonPPs [unfoldPeakPosnToMatch [$curPeak toProtonShift] [$toShiftAssign protonShift] $foldedAlongToProton $toProtonTop $toProtonBot $toProtonSignChanges]

	    foreach protonPP $unfoldedToProtonPPs {
		if {[unfoldedPeakPosnMatches $protonPP [$toShiftAssign protonShift] $tpt]} {

		    if {! [$curPeak hasShiftAssignmentNamed [$toShiftAssign name]]} {
			$curPeak addShiftAssignment [$toShiftAssign cget -this]
		    }

		    set curMatchSign [lindex $protonPP 1]
		    lappend matchingToShiftAssigns [list [$toShiftAssign name] $curMatchSign [lindex $protonPP 0]]
		}
	    }

	    rename $toShiftAssign ""
	}

	#
	# now generate all possible pairings of from and to shiftAssigns that match the 
	# current restraint, check that their expected signs match the peak's, 
	# and generate corresponding PeakAssignments
	#

	foreach fromMatch $matchingFromShiftAssigns {
	    set fromShiftAssign      [$pot shiftAssignmentNamed [lindex $fromMatch 0]]
	    ShiftAssignment -this $fromShiftAssign

	    set fromSign             [lindex $fromMatch 1]
	    set unfoldedFromProtonPP [lindex $fromMatch 2]
	    set unfoldedFromHeavyPP  [lindex $fromMatch 3]

	    foreach toMatch $matchingToShiftAssigns {
		set toShiftAssign      [$pot shiftAssignmentNamed [lindex $toMatch 0]]
		ShiftAssignment -this $toShiftAssign

		set toSign             [lindex $toMatch 1]
		set unfoldedToProtonPP [lindex $toMatch 2]

		set expectedSign [expr $unfoldedSign * $fromSign * $toSign]

		if {(($expectedSign == [sign [$curPeak intensity]]) || $ignoreSign)} {
			
		    #
		    # record the PeakAssignment in the Peak
		    #
		    
		    incr totNumPeakAssignments
	    
		    set curPAName [format "%s_%d" [$curPeak name] [$curPeak numPeakAssignments]]
		    set curPA [PeakAssignment -args $curPAName]
		    $curPA setNMono $nMono
		    $curPA setAveExp $aveExp
	    
		    $curPeak addPeakAssignment [$curPA cget -this]
	
		    $curPA setFromAssignment [$fromShiftAssign cget -this]
		    $curPA setToAssignment   [$toShiftAssign   cget -this]
		    
		    $curPA setUnfoldedFromProtonPeakPosition    $unfoldedFromProtonPP
		    $curPA setUnfoldedFromHeavyatomPeakPosition $unfoldedFromHeavyPP
		    $curPA setUnfoldedToProtonPeakPosition      $unfoldedToProtonPP

		    $curPA resetNumFiltersFailed 
		    $curPA setPreviousLikelihood 1.0
		    
		    $curPA -disown
		    rename $curPA ""
		}
		rename $toShiftAssign ""
	    }
	    rename $fromShiftAssign ""
	}

	#
	# if I found no peakAssignments, complain
	#
	
	if {! [$curPeak isAssigned]} {
	    $curPeak appendToNote "No 3d assignments found."	    
	    incr numUnassignablePeaks
	} else {
	    incr numAssignedPeaks
	}

	rename $curPeak ""
    }

    if {$remVar != ""} {
	upvar $remVar tempRem

	set newRem [format "Matching %d 3-D peaks with %d from and %d to shiftAssignments" \
			[llength $peakList] [llength $fromSAs] [llength $toSAs]]

	if {$useIndividualTols} {
	    appendLineToString newRem "Using each shiftAssign's individual tolerances.  Defaults are: "
	}

	appendLineToString newRem  "from proton tolerance:     $fromProtonTol ppm"
	appendLineToString newRem  "from heavyatom tolerance:  $fromHeavyTol ppm"
	appendLineToString newRem  "to proton tolerance:       $toProtonTol ppm"

	appendLineToString newRem  [format "%d peaks were assigned with a total of %d assignments (mean degeneracy: %f)" \
					$numAssignedPeaks $totNumPeakAssignments [expr $totNumPeakAssignments / max(1.,double($numAssignedPeaks))]]
	appendLineToString newRem  [format " %d peaks (%f%% of the total) were found to be without assignments." \
					$numUnassignablePeaks [expr 100.0 * $numUnassignablePeaks / max(1.,double($numUnassignablePeaks + $numAssignedPeaks))]]
	
	lappend tempRem $newRem

	set newRem [format "%d peaks were located outside of the spectrum" [llength $peaksOutsideSpectrum]]
	if {$verboseOutsideSpectrum} {
	    foreach outsideError $peaksOutsideSpectrum {
		appendLineToString newRem $outsideError
	    }
	}

	lappend tempRem $newRem
    }

    return ""
}


proc match4d args {

    set pot               [requiredFlagVal $args -pot]
    set fromProtonTol     [requiredFlagVal $args -fromProtonTolerancePPM]
    set fromProtonRange   [requiredFlagVal $args -fromProtonSpectralRangePPM]
    set fromHeavyTol      [requiredFlagVal $args -fromHeavyatomTolerancePPM]
    set fromHeavyRange    [requiredFlagVal $args -fromHeavyatomSpectralRangePPM]
    set toProtonTol       [requiredFlagVal $args -toProtonTolerancePPM]
    set toProtonRange     [requiredFlagVal $args -toProtonSpectralRangePPM]
    set toHeavyTol        [requiredFlagVal $args -toHeavyatomTolerancePPM]
    set toHeavyRange      [requiredFlagVal $args -toHeavyatomSpectralRangePPM]
    set peakList          [flagVal $args -peakList [$pot peaks]]
    set SAList            [flagVal $args -shiftAssignList \
			       [$pot shiftAssignments]]
    set remVar            [flagVal $args -remarksVariableName ""]
    set minSALikelihood       [flagVal $args -minShiftAssignmentLikelihood -1]
    set verboseOutsideSpectrum [flagExists $args -verboseOutsideSpectrum]

    set useIndividualTols     [flagExists $args -useIndividualTolerances]

    set fromProtonSignChanges [flagExists $args [list -signChangesUponFoldingFromProtonDimension \
						     -signChangesUponAliasingFromProtonDimension]]

    set fromHeavySignChanges  [flagExists $args [list -signChangesUponFoldingFromHeavyatomDimension \
						     -signChangesUponAliasingFromHeavyatomDimension]]

    set toProtonSignChanges   [flagExists $args [list -signChangesUponFoldingToProtonDimension \
						     -signChangesUponAliasingToProtonDimension]]

    set toHeavySignChanges    [flagExists $args [list -signChangesUponFoldingToHeavyatomDimension \
						     -signChangesUponAliasingToHeavyatomDimension]]

    set foldedAlongFromProton [flagExists $args -shiftsFoldedAlongFromProtonDimension]
    set foldedAlongFromHeavy  [flagExists $args -shiftsFoldedAlongFromHeavyatomDimension]
    set foldedAlongToProton   [flagExists $args -shiftsFoldedAlongToProtonDimension]
    set foldedAlongToHeavy    [flagExists $args -shiftsFoldedAlongToHeavyatomDimension]

    marvinPyth command {import pasd}
    set defAveExp [lindex \
		    [lindex \
			 [marvinPyth command {out=pasd.aveExp} \
			      {} {out}] \
			 0] 1]
    set aveExp                [flagVal $args -aveExp $defAveExp]

    set unfoldedSign 0

    if {[flagExists $args -nonFoldedPeaksArePositive]} {
	set unfoldedSign 1
    }

    if {[flagExists $args -nonFoldedPeaksAreNegative]} {
	set unfoldedSign -1
    }

    if {$unfoldedSign == 0} {
	error "Missing either -nonFoldedPeaksArePositive or -nonFoldedPeaksAreNegative flag"
    }

    set ignoreSign [flagExists $args -ignoreSign]

    #
    # extract spectral top/bottom from range flags
    #

    set fromProtonTop [max $fromProtonRange]
    set fromProtonBot [min $fromProtonRange]
    set fromHeavyTop  [max $fromHeavyRange]
    set fromHeavyBot  [min $fromHeavyRange]
    set toProtonTop   [max $toProtonRange]
    set toProtonBot   [min $toProtonRange]
    set toHeavyTop    [max $toHeavyRange]
    set toHeavyBot    [min $toHeavyRange]

    #
    # assemble lists of from and to shiftAssigns that can be used--ie., 
    # those that have all needed chemical shifts, and appropriate previousLikelihoods
    #

    set fromSAs [list]
    set toSAs   [list]

    foreach curSA $SAList {
	ShiftAssignment -this $curSA
	if {[$curSA isFrom] && [$curSA hasProtonShift] && [$curSA hasHeavyatomShift] && ([$curSA previousLikelihood] >= $minSALikelihood)} {
	    lappend fromSAs [$curSA cget -this]
	} elseif {[$curSA isTo] && [$curSA hasProtonShift] && [$curSA hasHeavyatomShift] && ([$curSA previousLikelihood] >= $minSALikelihood)} {
	    lappend toSAs [$curSA cget -this]
	}
	rename $curSA ""
    }

    #
    # Need to record unfolded positions of each shiftAssign for use in 
    # explicit inverse exceptions.  Do it here because I have the spectral 
    # ranges and so forth
    #
    
    recordUnfoldedPositions \
	-shiftAssignments                        [concat $fromSAs $toSAs] \
	-fromProtonSpectralRange                 $fromProtonRange \
	-fromHeavyatomSpectralRange              $fromHeavyRange \
	-toProtonSpectralRange                   $toProtonRange \
	-toHeavyatomSpectralRange                $toHeavyRange \
	-fromProtonSignChanges                   $fromProtonSignChanges \
	-fromHeavyatomSignChanges                $fromHeavySignChanges \
	-toProtonSignChanges                     $toProtonSignChanges \
	-toHeavyatomSignChanges                  $toHeavySignChanges \
	-shiftsFoldedAlongFromProtonDimension    $foldedAlongFromProton \
	-shiftsFoldedAlongFromHeavyatomDimension $foldedAlongFromHeavy \
	-shiftsFoldedAlongToProtonDimension      $foldedAlongToProton \
	-shiftsFoldedAlongToHeavyatomDimension   $foldedAlongToHeavy \


    #
    # remove any linked peakAssignments that exist
    #

    $pot resetLinkedPeakAssignmentNames

    #
    # For each peak, try unfolding its position to best match all 
    # from and to shiftAssignments.  
    #
    # Generate all possible pairings of from and to ShiftAssignments 
    # whose shifts are within tolerance of the peak's unfolded position. 
    #
    # Make sure the peak's sign matches the expected value (which depends 
    # on how many times it was unfolded to match the shiftAssignments).  
    #
    # Record those pairings of ShiftAssignments as assignments
    #

    set numPeaks [llength $peakList]
    set curPeakNum 0
    set numUnassignablePeaks 0
    set numAssignedPeaks 0
    set totNumPeakAssignments 0
    set peaksOutsideSpectrum [list]

    marvinPyth command {import pasd}
    set nMono [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.nMono} \
			   {} {out}] \
		      0] 1]

    foreach curPeak $peakList {

	Peak -this $curPeak

	updateUser [format "Checking for matches against peak %s (%d of %d) \r" \
			[$curPeak name] [incr curPeakNum] $numPeaks]

	#
	# remove any existing shift and PeakAssignments from the current peak
	#

	$curPeak removeAllPeakAssignments
	$curPeak removeAllShiftAssignments
	
	#
	# check that this peak's position is inside the spectrum
	#

	if {[peakOutsideSpectrum [$curPeak fromProtonShift] $fromProtonTop $fromProtonBot]} {

	    set err [format "from proton location at %f is outside of spectrum (%f .. %f)" \
			 [$curPeak fromProtonShift] $fromProtonTop $fromProtonBot]


	    $curPeak appendToNote $err
	    
	    updateUser [format "Peak %s %s \n" [$curPeak name] $err]

	    lappend peaksOutsideSpectrum [format "Peak %s %s" [$curPeak name] $err]
	}

	if {[peakOutsideSpectrum [$curPeak fromHeavyatomShift] $fromHeavyTop $fromHeavyBot]} {

	    set err [format "from heavyatom location at %f is outside of spectrum (%f .. %f)" \
			 [$curPeak fromHeavyatomShift] $fromHeavyTop $fromHeavyBot]

	    $curPeak appendToNote $err
	    
	    updateUser [format "Peak %s %s \n" [$curPeak name] $err]
	    
	    lappend peaksOutsideSpectrum [format "Peak %s %s" [$curPeak name] $err]
	}


	if {[peakOutsideSpectrum [$curPeak toProtonShift] $toProtonTop $toProtonBot]} {

	    set err [format "to proton location at %f is outside of spectrum (%f .. %f)" \
			 [$curPeak toProtonShift] $toProtonTop $toProtonBot]

	    $curPeak appendToNote $err
	    
	    updateUser [format "Peak %s %s \n" [$curPeak name] $err]

	    lappend peaksOutsideSpectrum [format "Peak %s %s" [$curPeak name] $err]
	}

	if {[peakOutsideSpectrum [$curPeak toHeavyatomShift] $toHeavyTop $toHeavyBot]} {

	    set err [format "to heavyatom location at %f is outside of spectrum (%f .. %f)" \
			 [$curPeak toHeavyatomShift] $toHeavyTop $toHeavyBot]

	    $curPeak appendToNote $err
	    
	    updateUser [format "Peak %s %s \n" [$curPeak name] $err]

	    lappend peaksOutsideSpectrum [format "Peak %s %s" [$curPeak name] $err]
	}

	#
	# find from shiftAssigns that match the restraint's from dimensions.
	# Record each one's name, expected sign, and unfolded peak positions
	#


	set matchingFromShiftAssigns [list]

	foreach fromShiftAssign $fromSAs {	    
	    ShiftAssignment -this $fromShiftAssign

	    if {$useIndividualTols && [$fromShiftAssign hasProtonTolerance]} {
		set fpt [$fromShiftAssign protonTolerance]
	    } else {
		set fpt $fromProtonTol
	    }

	    if {$useIndividualTols && [$fromShiftAssign hasHeavyatomTolerance]} {
		set fht [$fromShiftAssign heavyatomTolerance]
	    } else {
		set fht $fromHeavyTol
	    }

	    if {(! [$fromShiftAssign hasProtonShift]) || (! [$fromShiftAssign hasHeavyatomShift])} {
		rename $fromShiftAssign ""
		continue
	    }

	    set unfoldedFromProtonPPs [unfoldPeakPosnToMatch [$curPeak fromProtonShift]    [$fromShiftAssign protonShift]    $foldedAlongFromProton $fromProtonTop $fromProtonBot $fromProtonSignChanges] 
	    set unfoldedFromHeavyPPs  [unfoldPeakPosnToMatch [$curPeak fromHeavyatomShift] [$fromShiftAssign heavyatomShift] $foldedAlongFromHeavy  $fromHeavyTop  $fromHeavyBot  $fromHeavySignChanges] 

	    foreach protonPP $unfoldedFromProtonPPs {
		if {[unfoldedPeakPosnMatches $protonPP [$fromShiftAssign protonShift] $fpt]} {
		    foreach heavyPP $unfoldedFromHeavyPPs {
			if {[unfoldedPeakPosnMatches $heavyPP [$fromShiftAssign heavyatomShift] $fht]} {
			    
			    if {! [$curPeak hasShiftAssignmentNamed [$fromShiftAssign name]]} {
				$curPeak addShiftAssignment [$fromShiftAssign cget -this]
			    }

			    set curMatchSign [expr [lindex $protonPP 1] * [lindex $heavyPP 1]]
			    lappend matchingFromShiftAssigns [list [$fromShiftAssign name] $curMatchSign [lindex $protonPP 0] [lindex $heavyPP 0]]
			}
		    }
		}
	    }
	    
	    rename $fromShiftAssign ""
	}

		
	#
	# find to shiftAssigns that match the restraint's to dimensions
	# Record each one's name, expected sign, and unfolded peak positions
	#


	set matchingToShiftAssigns [list]

	foreach toShiftAssign $toSAs {	    
	    ShiftAssignment -this $toShiftAssign

	    if {$useIndividualTols && [$toShiftAssign hasProtonTolerance]} {
		set tpt [$toShiftAssign protonTolerance]
	    } else {
		set tpt $toProtonTol
	    }

	    if {$useIndividualTols && [$toShiftAssign hasHeavyatomTolerance]} {
		set tht [$toShiftAssign heavyatomTolerance]
	    } else {
		set tht $toHeavyTol
	    }

	    if {(! [$toShiftAssign hasProtonShift]) || (! [$toShiftAssign hasHeavyatomShift])} {
		rename $toShiftAssign ""
		continue
	    }

	    set unfoldedToProtonPPs [unfoldPeakPosnToMatch [$curPeak toProtonShift]    [$toShiftAssign protonShift]    $foldedAlongToProton $toProtonTop $toProtonBot $toProtonSignChanges] 
	    set unfoldedToHeavyPPs  [unfoldPeakPosnToMatch [$curPeak toHeavyatomShift] [$toShiftAssign heavyatomShift] $foldedAlongToHeavy  $toHeavyTop  $toHeavyBot  $toHeavySignChanges] 

	    foreach protonPP $unfoldedToProtonPPs {
		if {[unfoldedPeakPosnMatches $protonPP [$toShiftAssign protonShift] $tpt]} {
		    foreach heavyPP $unfoldedToHeavyPPs {
			if {[unfoldedPeakPosnMatches $heavyPP [$toShiftAssign heavyatomShift] $tht]} {
			    
			    if {! [$curPeak hasShiftAssignmentNamed [$toShiftAssign name]]} {
				$curPeak addShiftAssignment [$toShiftAssign cget -this]
			    }

			    set curMatchSign [expr [lindex $protonPP 1] * [lindex $heavyPP 1]]
			    lappend matchingToShiftAssigns [list [$toShiftAssign name] $curMatchSign [lindex $protonPP 0] [lindex $heavyPP 0]]
			}
		    }
		}
	    }
	    
	    rename $toShiftAssign ""
	}

	#
	# now generate all possible pairings of from and to shiftAssigns that match the 
	# current restraint, check that their expected signs match the peak's, 
	# and generate corresponding PeakAssignments
	#

	foreach fromMatch $matchingFromShiftAssigns {
	    set fromShiftAssign      [$pot shiftAssignmentNamed [lindex $fromMatch 0]]
	    ShiftAssignment -this $fromShiftAssign

	    set fromSign             [lindex $fromMatch 1]
	    set unfoldedFromProtonPP [lindex $fromMatch 2]
	    set unfoldedFromHeavyPP  [lindex $fromMatch 3]

	    foreach toMatch $matchingToShiftAssigns {
		set toShiftAssign      [$pot shiftAssignmentNamed [lindex $toMatch 0]]
		ShiftAssignment -this $toShiftAssign

		set toSign             [lindex $toMatch 1]
		set unfoldedToProtonPP [lindex $toMatch 2]
		set unfoldedToHeavyPP  [lindex $toMatch 3]

		set expectedSign [expr $unfoldedSign * $fromSign * $toSign]

		if {(($expectedSign == [sign [$curPeak intensity]]) || $ignoreSign)} {
			
		    #
		    # record the PeakAssignment in the peak
		    #
			
		    incr totNumPeakAssignments
			
		    set curPAName [format "%s_%d" [$curPeak name] [$curPeak numPeakAssignments]]
		    set curPA [PeakAssignment -args $curPAName]
		    $curPA setNMono $nMono
		    $curPA setAveExp $aveExp
			
		    $curPeak addPeakAssignment [$curPA cget -this]
			
		    $curPA setFromAssignment [$fromShiftAssign cget -this]
		    $curPA setToAssignment   [$toShiftAssign   cget -this]

		    $curPA setUnfoldedFromProtonPeakPosition    $unfoldedFromProtonPP
		    $curPA setUnfoldedFromHeavyatomPeakPosition $unfoldedFromHeavyPP
		    $curPA setUnfoldedToProtonPeakPosition      $unfoldedToProtonPP
		    $curPA setUnfoldedToHeavyatomPeakPosition   $unfoldedToHeavyPP

		    $curPA resetNumFiltersFailed 
		    $curPA setPreviousLikelihood 1.0

		    $curPA -disown
		    rename $curPA ""
		}
		rename $toShiftAssign ""
	    }
	    rename $fromShiftAssign ""
	}

	#
	# if I found no assignments, complain
	#
	
	if {! [$curPeak isAssigned]} {
	    $curPeak appendToNote "No 4d assignments found."	    
	    incr numUnassignablePeaks
	} else {
	    incr numAssignedPeaks
	}
	
	rename $curPeak ""
    }
    
    if {$remVar != ""} {
	upvar $remVar tempRem

	set newRem [format "Matching %d 4-D peaks with %d from and %d to shiftAssignments" \
			[llength $peakList] [llength $fromSAs] [llength $toSAs]]

	if {$useIndividualTols} {
	    appendLineToString newRem "Using each shiftAssign's individual tolerances.  Defaults are: "
	}

	appendLineToString newRem  "from proton tolerance:     $fromProtonTol ppm"
	appendLineToString newRem  "from heavyatom tolerance:  $fromHeavyTol ppm"
	appendLineToString newRem  "to proton tolerance:       $toProtonTol ppm"
	appendLineToString newRem  "to heavyatom tolerance:    $toHeavyTol ppm"

	appendLineToString newRem  [format "%d peaks were assigned with a total of %d assignments (mean degeneracy: %f)" \
					$numAssignedPeaks $totNumPeakAssignments [expr $totNumPeakAssignments / max(1.,double($numAssignedPeaks))]]
	appendLineToString newRem  [format " %d peaks (%f%% of the total) were found to be without assignments." \
					$numUnassignablePeaks [expr 100.0 * $numUnassignablePeaks / max(1.,double($numUnassignablePeaks + $numAssignedPeaks))]]
	
	lappend tempRem $newRem
	
	set newRem [format "%d peaks were located outside of the spectrum" [llength $peaksOutsideSpectrum]]
	if {$verboseOutsideSpectrum} {
	    foreach outsideError $peaksOutsideSpectrum {
		appendLineToString newRem $outsideError
	    }
	}

	lappend tempRem $newRem
    }

    return ""
}


#
# FIX: Note that this is incomplete.  It doesn't work if the diagonal has been folded
#

proc removeDiagonalPeaks args {

    set pot           [requiredFlagVal $args -pot]
    set peakList      [flagVal $args -peakList]
    set protTol       [flagVal $args -tolerance 0.001]
    set heavyTol      [flagVal $args -heavyatomTolerance -1]
    set remVar        [flagVal $args -remarksVariableName ""]

    if {$peakList == ""} {
	set peakList [$pot peaks]
    }

    set nRemoved 0
    set nStart [llength $peakList]

    set invsqrt2 [expr 1.0 / double(sqrt(2.0))]

    set rCount 0

    foreach p $peakList {

	Peak -this $p

	updateUser [format "Checking peak %s (%d of %d) \r" \
			[$p name] [incr rCount] [llength $peakList]]

	# compute distance from (Hi,Hj) to diagonal.  See 
	# ATNOS paper, eqn 21.  Their default cutoff is 
	# 0.7 ppm for 2D NOESYs, 0.6 ppm for 3D.  

	set diagDist [expr abs([$p fromProtonShift] - [$p toProtonShift]) * $invsqrt2]

	if {$diagDist < $protTol} {

	    set heavyOK 1
	    if {[$p hasFromHeavyatomShift] && [$p hasToHeavyatomShift] && ($heavyTol != -1) } {

		set heavyDiagDist [expr abs([$p fromHeavyatomShift] - [$p toHeavyatomShift]) * $invsqrt2]		
		if {$heavyDiagDist >= $heavyTol} {
		    set heavyOK 0
		}
	    }
	    
	    if {$heavyOK} {
		$pot removePeakNamed [$p name]
		$p -acquire
		incr nRemoved
	    }
	}

	rename $p ""
    }

    if {$remVar != ""} {
	upvar $remVar tempRem
	if {$heavyTol == -1} {
	    lappend tempRem [format "%d of %d peaks were dropped because of proximity < %f ppm to the diagonal." \
				 $nRemoved $nStart $protTol]
	} else {
	    lappend tempRem [format "%d of %d peaks were dropped because of proximity < %f ppm (proton) and %f ppm (heavyatom) to the diagonal." \
				 $nRemoved $nStart $protTol $heavyTol]
	}
    }
    return ""
}

proc findSolventLine args {

    set anNOEpotential  [flagVal $args -pot]
    set aPeakList       [flagVal $args -peakList]
    set range           [flagVal $args -shiftRange [list 4.5 4.9]]
    set binWidth        [flagVal $args -binWidth 0.01]
    set cutoff          [flagVal $args -cutoff 0.20]
    set minNumPeaks     [flagVal $args -minNumPeaks 500]
    set remVar          [flagVal $args -remarksVariableName ""]
    set useTo           [flagExists $args -toProtonDimension]
    
    if {($anNOEpotential == "") && ($aPeakList == "")} {
	error "Neither -pot or -peakList defined in call to findSolventLine"
    }
    
    if {$anNOEpotential != ""} {
	set aPeakList [$anNOEpotential peaks]
    }
    
    set minVal [min $range]
    set maxVal [max $range]
    
    #
    # gather proton shifts that are in the relevant range
    #
    
    set shiftsInRange [list]
    foreach p $aPeakList {
	
	Peak -this $p
	
	if {$useTo} {
	    set curVal [$p toProtonShift]
	} else {
	    set curVal [$p fromProtonShift]
	}
	
	if {($curVal > $minVal) && ($curVal < $maxVal)} {
	    lappend shiftsInRange $curVal
	}
	
	rename $p ""
    }

    #
    # make sure there's a reasonable number of peaks in 
    # the relevant range
    #

    if {[llength $shiftsInRange] < $minNumPeaks} {

	error "No solvent line detected."
    }

    
    #
    # determine fraction of shifts in each bin
    #
    
    set nBins [expr ceil(($maxVal - $minVal) / double($binWidth))]
    
    for {set binCount 0} {$binCount < $nBins} {incr binCount} {
	set numInBin($binCount) 0
    }
    
    foreach shift $shiftsInRange {
	incr numInBin([whichBin $minVal $maxVal $nBins $shift])
    }
    
    set numShiftsInRange [llength $shiftsInRange]
    
    for {set binCount 0} {$binCount < $nBins} {incr binCount} {
	set fracInBin($binCount) [expr $numInBin($binCount) / double($numShiftsInRange)]
    }
    
    #
    # find the max bin
    #
    
    set maxFrac -1
    
    for {set binCount 0} {$binCount < $nBins} {incr binCount} {
	if {$fracInBin($binCount) > $maxFrac} {
	    set maxFrac $fracInBin($binCount)
	    set maxBin $binCount
	}
    }

    set binRemark [format "Largest bin has %f%% of the total %d peaks in the range %f .. %f ppm" \
		       [expr $maxFrac * 100.0] $numShiftsInRange $minVal $maxVal]
    
    #
    # make sure its fraction is high enough
    #
    
    if {$maxFrac < $cutoff} {
	
	if {$remVar != ""} {
	    
	    if {$useTo} {
		set newRem "No solvent line detected in to proton dimension"
	    } else {
		set newRem "No solvent line detected in from proton dimension"
	    }

	    appendLineToString newRem $binRemark
	    
	    upvar $remVar tempRem
	    lappend tempRem $newRem
	}

	error "No solvent line detected."
    }
    
    #
    # return the center of the max bin
    #
    
    set actBinWidth [expr ($maxVal - $minVal) / double($nBins)]
    set center [expr $minVal + ($maxBin + 0.5) * $actBinWidth]
    
    if {$remVar != ""} {

	if {$useTo} {
	    set newRem "Solvent line detected at $center ppm in to proton dimension"
	} else {
	    set newRem "Solvent line detected at $center ppm in from proton dimension"
	}
	appendLineToString newRem $binRemark

	upvar $remVar tempRem
	lappend tempRem $newRem
    }
    
    return $center
}




proc removeSolventPeaks args {

    set peakList       [requiredFlagVal $args -peakList]
    set pot            [requiredFlagVal $args -pot]
    set tol            [flagVal $args -tolerance 0.02]
    set cutoff         [flagVal $args -cutoff 0.20]
    set autoRange      [flagVal $args -autoDetectRange [list 4.5 4.9]]
    set solShift       [flagVal $args -solventShift 4.76]
    set remVar         [flagVal $args -remarksVariableName ""]
    set useTo          [flagExists $args -toProtonDimension]

    #
    # I need to move this here to allow remarks to be written by 
    # the findSolventLine proc
    #

    if {$remVar != ""} {
	upvar $remVar tempRem
    }

    global errorInfo

    #
    # try to find the solvent line automatically.  
    # Throws an error if it can't find one.
    #

    if {[flagExists $args -autoDetect]} {

	if {$useTo} {

	    if {[catch {set solShift [findSolventLine -peakList $peakList \
					  -shiftRange $autoRange \
					  -cutoff $cutoff \
					  -toProtonDimension \
					  -remarksVariableName tempRem]}]} {
		set errorInfo ""
		return 
	    }

	} else {

	    
	    if {[catch {set solShift [findSolventLine -peakList $peakList \
					  -shiftRange $autoRange \
					  -cutoff $cutoff \
					  -remarksVariableName tempRem]}]} {

		set errorInfo ""
		return
	    }
	}

	updateUser "Found solvent line at $solShift ppm.\n"
    }
    
    #
    # eliminate all peaks in the solvent band
    #

    
    set nStart [llength $peakList]
    set nRemoved 0

    foreach p $peakList {

	Peak -this $p

	if {$useTo} {
	    set curShift [$p toProtonShift]
	} else {
	    set curShift [$p fromProtonShift]
	}

	set delta [expr abs($curShift - $solShift)]
	if {$delta < $tol} {

	    # cur peak is within solvent band.
	    
	    $pot removePeakNamed [$p name]
	    $p -acquire
	    incr nRemoved
	}

	rename $p ""
    }

    if {$remVar != ""} {

	set newRem ""

	if {$useTo} {
	    appendLineToString newRem "Checking to proton dimension for solvent peaks."
	} else {
	    appendLineToString newRem "Checking from proton dimension for solvent peaks."
	}
	
	appendLineToString newRem \
	    [format "%d of %d peaks were dropped because of proximity to the solvent (%f +/- %f ppm)." \
		 $nRemoved $nStart $solShift $tol]
	
	lappend tempRem $newRem
    }
    return ""
}




proc generateDistanceBounds args {

    set peakList  [requiredFlagVal $args -peakList]

    set intensityList [list]
    foreach p $peakList {

	Peak -this $p
	lappend intensityList [$p intensity]
	rename $p ""
    }

    set distanceBins [getDistanceBinsFromIntensities $intensityList]
       
    foreach p $peakList {

        Peak -this $p

	set bounds [estimateBounds [$p intensity] $distanceBins]
	$p setLowBound [lindex $bounds 0]
	$p setUpBound [lindex $bounds 1]

	rename $p ""
    }
}



proc correctUpBoundsForMethyls args {

    #
    # This isn't necessary.  Could just use the methyl flag on each 
    # shiftAssignment and the methyl correction defined for the inverse
    # restraints in MarvinNOEPotential
    #

    set anNOEpotential  [flagVal $args -pot]
    set aPeakList       [flagVal $args -peakList]
    
    if {($anNOEpotential == "") && ($aPeakList == "")} {
	error "Neither -pot or -peakList defined in call to correctUpBoundsForMethyls"
    }
    
    if {$anNOEpotential != ""} {
	set aPeakList [$anNOEpotential peaks]
    }
    
    foreach curPeak $aPeakList {
	
	Peak -this $curPeak
	
	foreach curPA [$curPeak peakAssignments] {
	    
	    PeakAssignment -this $curPA

	    set fromAssign [$curPA fromAssignment]
	    set toAssign   [$curPA toAssignment]
	    
	    ShiftAssignment -this $fromAssign
	    ShiftAssignment -this $toAssign

	    set nMethyls [expr [$fromAssign isMethyl] + [$toAssign isMethyl]]

	    if {$nMethyls > 0} {
		
		set correct [expr [$curPA upBoundCorrection] + (0.5 * $nMethyls)]
		$curPA setUpBoundCorrection $correct
		
		$curPA appendToNote \
		    [format "upper bound increased because selections involve %d methyls" $nMethyls]
	    }

	    rename $fromAssign ""
	    rename $toAssign ""
	    rename $curPA ""
	}

	rename $curPeak ""
    }
}


#
# Given a peak's intensity and a list 
# of (min intensity, upper bound) pairs,
# return the proper lower and upper bounds
#

proc estimateBounds {intensity distanceBins} {

    set intensity [expr abs($intensity)]

    foreach bin $distanceBins {

	if {$intensity > [lindex $bin 0]} {
	    return [lindex $bin 1]
	}
    }

    # we only get here if this intensity is smaller 
    # than the smallest bin.

    error "Intensity too small to estimate bound"
}


proc getDistanceBinsFromIntensities {a} {

    #
    # this can be called with either the name of a peak array
    # (for use in processPippAssignedNOEFile) or with a list
    # of intensities (for use with AENEAS)
    #

    if {[llength $a] == 1} {

	upvar $a peaks
	
	set intensities [list]
	
	foreach peakID [array names peaks] {
	    
	    set curPeak $peaks($peakID)
	    set firstAssign [lindex $curPeak 0]
	    set curIntensity [lindex $firstAssign 1]
	    lappend intensities [expr abs($curIntensity)]
	}

    } else {

	set intensities [list]
	
	foreach curIntensity $a {
	    lappend intensities [expr abs ($curIntensity)]
	}
    }
    
    set sortedIntensities [lsort -real -decreasing $intensities]
    
    # 
    # uses Dan Garrett's "algorithm" for now:
    # 1.  Sort intensities list
    # 2.  largest 20% of peaks become "strong" (1.8 - 2.7)
    # 3.  next 30% become "medium" (1.8 - 3.3)
    # 4.  next 30% into "weak" (1.8 - 5.0)
    # 5.  last 20% into "very weak" (1.8 - 6.0)
    #

    set numIntensities [llength $sortedIntensities]
    
    set lowestStrongPeak   [lindex $sortedIntensities \
				[expr round($numIntensities * 0.2)]]
    set lowestMedPeak      [lindex $sortedIntensities \
				[expr round($numIntensities * 0.5)]]
    set lowestWeakPeak     [lindex $sortedIntensities \
				[expr round($numIntensities * 0.8)]]
    set lowestVeryWeakPeak [lindex $sortedIntensities $numIntensities]

    marvinPyth command "import pasd"
    set b0 [marvinPyth command {o1,o2=pasd.distanceBins[0]} {} {o1 o2}]
    set b1 [marvinPyth command {o1,o2=pasd.distanceBins[1]} {} {o1 o2}]
    set b2 [marvinPyth command {o1,o2=pasd.distanceBins[2]} {} {o1 o2}]
    set b3 [marvinPyth command {o1,o2=pasd.distanceBins[3]} {} {o1 o2}]

    set b00 [lindex [lindex $b0 0] 1]
    set b01 [lindex [lindex $b0 1] 1]
    set b10 [lindex [lindex $b1 0] 1]
    set b11 [lindex [lindex $b1 1] 1]
    set b20 [lindex [lindex $b2 0] 1]
    set b21 [lindex [lindex $b2 1] 1]
    set b30 [lindex [lindex $b3 0] 1]
    set b31 [lindex [lindex $b3 1] 1]

    set retVal [list]
    lappend retVal [list $lowestStrongPeak   [list $b00 $b01]]
    lappend retVal [list $lowestMedPeak      [list $b10 $b11]]
    lappend retVal [list $lowestWeakPeak     [list $b20 $b21]]
    lappend retVal [list $lowestVeryWeakPeak [list $b30 $b31]]

    return $retVal
}









#
# Given a list of shift entries, find atoms with no assignment
#

proc reportUnassignedAtoms args {

    set shiftList           [requiredFlagVal $args -shiftList]
    set flexibleSelection   [flagVal $args -flexibleRegionSelection ""]
    set remVar              [flagVal $args -remarksVariableName ""]

    marvinPyth command "import pasd"
    marvinPyth command "localShifts=pasd.convertShiftsFromTCL('$shiftList')"
    set cmd "ret=pasd.findUnassignedAtoms(localShifts,
                                 flexibleRegion=flexibleSelection,
                                  tclOutput=True)"
    set ret       [marvinPyth command $cmd \
			   [list [list flexibleSelection $flexibleSelection]] \
			   {ret}]
    set msg [lindex [lindex \
		      [marvinPyth command $cmd \
			   [list [list flexibleSelection $flexibleSelection]] \
			   {ret}] \
			 0] 1]
    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem $msg
    }

    if {$flexibleSelection != ""} {
	rename $flexibleSelection ""
    }
}



#
# Procs for creating ShiftAssignments
#


proc createShiftAssignments args {

    set shiftList              [requiredFlagVal $args -shiftList]
    set pot                    [requiredFlagVal $args -pot]
    set fromProtonSelStr       [requiredFlagVal $args \
				    -fromProtonSelectionString]
    set toProtonSelStr         [requiredFlagVal $args \
				    -toProtonSelectionString]
    set flexibleRegionSelStr   [flagVal $args \
				    -flexibleRegionSelectionString "NONE"]
    set fromHeavySelStr        [flagVal $args \
				    -fromHeavyatomSelectionString "NONE"]
    set toHeavySelStr          [flagVal $args \
				    -toHeavyatomSelectionString "NONE"]
    set fromProtonSolventRange [flagVal $args \
				    -fromProtonDimensionSolventRange "NONE"]
    set toProtonSolventRange   [flagVal $args \
				    -toProtonDimensionSolventRange "NONE"]
    set namePrefix             [flagVal $args -namePrefix ""]
    set remVar                 [flagVal $args -remarksVariableName ""]

    #
    # Convert all the shift list entries from strings to AtomSels
    #
    # ONLY CONVERT THE FIRST SELECTION OF EACH ENTRY!  because with
    # PIPP, I only get multiple selections from HB1 | HB2 style
    # entries.  Since these always come in pairs, I end up doubling
    # the number of shift assigns from these non-stereo methylene
    # entries
    #
    # This will eventually break processing of more sophisticated
    # shift tables' data!
    #

    #
    # complain both here and in the report about empty atom selections
    #

    updateUser "Converting shift list to atomSels \r"
    
    set badSelList [list]
    set shiftSelectionList [list]

    foreach entry $shiftList {
	set curShift [lindex $entry 0]
	set curAtomSel [AtomSel -args [lindex $entry 1]]
	if { $curShift == {} } {
	    updateUser \
		[format "Shift table entry for %s specifies no chemical shift value\n" \
		     [$curAtomSel string]]
	} else {
	    if {[$curAtomSel size] == 0} {
		lappend badSelList [list $curShift [$curAtomSel string]]
		updateUser [format "Shift table entry at %f ppm selects no atoms: %s \n" $curShift [$curAtomSel string]]
	    } else {
		lappend shiftSelectionList [list $curShift $curAtomSel]
	    }
	}
    }
 
    #
    # go through the shiftList now, caching all the shiftEntries
    # that match each chemical type of interest
    #

    updateUser "Selecting shift entries that match each atom type \r"
    
    set fromProtonShiftList [shiftEntriesThatMatchChemicalType \
				 $shiftSelectionList $fromProtonSelStr]
    set toProtonShiftList   [shiftEntriesThatMatchChemicalType \
				 $shiftSelectionList $toProtonSelStr]
    
    if {$fromHeavySelStr != "NONE"} {
	set fromHeavyShiftList  [shiftEntriesThatMatchChemicalType $shiftSelectionList $fromHeavySelStr]
    } 
    
    if {$toHeavySelStr != "NONE"} {
	set toHeavyShiftList  [shiftEntriesThatMatchChemicalType $shiftSelectionList $toHeavySelStr]
    }
    

    updateUser "Caching bondArray \r"
    cacheBondArray theBondArray
    
    #
    # generate the from and to ShiftAssignments
    #

    set numShiftAssignsInSolventRange 0
    set assCount -1

    updateUser "Generating from ShiftAssignments \r"

    set nFromShiftAssignments 0

    foreach ps $fromProtonShiftList {
	set protonShift [lindex $ps 0]
	set protonSel   [lindex $ps end]
	set compressedProtonSelString [compressedSelectionString [$protonSel string]]

	if {$fromProtonSolventRange != "NONE"} {

	    set top [max $fromProtonSolventRange]
	    set bot [min $fromProtonSolventRange]

	    if {($protonShift < $top) && ($protonShift > $bot)} {
		incr numShiftAssignsInSolventRange
		continue
	    }
	}

	if {$fromHeavySelStr == "NONE"} {
	    
	    set assName [format "%s%d_from%s" $namePrefix [incr assCount] $compressedProtonSelString]
	    set newAss  [ShiftAssignment -args $assName]

	    updateUser "Generating shiftAssignment $assName \r"
				
	    $newAss setProtonSelectionString    [$protonSel string]
	    $newAss setProtonShift              $protonShift
	    $newAss setIsFrom
	    $newAss setPreviousLikelihood       1.0
	    
	    $pot addShiftAssignment [$newAss cget -this]
	    $newAss -disown
	    rename $newAss ""

	    incr nFromShiftAssignments
	    
	} else {

	    foreach hs $fromHeavyShiftList {
		set heavyShift [lindex $hs 0]
		set heavySel   [lindex $hs end]
		
		if {[selectionsAreBonded $protonSel $heavySel theBondArray]} {
		    
		    set assName [format "%s%d_from%s" \
				     $namePrefix [incr assCount] \
				     $compressedProtonSelString]
		    set newAss  [ShiftAssignment -args $assName]

		    updateUser "Generating shiftAssignment $assName \r"
		    
		    $newAss setProtonSelectionString    [$protonSel string]
		    $newAss setHeavyatomSelectionString [$heavySel string]
		    $newAss setProtonShift              $protonShift
		    $newAss setHeavyatomShift           $heavyShift
		    $newAss setIsFrom
		    $newAss setPreviousLikelihood       1.0
		    
		    $pot addShiftAssignment [$newAss cget -this]
		    $newAss -disown
		    rename $newAss ""

		    incr nFromShiftAssignments

		}
	    }
	}
    }

    
    updateUser "Generating to ShiftAssignments \r"

    set nToShiftAssignments 0

    foreach ps $toProtonShiftList {
	set protonShift [lindex $ps 0]
	set protonSel   [lindex $ps end]
	set compressedProtonSelString [compressedSelectionString [$protonSel string]]

	if {$toProtonSolventRange != "NONE"} {

	    set top [max $toProtonSolventRange]
	    set bot [min $toProtonSolventRange]

	    if {($protonShift < $top) && ($protonShift > $bot)} {
		incr numShiftAssignsInSolventRange
		continue
	    }
	}

	if {$toHeavySelStr == "NONE"} {
	    
	    set assName [format "%s%d_to%s" $namePrefix [incr assCount] $compressedProtonSelString]
	    set newAss  [ShiftAssignment -args $assName]

	    updateUser "Generating shiftAssignment $assName \r"
				
	    $newAss setProtonSelectionString    [$protonSel string]
	    $newAss setProtonShift              $protonShift
	    $newAss setIsTo
	    $newAss setPreviousLikelihood       1.0
	    
	    $pot addShiftAssignment [$newAss cget -this]
	    $newAss -disown
	    rename $newAss ""

	    incr nToShiftAssignments 
	    
	} else {

	    foreach hs $toHeavyShiftList {
		set heavyShift [lindex $hs 0]
		set heavySel   [lindex $hs end]
		
		if {[selectionsAreBonded $protonSel $heavySel theBondArray]} {
		    
		    set assName [format "%s%d_to%s" $namePrefix [incr assCount] $compressedProtonSelString]
		    set newAss  [ShiftAssignment -args $assName]

		    updateUser "Generating shiftAssignment $assName \r"
		    
		    $newAss setProtonSelectionString    [$protonSel string]
		    $newAss setHeavyatomSelectionString [$heavySel string]
		    $newAss setProtonShift              $protonShift
		    $newAss setHeavyatomShift           $heavyShift
		    $newAss setIsTo
		    $newAss setPreviousLikelihood       1.0
		    
		    $pot addShiftAssignment [$newAss cget -this]
		    $newAss -disown
		    rename $newAss ""

		    incr nToShiftAssignments 
		}
	    }
	}
    }

    #
    # remove shiftAssigns whose selections overlap the flexible region selection
    #

    if {$flexibleRegionSelStr != "NONE"} {
	set flexibleSelection [AtomSel -args $flexibleRegionSelStr]
	set SAsToRemove [list]

	foreach curSA [$pot shiftAssignments] {
	    ShiftAssignment -this $curSA
	    
	    if {[$flexibleSelection intersects [$curSA protonSelection]]} {
		lappend SAsToRemove [$curSA cget -this]
	    }
	    
	    if {[$curSA hasHeavyatomSelection]} {
		if {[$flexibleSelection intersects [$curSA heavyatomSelection]]} {
		    lappend SAsToRemove [$curSA cget -this]
		}
	    }

	    rename $curSA ""
	}

	set SAsToRemove [lsort -unique $SAsToRemove]
	set nRemovedForFlex [llength $SAsToRemove]

	foreach curSA $SAsToRemove {
	    ShiftAssignment -this $curSA
	    $pot removeShiftAssignmentNamed [$curSA name]
	    $curSA -acquire
	    rename $curSA ""
	}
    }

    #
    # clean up
    #

    array unset theBondArray
    
    foreach entry $shiftSelectionList {
	foreach sel [lrange $entry 1 end] {
	    rename $sel ""
	}
    }

    if {$remVar != ""} {
	upvar $remVar tempRem

	if {[llength $badSelList] > 0} {
	    set rpt [format "%d shift table entries select zero atoms:" [llength $badSelList]]
	    foreach elem $badSelList {
		appendLineToString rpt [format "%f ppm, selection = %s" [lindex $elem 0] [lindex $elem 1]]
	    }
	    lappend tempRem $rpt
	}

	set rpt [format "%d from shiftAssigns were created:" $nFromShiftAssignments]
	appendLineToString rpt [format "   protons selected by %s" $fromProtonSelStr]
	if {$fromHeavySelStr != "NONE"} {
	    appendLineToString rpt [format "   heavyatoms selected by %s" $fromHeavySelStr]
	}

	lappend tempRem $rpt

	set rpt [format "%d to shiftAssigns were created:" $nToShiftAssignments]
	appendLineToString rpt [format "   protons selected by %s" $toProtonSelStr]
	if {$toHeavySelStr != "NONE"} {
	    appendLineToString rpt [format "   heavyatoms selected by %s" $toHeavySelStr]
	}

	lappend tempRem $rpt

	set rpt [format "Solvent range for from proton dimension was %s" $fromProtonSolventRange]
	appendLineToString rpt [format "Solvent range for to proton dimension was %s" $toProtonSolventRange]
	appendLineToString rpt [format "%d shiftAssignments were within one or the other range and were removed" \
				     $numShiftAssignsInSolventRange]

	lappend tempRem $rpt

	if {$flexibleRegionSelStr != "NONE"} {
	    set rpt [format "Known flexible region is %s" $flexibleRegionSelStr]
	    appendLineToString rpt [format "%d shiftAssignments were removed for having proton or heavyatom selections that overlapped with that region" $nRemovedForFlex]
	    lappend tempRem $rpt
	}
    }
    return ""
}



#
# a quick check for equivalent peakAssigns
#

proc findEquivPeakAssigns args {

    set pot [requiredFlagVal $args -pot] 

    foreach curPeak [$pot peaks] {
	Peak -this $curPeak
	foreach curPA [$curPeak peakAssignments] {
	    PeakAssignment -this $curPA
	    set names [join [list [$curPA fromAssignmentName] [$curPA toAssignmentName]]]
	    lappend SAnames($names) [$curPA name]
	    rename $curPA ""
	} 
	rename $curPeak ""
    }

    set nEquivSets 0
    set nRedundantPAs 0

    foreach n [array names SAnames] {
	if {[llength $SAnames($n)] > 1} {
	    set fromSA [lindex [split $n] 0]
	    set toSA   [lindex [split $n] 1]
	    puts [format "Equivalent peak Assigns using SAs %s" $n]
	    foreach elem $SAnames($n) {
		puts [format "   peakAssign %s" $elem]
	    }
	    incr nEquivSets
	    incr nRedundantPAs [expr [llength $SAnames($n)] -1]
	}
    }

    puts [format "total number of sets of equivalent PAs: %d" $nEquivSets]
    puts [format "total number of redundant PAs: %d" $nRedundantPAs]
}


proc initializePeakAssignmentNumFiltersFailed args {

    set pot [requiredFlagVal $args -pot]

    foreach curPeak [$pot peaks] {
	Peak -this $curPeak
	foreach curPA [$curPeak peakAssignments] {
	    PeakAssignment -this $curPA
	    $curPA resetNumFiltersFailed
	    rename $curPA ""
	}
	rename $curPeak ""
    }
}







proc createExplicitInverseExceptions args {

    set pot                    [requiredFlagVal $args -pot]
    set maxFiltersFailed       [flagVal $args -failedFiltersCutoff 9999]
    set diagTol                [flagVal $args -diagonalTolerance 0.001]
    set remVar                 [flagVal $args -remarksVariableName ""]
    set fromProtonSolventRange [flagVal $args -fromProtonSolventRange]
    set   toProtonSolventRange [flagVal $args   -toProtonSolventRange]
    
    set exceptionList [list]

    #
    # read out explicit exceptions that were already generated (by netFilter)
    #

    foreach {fromname toname} [$pot explicitInverseExceptions] {
	set curExcept [join [list $fromname $toname]]
	lappend exceptionList $curExcept
    }

    set nStartingExceptions [llength $exceptionList]

    #
    # generate exceptions corresponding to peakAssignments with <= maxFiltersFailed
    #

    foreach curPeak [$pot peaks] {
	Peak -this $curPeak
	foreach curPA [$curPeak peakAssignments] {
	    PeakAssignment -this $curPA
	    if {[$curPA numFiltersFailed] <= $maxFiltersFailed} {
		set curExcept [join [list [$curPA fromAssignmentName] [$curPA toAssignmentName]]]
		lappend exceptionList $curExcept
	    }
	    rename $curPA ""
	}
	rename $curPeak ""
    }

    #
    # generate exceptions corresponding to peaks that would appear near the diagonal
    #

    set fromShiftList [list]
    set   toShiftList [list]

    foreach curSA [$pot shiftAssignments] {

	ShiftAssignment -this $curSA

	if {[$curSA isFrom]} {
	    lappend fromShiftList [$curSA cget -this]
	} else {
	    lappend toShiftList [$curSA cget -this]
	}

	rename $curSA ""
    }

    set invsqrt2 [expr 1.0 / double(sqrt(2.0))]

    foreach curFromSA $fromShiftList {
	ShiftAssignment -this $curFromSA
	if {[$curFromSA hasFoldedProtonShift]} {
	    
	    foreach curToSA $toShiftList {
		ShiftAssignment -this $curToSA 
		if {[$curToSA hasFoldedProtonShift]} {
		    set diagDist [expr abs([$curFromSA foldedProtonShift] - [$curToSA foldedProtonShift]) * $invsqrt2]
		    if {$diagDist < $diagTol} {
			lappend exceptionList [join [list [$curFromSA name] [$curToSA name]]]
		    }
		}
		rename $curToSA ""
	    }
	}
	rename $curFromSA ""
    }

    #
    # self-interaction exceptions
    #

    foreach curFromSA $fromShiftList {
	ShiftAssignment -this $curFromSA
	set curFromSel [$curFromSA protonSelection]
	AtomSel -this $curFromSel
	foreach curToSA $toShiftList {
	    ShiftAssignment -this $curToSA 
	    
	    if {[$curFromSel intersects [$curToSA protonSelection]]} {
		lappend exceptionList [join [list [$curFromSA name] [$curToSA name]]]
	    }

	    rename $curToSA ""
	}
	rename $curFromSel ""
	rename $curFromSA ""
    }

    #
    # solvent exceptions
    #

    if {[llength $fromProtonSolventRange] > 0} {

	set minSolv [min $fromProtonSolventRange]
	set maxSolv [max $fromProtonSolventRange]

	foreach curFromSA $fromShiftList {
	    ShiftAssignment -this $curFromSA
	    if {[$curFromSA hasFoldedProtonShift]} {
		if {([$curFromSA foldedProtonShift] > $minSolv) &&
		    ([$curFromSA foldedProtonShift] < $maxSolv)} {
		    
		    #
		    # this SA is within the solvent range.  Create inverse exceptions 
		    # from it to everything.
		    #

		    foreach curToSA $toShiftList {
			ShiftAssignment -this $curToSA 
			lappend exceptionList [join [list [$curFromSA name] [$curToSA name]]]
			rename $curToSA ""
		    }
		}
	    }
	    rename $curFromSA ""
	}
    }

    if {[llength $toProtonSolventRange] > 0} {

	set minSolv [min $toProtonSolventRange]
	set maxSolv [max $toProtonSolventRange]

	foreach curToSA $toShiftList {
	    ShiftAssignment -this $curToSA
	    if {[$curToSA hasFoldedProtonShift]} {
		if {([$curToSA foldedProtonShift] > $minSolv) &&
		    ([$curToSA foldedProtonShift] < $maxSolv)} {

		    #
		    # this SA is within the solvent range.  Create inverse exceptions
		    # to it from everything
		    #

		    foreach curFromSA $fromShiftList {
			ShiftAssignment -this $curFromSA 
			lappend exceptionList [join [list [$curFromSA name] [$curToSA name]]]
			rename $curFromSA ""
		    }
		}
	    }
	    rename $curToSA ""
	}
    }
  
    
    #
    # for every exception that includes a shiftAssign that has an overlapping
    # proton selection with another shiftAssign
    # (eg., two expanded stereoassigns), copy it to make the other SA an
    # exception too
    #
    # from SAs first, then to SAs
    #

    foreach curFromSA $fromShiftList {
	ShiftAssignment -this $curFromSA

	set curSel [$curFromSA protonSelection] 
	AtomSel -this $curSel

	set overlaps([$curFromSA name]) [list]

	set equivs [shiftAssignmentsWithOverlappingSelections -shiftAssignments  $fromShiftList -target $curFromSA]
	foreach equivName $equivs {
	    if {$equivName != [$curFromSA name]} {
		lappend overlaps([$curFromSA name]) $equivName
	    }
	}
	rename $curSel ""
	rename $curFromSA ""
    }

    foreach curToSA $toShiftList {
	ShiftAssignment -this $curToSA

	set curSel [$curToSA protonSelection] 
	AtomSel -this $curSel

	set overlaps([$curToSA name]) [list]

	set equivs [shiftAssignmentsWithOverlappingSelections -shiftAssignments  $toShiftList -target $curToSA]
	foreach equivName $equivs {
	    if {$equivName != [$curToSA name]} {
		lappend overlaps([$curToSA name]) $equivName
	    }
	}
	rename $curSel ""
	rename $curToSA ""
    }
    
    foreach elem $exceptionList {
	set curFrom [lindex [split $elem] 0]
	set curTo   [lindex [split $elem] 1]

	foreach equiv $overlaps($curFrom) {
	    lappend exceptionList [join [list $equiv $curTo]]
	}
    }

    foreach elem $exceptionList {
	set curFrom [lindex [split $elem] 0]
	set curTo   [lindex [split $elem] 1]

	foreach equiv $overlaps($curTo) {
	    lappend exceptionList [join [list $curFrom $equiv]]
	}
    }

    set exceptionList [lsort -unique $exceptionList]
 
    $pot removeAllExplicitInverseExceptions

    foreach elem $exceptionList {
	set curFrom [lindex [split $elem] 0]
	set curTo   [lindex [split $elem] 1]

	$pot addExplicitInverseException $curFrom $curTo
    }

    if {$remVar != ""} {
	upvar $remVar tempRem

	set newRem [format "In createExplicitInverseExceptions, began with %d exceptions, and added %d more" $nStartingExceptions [expr [llength $exceptionList] - $nStartingExceptions]]
	appendLineToString newRem [format "These exceptions represent %f%% of the total possible inverse exceptions" [expr 100.0 * [llength $exceptionList] / double([llength $fromShiftList] * [llength $toShiftList])]]
	lappend tempRem $newRem
    }

    return ""
}



    

proc correctShiftOffset args {

    #
    # Find all the symmetric, intraresidue peak pairs.
    #
    # Determine the most common difference between the peak position
    # on the from proton dimension and the peak position 
    # of the same SA on the to dimension.
    #
    # If the offset is larger than a cutoff and the number of peaks 
    # used in determining it is larger than a cutoff, apply it to all the 
    # peak locations along the indicated dimension and throw an exception
    # (to let the calling proc know that peak locations have changed)
    #
    # If the correction is smaller than the cutoff, do nothing.
    #

    set pot                  [requiredFlagVal $args -pot]
    set shiftAssigns         [flagVal $args -shiftAssignments [$pot shiftAssignments]]
    set peaks                [flagVal $args -peaks [$pot peaks]]
    set remVar               [flagVal $args -remarksVariableName ""]
    set correctToDimension   [flagExists $args -correctToProtonDimension]
    set correctFromDimension [flagExists $args -correctFromProtonDimension]
    set isVerbose            [flagExists $args -verbose]
    set binWidth             [flagVal $args -binWidth 0.001]
    set minNumSymmetricPeaks [flagVal $args -minNumSymmetricPeaks 100]
    set minOffsetToApply     [flagVal $args -minOffsetToApply 0.015]

    #
    # for each shiftAssignment, record its opposite sense partner name
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
    # for each intraresidue PeakAssignment,
    # record Peak, PeakAssign, from, and to shiftAssign names in a list,
    # with an array index based on the names of its from and to shiftAssigns
    #

    set peakAssData [list]

    set count 0

    set bbnSeqSel [AtomSel -args "(name hn or name ha or name ha* or name hb or name hb*)"]

    foreach p $peaks {

	updateUser [format "Recording data for peak %d of %d \r" \
			[incr count] [llength $peaks]]

	Peak -this $p

	foreach pa [$p peakAssignments] {

	    PeakAssignment -this $pa
	    
	    if {([$pa isIntraresidue]) ||
		([$pa isSequential] &&
		 [$bbnSeqSel intersects [$pa fromProtonSelection]] &&
		 [$bbnSeqSel intersects [$pa   toProtonSelection]])} {
		
		if {[$pa hasUnfoldedFromProtonPeakPosition] &&
		    [$pa hasUnfoldedToProtonPeakPosition] && 
		    [$p hasIntensity]} {

		    set fromProtonPP [$pa unfoldedFromProtonPeakPosition]
		    set toProtonPP   [$pa unfoldedToProtonPeakPosition]
		    set intensity    [$p intensity]

		    lappend peakAssData [list [$p name] [$pa name] [$pa fromAssignmentName] [$pa toAssignmentName] \
					     $fromProtonPP $toProtonPP $intensity]
		
		    lappend peakAssIndex([join [list [$pa fromAssignmentName] [$pa toAssignmentName]]]) \
			[expr [llength $peakAssData] - 1]
		
		}
	    }
	    rename $pa ""
	}
	rename $p ""
    }
    
    rename $bbnSeqSel ""
    
    #
    # For each intraresidue/bbn-seq PeakAssignment, see if a symmetry match exists.
    # If so, record the difference between the peak location on the from and to dimensions
    #

    set count 0
    set numExamples 0
    set totDeltaShift 0

    set deltas [list]
    set peaksUsed [list]

    foreach curData $peakAssData {

	updateUser [format "Examining peak assignment %d of %d for symmetry \r" \
			[incr count] [llength $peakAssData]]

	set curPeakName       [lindex $curData 0]
	set curPAName         [lindex $curData 1]
	set curFromAssignName [lindex $curData 2]
	set curToAssignName   [lindex $curData 3]
	set curFromProtonPP   [lindex $curData 4]
	set curToProtonPP     [lindex $curData 5]
	set curIntensity      [lindex $curData 6]

	set partnerFromAssignName $toFromPartner($curToAssignName)
	set partnerToAssignName   $toFromPartner($curFromAssignName)


	# 
	# if there are no symmetry partners (based on the partner shiftAssign), then skip this entry
	#

	if {! [info exists peakAssIndex([join [list $partnerFromAssignName $partnerToAssignName]])]} {

	    continue
	}
	

	#
	# evaluate each symmetry partner's  delta shift, making sure that a peak can't be its own symmetry partner
	# In cases where a peak has more than one possible symmetry partner, include all partners' 
	# shift offset estimates, since one will be right, and the others will be random.
	#

	set matchIndexes $peakAssIndex([join [list $partnerFromAssignName $partnerToAssignName]])

	foreach matchIndex $matchIndexes {

	    set curMatchData       [lindex $peakAssData $matchIndex]
	    set curMatchPeakName   [lindex $curMatchData 0]
	    set curMatchFromAssignName [lindex $curMatchData 2]
	    set curMatchToAssignName   [lindex $curMatchData 3]
	    set curMatchToProtonPP [lindex $curMatchData 5]
	    set curMatchIntensity  [lindex $curMatchData 6]

	    if {$curMatchPeakName != $curPeakName} {

		set curDeltaShift [expr $curFromProtonPP - $curMatchToProtonPP]
		lappend deltas $curDeltaShift
		lappend peaksUsed $curPeakName
		lappend peaksUsed $curMatchPeakName
	    }
	}
    }

    set shiftOffset 0

    if {[llength $deltas] >= $minNumSymmetricPeaks} {

	#
	# divide up range of shift offset estimates into bins $binWidth wide.
	# Find the bin with the greatest number of examples.
	# Return the center of that bin as the estimated offset
	#
	
	set minOffset [min $deltas]
	set maxOffset [max $deltas]
	set nBins [expr int(ceil(($maxOffset - $minOffset) / double($binWidth)))]
	set actBinWidth [expr ($maxOffset - $minOffset) / double($nBins)]
	
	for {set count 0} {$count < $nBins} {incr count} {
	    set nExamples($count) 0
	}
	
	foreach elem $deltas {
	    set curBin [whichBin $minOffset $maxOffset $nBins $elem]
	    incr nExamples($curBin)
	}
	
	set maxExamples -1
	set maxBinCenter -1
	
	foreach binNum [array names nExamples] {
	    if {$nExamples($binNum) > $maxExamples} {
		set maxExamples $nExamples($binNum)
		set maxBinCenter [binCenter $minOffset $maxOffset $nBins $binNum]
	    }
	}
	
	set shiftOffset $maxBinCenter
	
	#
	# if the offset is large enough to bother with, correct every peak's position
	#
	
	if {[expr abs($shiftOffset)] >= [expr abs($minOffsetToApply)]} {

	    if {$correctToDimension && $correctFromDimension} {
	    
		foreach curPeak $peaks {
		    Peak -this $curPeak 
		    
		    $curPeak setFromProtonShift [expr [$curPeak fromProtonShift] - ($shiftOffset * 0.5)]
		    $curPeak setToProtonShift   [expr [$curPeak toProtonShift]   + ($shiftOffset * 0.5)]
		    
		    rename $curPeak ""
		}
		
	    } elseif {$correctToDimension} {
		
		foreach curPeak $peaks {
		    Peak -this $curPeak 
		    
		    $curPeak setToProtonShift [expr [$curPeak toProtonShift] + $shiftOffset]
		    
		    rename $curPeak ""
		}
		
	    } elseif {$correctFromDimension} {
		
		foreach curPeak $peaks {
		    Peak -this $curPeak 
		    
		    $curPeak setFromProtonShift [expr [$curPeak fromProtonShift] - $shiftOffset]
		    
		    rename $curPeak ""
		}
	    }
	}
    }

    #
    # summarize results
    #

    if {$remVar != ""} {

	set curRpt "correctShiftOffset:  "

	if {[llength $deltas] < $minNumSymmetricPeaks} {
	    appendLineToString curRpt [format "Skipping.  Not enough symmetric peaks found (need %d, have %d)" $minNumSymmetricPeaks [llength $deltas]]
	} else {

	    if {$isVerbose} {
		foreach elem $deltas {
		    appendLineToString curRpt [format "shift offset estimate: %f" $elem]
		}
	    }
	    
	    appendLineToString curRpt [histogram -data $deltas -title "Estimates of shift offset correction"]
	    
	    if {[expr abs($shiftOffset)] < [expr abs($minOffsetToApply)]} {
		appendLineToString curRpt [format "Shift offset of %f ppm not applied, because it's smaller than minimum offset %f ppm" \
					       $shiftOffset $minOffsetToApply]
	    } else {

		if {$correctToDimension && $correctFromDimension} {
		    appendLineToString curRpt "Correcting both from and to proton dimensions of the peak locations, by applying half the correction to each"
		} elseif {$correctToDimension} {
		    appendLineToString curRpt "Correcting to proton dimension of the peak locations"
		} elseif {$correctFromDimension} {
		    appendLineToString curRpt "Correcting from proton dimension of the peak locations"
		} else {
		    appendLineToString curRpt "Not applying any correction to the peak locations"
		}
	    
		appendLineToString curRpt [format "Shift offset estimates ranged from %f to %f ppm.  Divided that range into %d bins %f ppm wide to find most-populated" \
					       $maxOffset $minOffset $nBins $actBinWidth]
	    
		appendLineToString curRpt [format "Best correction is %f ppm, based on %d symmetric intraresidue or backbone-sequential peakAssignments on %d peaks" \
					       $shiftOffset [llength $deltas] [llength [lsort -unique $peaksUsed]]]
		
	    }
	}

	upvar $remVar tempRem
	lappend tempRem $curRpt
    }

    #
    # throw an error if I applied a shift offset, to let the caller know that
    # I changed the peak locations
    #

    if {$shiftOffset >= $minOffsetToApply} {
	error "Shift offset detected--peak locations changed"
    }

    return ""
}













proc initializeLikelihoodsFromFilters args {

    set pot              [requiredFlagVal $args -pot]
    set maxFiltersFailed [flagVal $args -maxFiltersFailed 0]
    
    foreach curPeak [$pot peaks] {
	Peak -this $curPeak
	foreach curPA [$curPeak peakAssignments] {
	    PeakAssignment -this $curPA

	    if {[$curPA numFiltersFailed] > $maxFiltersFailed} {
		$curPA setPreviousLikelihood 0
	    } else {
		$curPA setPreviousLikelihood 1
	    }
	    
	    rename $curPA ""
	} 
	rename $curPeak ""
    }

    return ""
}


proc standard2dInitMatch args {

    #
    # gathers together the standard initial matching script for 2D spectra
    #
    
    set pot                  [requiredFlagVal $args -pot]
    set fromProtonRange      [requiredFlagVal $args -fromProtonSpectralRangePPM]
    set toProtonRange        [requiredFlagVal $args -toProtonSpectralRangePPM]
    set basePhase            [requiredFlagVal $args -basePhase]
    set refStructFileName    [flagVal $args -referenceStructureFile]
    set exceptionsFileName   [flagVal $args -exceptionsFileName]
    set peakFileName         [flagVal $args -peakFileName]
    set saFileName           [flagVal $args -shiftAssignmentsFileName]
    set peakRemVar           [requiredFlagVal $args -peakRemarksVariableName]
    set saRemVar             [requiredFlagVal $args -saRemarksVariableName]
    set writeEachStage       [flagExists $args -writeEachStage]
    set allWeakPeaks         [flagExists $args -allWeakPeaks]
    set fromProtonSolventRange [flagVal $args -fromProtonSolventRange "NONE"]
    set toProtonSolventRange   [flagVal $args -toProtonSolventRange "NONE"]
    set fromProtonBroadTol   [flagVal $args -fromProtonBroadTolerancePPM 0.075]
    set toProtonBroadTol     [flagVal $args   -toProtonBroadTolerancePPM 0.075]
    set fromProtonTightTol   [flagVal $args -fromProtonTightTolerancePPM 0.02]
    set toProtonTightTol     [flagVal $args   -toProtonTightTolerancePPM 0.02]
    set doStripeCorrection   [expr ! [flagExists $args  -noStripeCorrection]]
    set doCorrectShiftOffset [expr ! [flagExists $args  -noCorrectShiftOffset]]

    marvinPyth command {import pasd}
    set defAveExp [lindex \
		    [lindex \
			 [marvinPyth command {out=pasd.aveExp} \
			      {} {out}] \
			 0] 1]
    set aveExp               [flagVal $args -aveExp $defAveExp]


    upvar $peakRemVar peakRemarks
    upvar $saRemVar   saRemarks

    global errorInfo
    


    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "fromProton" \
	-minVal [max $fromProtonRange] \
	-maxVal 99999 \
	-remarksVariableName saRemarks

    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "fromProton" \
	-minVal -99999 \
	-maxVal [min $fromProtonRange] \
	-remarksVariableName saRemarks

    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "toProton" \
	-minVal [max $toProtonRange] \
	-maxVal 99999 \
	-remarksVariableName saRemarks

    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "toProton" \
	-minVal -99999 \
	-maxVal [min $toProtonRange] \
	-remarksVariableName saRemarks    

    removeSolventPeaks \
 	-peakList [$pot peaks] \
 	-pot $pot \
 	-autoDetect \
 	-fromProtonDimension \
 	-tolerance 0.05 \
	-solventRangeVariableName fromProtonSolventRange \
 	-remarksVariableName peakRemarks

    removeSolventPeaks \
 	-peakList [$pot peaks] \
 	-pot $pot \
 	-autoDetect \
 	-toProtonDimension \
 	-tolerance 0.05 \
	-solventRangeVariableName toProtonSolventRange \
 	-remarksVariableName peakRemarks


    consolidateIdenticalStereopartners \
	-pot $pot \
	-remarksVariableName saRemarks

    recordStereoPartnersForShiftAssignments \
    	-pot $pot \
	-remarksVariableName saRemarks

    expandStereospecificShiftAssignments \
	-shiftAssignments [$pot shiftAssignments] \
	-remarksVariableName saRemarks
    
    recordToFromPartnersForShiftAssignments \
	-pot $pot \
	-remarksVariableName saRemarks
    
    markMethylShiftAssignments \
	-shiftAssignments [$pot shiftAssignments] \
	-remarksVariableName saRemarks

    if {$basePhase == "negative"} {
	set basePhaseFlag "-nonFoldedPeaksAreNegative"
    } else {
	set basePhaseFlag "-nonFoldedPeaksArePositive"
    }

    #
    # broad-tolerance match
    #
    
    match2d \
 	-peakList [$pot peaks] \
 	-fromProtonTolerancePPM $fromProtonBroadTol \
 	-fromProtonSpectralRangePPM  $fromProtonRange \
 	-toProtonTolerancePPM $toProtonBroadTol \
 	-toProtonSpectralRangePPM $toProtonRange \
 	-pot $pot \
 	-remarksVariableName peakRemarks \
 	-shiftsAliasedAlongFromProtonDimension \
 	-shiftsAliasedAlongToProtonDimension \
 	-verboseOutsideSpectrum \
 	$basePhaseFlag \
	-aveExp $aveExp
    
    #
    # correct the shift offset between the proton dimensions,
    # yielding updated peak locations
    #
    
    set needsRematch 0
    
    if { $doCorrectShiftOffset } {
	if { [catch {correctShiftOffset \
			 -pot $pot  \
			 -correctToProtonDimension \
			 -remarksVariableName peakRemarks} msg] } {
	    
	    if {$msg == "Shift offset detected--peak locations changed"} {
		set needsRematch 1
		set errorInfo ""
	    } else {
		set savedInfo $errorInfo
		error $savedInfo
	    }
	}
    }
    
    
    #
    # clean up the peak list based on global spectral features
    #
    
    removeDiagonalPeaks \
 	-peakList [$pot peaks] \
 	-pot $pot \
 	-tolerance [expr 0.25 * ($fromProtonTightTol + $fromProtonTightTol)] \
 	-remarksVariableName peakRemarks
        

    #
    # Optional broad-tolerance rematch, if a shift offset was detected
    #
    
    if {$needsRematch} {

 	match2d \
 	    -peakList [$pot peaks] \
 	    -fromProtonTolerancePPM $fromProtonBroadTol \
 	    -fromProtonSpectralRangePPM  $fromProtonRange \
 	    -toProtonTolerancePPM $toProtonBroadTol \
 	    -toProtonSpectralRangePPM $toProtonRange \
 	    -pot $pot \
 	    -remarksVariableName peakRemarks \
 	    -shiftsAliasedAlongFromProtonDimension \
 	    -shiftsAliasedAlongToProtonDimension \
 	    $basePhaseFlag \
	    -aveExp $aveExp
    }
    
    #
    # generate the distance bounds
    #
    
    if {$allWeakPeaks} {
 	foreach curPeak [$pot peaks] {
 	    Peak -this $curPeak
 	    $curPeak setUpBound 5.0
 	    $curPeak setLowBound 1.8
 	    rename $curPeak ""
 	}
    } else {
 	generateDistanceBounds -peakList [$pot peaks]
    }
    
    correctUpBoundsForMethyls -peakList [$pot peaks]

    
    #
    # generate inverse exceptions for statistics gathering
    #
    
    createExplicitInverseExceptions \
	-pot $pot \
	-remarksVariableName saRemarks \
	-failedFiltersCutoff -1 \
	-fromProtonSolventRange $fromProtonSolventRange \
	-toProtonSolventRange $toProtonSolventRange
    
    if {$writeEachStage} {
	
	if {$exceptionsFileName != ""} {
	    writeExplicitInverseExceptions -pot $pot -filename [format "%s.stage1" $exceptionsFileName]
	}
	
	if {$saFileName != ""} {
	    writeShiftAssignments \
		-fileName [format "%s.stage1" $saFileName] \
		-shiftAssignments [$pot shiftAssignments]
	}
	
	if {$peakFileName != ""} {
	    writeMarvinPeaks \
		-fileName [format "%s.stage1" $peakFileName] \
		-pot $pot 
	}
    }

    updateUser [format "\nAfter stage1 write\n"]

    #
    # record performance for posterity
    #
    
    initialShiftAssignmentAnalysis \
 	-pot $pot \
 	-referenceStructureFile $refStructFileName \
 	-minLikelihood 0.9 \
 	-description "Broad-tolerance match" \
 	-remarksVariableName saRemarks

    initialPeakAnalysis \
 	-pot $pot \
 	-violCutoff 0.5 \
 	-minLikelihood 0.9 \
 	-referenceStructureFile $refStructFileName \
 	-description "Broad-tolerance match" \
 	-remarksVariableName peakRemarks

    updateUser [format "\nAfter SA/peak analysis \n"]

    #
    # remove exceptions
    #
    
    $pot removeAllExplicitInverseExceptions
    
    #
    # correct the shift assignments' chemical shifts by detecting stripes,
    # yielding updated shift assignment values
    #
    
    updateUser [format "\nStarting newStripe\n"]

    if {$doStripeCorrection} {
	stripeCorrection \
	    -pot $pot \
	    -remarksVariableName saRemarks
    }



    #
    # Tight tolerance re-match
    #
    
    match2d \
	-peakList [$pot peaks] \
	-fromProtonTolerancePPM $fromProtonTightTol \
	-fromProtonSpectralRangePPM  $fromProtonRange \
	-toProtonTolerancePPM $toProtonTightTol \
	-toProtonSpectralRangePPM $toProtonRange \
	-pot $pot \
	-remarksVariableName peakRemarks \
	-shiftsAliasedAlongFromProtonDimension \
	-shiftsAliasedAlongToProtonDimension \
	$basePhaseFlag \
	-aveExp $aveExp
    
    #
    # generate the distance bounds
    #
    
    if {$allWeakPeaks} {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    $curPeak setUpBound 5.0
	    $curPeak setLowBound 1.8
	    rename $curPeak ""
	}
    } else {
	generateDistanceBounds -peakList [$pot peaks]
    }

    correctUpBoundsForMethyls -peakList [$pot peaks]
    
    #
    # generate exceptions for statistics gathering
    #

    createExplicitInverseExceptions \
	-pot $pot \
	-remarksVariableName saRemarks \
	-failedFiltersCutoff -1 \
	-fromProtonSolventRange $fromProtonSolventRange \
	-toProtonSolventRange $toProtonSolventRange

    if {$writeEachStage} {

	if {$exceptionsFileName != ""} {
	    writeExplicitInverseExceptions -pot $pot -filename [format "%s.stage2" $exceptionsFileName]
	}

	if {$saFileName != ""} {
	    writeShiftAssignments \
		-fileName [format "%s.stage2" $saFileName] \
		-shiftAssignments [$pot shiftAssignments]
	}
	
	if {$peakFileName != ""} {
	    writeMarvinPeaks \
		-fileName [format "%s.stage2" $peakFileName] \
		-pot $pot 
	}
    }


    #
    # record performance of tight-tolerance match for posterity
    #

    initialShiftAssignmentAnalysis \
	-pot $pot \
	-referenceStructureFile $refStructFileName \
	-minLikelihood 0.9 \
	-description "Tight-tolerance match" \
	-remarksVariableName saRemarks

    initialPeakAnalysis \
	-pot $pot \
	-violCutoff 0.5 \
	-minLikelihood 0.9 \
	-referenceStructureFile $refStructFileName \
	-description "Tight-tolerance match" \
	-remarksVariableName peakRemarks


    if {$exceptionsFileName != ""} {
	writeExplicitInverseExceptions -pot $pot -filename $exceptionsFileName
    }

    #
    # report statistics
    #
    
    initialShiftAssignmentAnalysis \
	-pot $pot \
	-referenceStructureFile $refStructFileName \
	-minLikelihood 0.9 \
	-description "end of initial match" \
	-remarksVariableName saRemarks
     
    if {$saFileName != ""} {
	writeShiftAssignments \
	    -fileName $saFileName \
	    -shiftAssignments [$pot shiftAssignments] \
	    -remarks $saRemarks
    }
        
    initialPeakAnalysis \
	-pot $pot \
	-violCutoff 0.5 \
	-minLikelihood 0.9 \
	-description "end of initial match" \
	-referenceStructureFile $refStructFileName \
	-remarksVariableName peakRemarks
     
    if {$peakFileName != ""} {
	writeMarvinPeaks \
	    -fileName $peakFileName \
	    -pot $pot \
	    -remarks $peakRemarks
    }
    

    return [list $peakRemarks $saRemarks]
}

proc standard3dInitMatch args {

    #
    # gathers together the standard initial matching script for 3D spectra
    #
    
    set pot                  [requiredFlagVal $args -pot]
    set fromProtonRange      [requiredFlagVal $args -fromProtonSpectralRangePPM]
    set fromHeavyatomRange   [requiredFlagVal $args -fromHeavyatomSpectralRangePPM]
    set toProtonRange        [requiredFlagVal $args -toProtonSpectralRangePPM]
    set basePhase            [requiredFlagVal $args -basePhase]
    set refStructFileName    [flagVal $args -referenceStructureFile]
    set exceptionsFileName   [flagVal $args -exceptionsFileName]
    set peakFileName         [flagVal $args -peakFileName]
    set saFileName           [flagVal $args -shiftAssignmentsFileName]
    set peakRemVar           [requiredFlagVal $args -peakRemarksVariableName]
    set saRemVar             [requiredFlagVal $args -saRemarksVariableName]
    set writeEachStage       [flagExists $args -writeEachStage]
    set allWeakPeaks         [flagExists $args -allWeakPeaks]
    set fromProtonSolventRange [flagVal $args -fromProtonSolventRange "NONE"]
    set toProtonSolventRange   [flagVal $args -toProtonSolventRange "NONE"]
    set fromProtonBroadTol   [flagVal $args -fromProtonBroadTolerancePPM 0.075]
    set toProtonBroadTol     [flagVal $args   -toProtonBroadTolerancePPM 0.075]
    set fromHeavyBroadTol    [flagVal $args -fromHeavyatomBroadTolerancePPM 0.75]
    set fromProtonTightTol   [flagVal $args -fromProtonTightTolerancePPM 0.02]
    set toProtonTightTol     [flagVal $args   -toProtonTightTolerancePPM 0.02]
    set fromHeavyTightTol    [flagVal $args -fromHeavyatomTightTolerancePPM 0.2]
    set doStripeCorrection   [expr ! [flagExists $args  -noStripeCorrection]]
    set doCorrectShiftOffset [expr ! [flagExists $args  -noCorrectShiftOffset]]

    marvinPyth command {import pasd}
    set defAveExp [lindex \
		    [lindex \
			 [marvinPyth command {out=pasd.aveExp} \
			      {} {out}] \
			 0] 1]
    set aveExp                [flagVal $args -aveExp $defAveExp]


    upvar $peakRemVar peakRemarks
    upvar $saRemVar   saRemarks

    global errorInfo
    


    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "fromProton" \
	-minVal [max $fromProtonRange] \
	-maxVal 99999 \
	-remarksVariableName saRemarks

    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "fromProton" \
	-minVal -99999 \
	-maxVal [min $fromProtonRange] \
	-remarksVariableName saRemarks

    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "toProton" \
	-minVal [max $toProtonRange] \
	-maxVal 99999 \
	-remarksVariableName saRemarks

    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "toProton" \
	-minVal -99999 \
	-maxVal [min $toProtonRange] \
	-remarksVariableName saRemarks    

    removeSolventPeaks \
 	-peakList [$pot peaks] \
 	-pot $pot \
 	-autoDetect \
 	-fromProtonDimension \
 	-tolerance 0.05 \
	-solventRangeVariableName fromProtonSolventRange \
 	-remarksVariableName peakRemarks

    removeSolventPeaks \
 	-peakList [$pot peaks] \
 	-pot $pot \
 	-autoDetect \
 	-toProtonDimension \
 	-tolerance 0.05 \
	-solventRangeVariableName toProtonSolventRange \
 	-remarksVariableName peakRemarks


    consolidateIdenticalStereopartners \
	-pot $pot \
	-remarksVariableName saRemarks

    recordStereoPartnersForShiftAssignments \
    	-pot $pot \
	-remarksVariableName saRemarks

    expandStereospecificShiftAssignments \
	-shiftAssignments [$pot shiftAssignments] \
	-remarksVariableName saRemarks
    
    recordToFromPartnersForShiftAssignments \
	-pot $pot \
	-remarksVariableName saRemarks
    
    markMethylShiftAssignments \
	-shiftAssignments [$pot shiftAssignments] \
	-remarksVariableName saRemarks

    if {$basePhase == "negative"} {
	set basePhaseFlag "-nonFoldedPeaksAreNegative"
    } else {
	set basePhaseFlag "-nonFoldedPeaksArePositive"
    }

    #
    # broad-tolerance match
    #
    
    match3d \
 	-peakList [$pot peaks] \
 	-fromProtonTolerancePPM $fromProtonBroadTol \
 	-fromProtonSpectralRangePPM  $fromProtonRange \
 	-fromHeavyatomTolerancePPM $fromHeavyBroadTol \
 	-fromHeavyatomSpectralRangePPM $fromHeavyatomRange \
 	-toProtonTolerancePPM $toProtonBroadTol \
 	-toProtonSpectralRangePPM $toProtonRange \
 	-pot $pot \
 	-remarksVariableName peakRemarks \
 	-signChangesUponFoldingFromHeavyatomDimension\
 	-shiftsAliasedAlongFromHeavyatomDimension \
 	-shiftsAliasedAlongFromProtonDimension \
 	-shiftsAliasedAlongToProtonDimension \
 	-verboseOutsideSpectrum \
 	$basePhaseFlag \
	-aveExp $aveExp
    
    #
    # correct the shift offset between the proton dimensions,
    # yielding updated peak locations
    #
    
    set needsRematch 0
    
    if { $doCorrectShiftOffset } {
	if { [catch {correctShiftOffset \
			 -pot $pot  \
			 -correctToProtonDimension \
			 -remarksVariableName peakRemarks} msg] } {
	
	    if {$msg == "Shift offset detected--peak locations changed"} {
		set needsRematch 1
		set errorInfo ""
	    } else {
		set savedInfo $errorInfo
		error $savedInfo
	    }
	}
    }
    
    
    #
    # clean up the peak list based on global spectral features
    #
    
    removeDiagonalPeaks \
 	-peakList [$pot peaks] \
 	-pot $pot \
 	-tolerance [expr 0.25 * ($fromProtonTightTol + $fromProtonTightTol)] \
 	-remarksVariableName peakRemarks
        

    #
    # Optional broad-tolerance rematch, if a shift offset was detected
    #
    
    if {$needsRematch} {

 	match3d \
 	    -peakList [$pot peaks] \
 	    -fromProtonTolerancePPM $fromProtonBroadTol \
 	    -fromProtonSpectralRangePPM  $fromProtonRange \
 	    -fromHeavyatomTolerancePPM $fromHeavyBroadTol \
 	    -fromHeavyatomSpectralRangePPM $fromHeavyatomRange \
 	    -toProtonTolerancePPM $toProtonBroadTol \
 	    -toProtonSpectralRangePPM $toProtonRange \
 	    -pot $pot \
 	    -remarksVariableName peakRemarks \
 	    -signChangesUponFoldingFromHeavyatomDimension\
 	    -shiftsAliasedAlongFromHeavyatomDimension \
 	    -shiftsAliasedAlongFromProtonDimension \
 	    -shiftsAliasedAlongToProtonDimension \
 	    $basePhaseFlag \
	    -aveExp $aveExp
    }
    
    #
    # generate the distance bounds
    #
    
    if {$allWeakPeaks} {
 	foreach curPeak [$pot peaks] {
 	    Peak -this $curPeak
 	    $curPeak setUpBound 5.0
 	    $curPeak setLowBound 1.8
 	    rename $curPeak ""
 	}
    } else {
 	generateDistanceBounds -peakList [$pot peaks]
    }
    
    correctUpBoundsForMethyls -peakList [$pot peaks]

    
    #
    # generate inverse exceptions for statistics gathering
    #
    
    createExplicitInverseExceptions \
	-pot $pot \
	-remarksVariableName saRemarks \
	-failedFiltersCutoff -1 \
	-fromProtonSolventRange $fromProtonSolventRange \
	-toProtonSolventRange $toProtonSolventRange
    
    if {$writeEachStage} {
	
	if {$exceptionsFileName != ""} {
	    writeExplicitInverseExceptions -pot $pot -filename [format "%s.stage1" $exceptionsFileName]
	}
	
	if {$saFileName != ""} {
	    writeShiftAssignments \
		-fileName [format "%s.stage1" $saFileName] \
		-shiftAssignments [$pot shiftAssignments]
	}
	
	if {$peakFileName != ""} {
	    writeMarvinPeaks \
		-fileName [format "%s.stage1" $peakFileName] \
		-pot $pot 
	}
    }

    updateUser [format "\nAfter stage1 write\n"]

    #
    # record performance for posterity
    #
    
    initialShiftAssignmentAnalysis \
 	-pot $pot \
 	-referenceStructureFile $refStructFileName \
 	-minLikelihood 0.9 \
 	-description "Broad-tolerance match" \
 	-remarksVariableName saRemarks

    initialPeakAnalysis \
 	-pot $pot \
 	-violCutoff 0.5 \
 	-minLikelihood 0.9 \
 	-referenceStructureFile $refStructFileName \
 	-description "Broad-tolerance match" \
 	-remarksVariableName peakRemarks

    updateUser [format "\nAfter SA/peak analysis \n"]

    #
    # remove exceptions
    #
    
    $pot removeAllExplicitInverseExceptions
    
    #
    # correct the shift assignments' chemical shifts by detecting stripes,
    # yielding updated shift assignment values
    #
    
    updateUser [format "\nStarting newStripe\n"]

    if {$doStripeCorrection} {
	stripeCorrection \
	    -pot $pot \
	    -remarksVariableName saRemarks
    }


#      stripeCorrection\
#   	-pot $pot \
#   	-stripeGenerationProtonTolerance 0.03 \
#   	-stripeGenerationHeavyatomTolerance 0.3 \
#   	-estimateMissingTargets \
#   	-remarksVariableName peakRemarks \
#   	-geminalFilter \
#   	-toFromFilter \
#   	-toFromShare \
#   	-useBackboneSequential 

    #
    # Tight tolerance re-match
    #

    match3d \
	-peakList [$pot peaks] \
	-fromProtonTolerancePPM $fromProtonTightTol \
	-fromProtonSpectralRangePPM  $fromProtonRange \
	-fromHeavyatomTolerancePPM $fromHeavyTightTol \
	-fromHeavyatomSpectralRangePPM $fromHeavyatomRange \
	-toProtonTolerancePPM $toProtonTightTol \
	-toProtonSpectralRangePPM $toProtonRange \
	-pot $pot \
	-remarksVariableName peakRemarks \
	-signChangesUponFoldingFromHeavyatomDimension\
	-shiftsAliasedAlongFromHeavyatomDimension \
	-shiftsAliasedAlongFromProtonDimension \
	-shiftsAliasedAlongToProtonDimension \
	$basePhaseFlag \
	-aveExp $aveExp
    
    #
    # generate the distance bounds
    #
    
    if {$allWeakPeaks} {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    $curPeak setUpBound 5.0
	    $curPeak setLowBound 1.8
	    rename $curPeak ""
	}
    } else {
	generateDistanceBounds -peakList [$pot peaks]
    }

    correctUpBoundsForMethyls -peakList [$pot peaks]
    
    #
    # generate exceptions for statistics gathering
    #

    createExplicitInverseExceptions \
	-pot $pot \
	-remarksVariableName saRemarks \
	-failedFiltersCutoff -1 \
	-fromProtonSolventRange $fromProtonSolventRange \
	-toProtonSolventRange $toProtonSolventRange

    if {$writeEachStage} {

	if {$exceptionsFileName != ""} {
	    writeExplicitInverseExceptions -pot $pot -filename [format "%s.stage2" $exceptionsFileName]
	}

	if {$saFileName != ""} {
	    writeShiftAssignments \
		-fileName [format "%s.stage2" $saFileName] \
		-shiftAssignments [$pot shiftAssignments]
	}
	
	if {$peakFileName != ""} {
	    writeMarvinPeaks \
		-fileName [format "%s.stage2" $peakFileName] \
		-pot $pot 
	}
    }


    #
    # record performance of tight-tolerance match for posterity
    #

    initialShiftAssignmentAnalysis \
	-pot $pot \
	-referenceStructureFile $refStructFileName \
	-minLikelihood 0.9 \
	-description "Tight-tolerance match" \
	-remarksVariableName saRemarks

    initialPeakAnalysis \
	-pot $pot \
	-violCutoff 0.5 \
	-minLikelihood 0.9 \
	-referenceStructureFile $refStructFileName \
	-description "Tight-tolerance match" \
	-remarksVariableName peakRemarks


    if {$exceptionsFileName != ""} {
	writeExplicitInverseExceptions -pot $pot -filename $exceptionsFileName
    }

    #
    # report statistics
    #
    
    initialShiftAssignmentAnalysis \
	-pot $pot \
	-referenceStructureFile $refStructFileName \
	-minLikelihood 0.9 \
	-description "end of initial match" \
	-remarksVariableName saRemarks
     
    if {$saFileName != ""} {
	writeShiftAssignments \
	    -fileName $saFileName \
	    -shiftAssignments [$pot shiftAssignments] \
	    -remarks $saRemarks
    }
        
    initialPeakAnalysis \
	-pot $pot \
	-violCutoff 0.5 \
	-minLikelihood 0.9 \
	-description "end of initial match" \
	-referenceStructureFile $refStructFileName \
	-remarksVariableName peakRemarks
     
    if {$peakFileName != ""} {
	writeMarvinPeaks \
	    -fileName $peakFileName \
	    -pot $pot \
	    -remarks $peakRemarks
    }
    

    return [list $peakRemarks $saRemarks]
}




proc standard4dInitMatch args {

    #
    # gathers together the standard initial matching script for 4D spectra
    #
    
    set pot                  [requiredFlagVal $args -pot]
    set fromProtonRange      [requiredFlagVal $args -fromProtonSpectralRangePPM]
    set fromHeavyatomRange   [requiredFlagVal $args -fromHeavyatomSpectralRangePPM]
    set toProtonRange        [requiredFlagVal $args -toProtonSpectralRangePPM]
    set toHeavyatomRange     [requiredFlagVal $args -toHeavyatomSpectralRangePPM]
    set basePhase            [requiredFlagVal $args -basePhase]
    set refStructFileName    [flagVal $args -referenceStructureFile]
    set exceptionsFileName   [requiredFlagVal $args -exceptionsFileName]
    set peakFileName         [flagVal $args -peakFileName ""]
    set saFileName           [flagVal $args -shiftAssignmentsFileName ""]
    set peakRemVar           [requiredFlagVal $args -peakRemarksVariableName]
    set saRemVar             [requiredFlagVal $args -saRemarksVariableName]
    set writeEachStage       [flagExists $args -writeEachStage]
    set allWeakPeaks         [flagExists $args -allWeakPeaks]
    set fromProtonSolventRange [flagVal $args -fromProtonSolventRange "NONE"]
    set toProtonSolventRange   [flagVal $args -toProtonSolventRange "NONE"]
    set fromProtonBroadTol   [flagVal $args -fromProtonBroadTolerancePPM 0.075]
    set toProtonBroadTol     [flagVal $args   -toProtonBroadTolerancePPM 0.075]
    set fromHeavyBroadTol    [flagVal $args -fromHeavyatomBroadTolerancePPM 0.75]
    set toHeavyBroadTol      [flagVal $args   -toHeavyatomBroadTolerancePPM 0.75]
    set fromProtonTightTol   [flagVal $args -fromProtonTightTolerancePPM 0.02]
    set toProtonTightTol     [flagVal $args   -toProtonTightTolerancePPM 0.02]
    set fromHeavyTightTol    [flagVal $args -fromHeavyatomTightTolerancePPM 0.2]
    set toHeavyTightTol      [flagVal $args   -toHeavyatomTightTolerancePPM 0.2]
    set aveExp               [flagVal $args -aveExp 6]
    set doStripeCorrection   [expr ! [flagExists $args  -noStripeCorrection]]
    set doCorrectShiftOffset [expr ! [flagExists $args  -noCorrectShiftOffset]]


    marvinPyth command {import pasd}
    set defAveExp [lindex \
		    [lindex \
			 [marvinPyth command {out=pasd.aveExp} \
			      {} {out}] \
			 0] 1]
    set aveExp               [flagVal $args -aveExp $defAveExp]


    upvar $peakRemVar peakRemarks
    upvar $saRemVar   saRemarks
    
    global errorInfo


    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "fromProton" \
	-minVal [max $fromProtonRange] \
	-maxVal 99999 \
	-remarksVariableName saRemarks

    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "fromProton" \
	-minVal -99999 \
	-maxVal [min $fromProtonRange] \
	-remarksVariableName saRemarks

    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "toProton" \
	-minVal [max $toProtonRange] \
	-maxVal 99999 \
	-remarksVariableName saRemarks

    dropShiftAssignsInShiftRange \
	-pot $pot \
	-dimension "toProton" \
	-minVal -99999 \
	-maxVal [min $toProtonRange] \
	-remarksVariableName saRemarks

    removeSolventPeaks \
 	-peakList [$pot peaks] \
 	-pot $pot \
 	-autoDetect \
 	-toProtonDimension \
	-solventRangeVariableName toProtonSolventRange \
 	-tolerance 0.05 \
 	-remarksVariableName peakRemarks
    
    removeSolventPeaks \
 	-peakList [$pot peaks] \
 	-pot $pot \
 	-autoDetect \
 	-fromProtonDimension \
	-solventRangeVariableName fromProtonSolventRange \
 	-tolerance 0.05 \
 	-remarksVariableName peakRemarks

    consolidateIdenticalStereopartners \
	-pot $pot \
	-remarksVariableName saRemarks

    recordStereoPartnersForShiftAssignments \
    	-pot $pot \
	-remarksVariableName saRemarks

    expandStereospecificShiftAssignments \
	-shiftAssignments [$pot shiftAssignments] \
	-remarksVariableName saRemarks
    
    recordToFromPartnersForShiftAssignments \
	-pot $pot \
	-remarksVariableName saRemarks
    
    markMethylShiftAssignments \
	-shiftAssignments [$pot shiftAssignments] \
	-remarksVariableName saRemarks


    if {$basePhase == "negative"} {
	set basePhaseFlag "-nonFoldedPeaksAreNegative"
    } else {
	set basePhaseFlag "-nonFoldedPeaksArePositive"
    }
    
    #
    # broad-tolerance match
    #
    
    match4d \
 	-peakList [$pot peaks] \
 	-fromProtonTolerancePPM $fromProtonBroadTol \
 	-fromProtonSpectralRangePPM  $fromProtonRange \
 	-fromHeavyatomTolerancePPM $fromHeavyBroadTol \
 	-fromHeavyatomSpectralRangePPM $fromHeavyatomRange \
 	-toProtonTolerancePPM $toProtonBroadTol \
 	-toProtonSpectralRangePPM $toProtonRange \
 	-toHeavyatomTolerancePPM $toHeavyBroadTol \
 	-toHeavyatomSpectralRangePPM $toHeavyatomRange \
 	-pot $pot \
 	-remarksVariableName peakRemarks \
 	-signChangesUponFoldingFromHeavyatomDimension\
 	-signChangesUponFoldingToHeavyatomDimension\
 	-shiftsAliasedAlongFromHeavyatomDimension \
 	-shiftsAliasedAlongToHeavyatomDimension \
 	-shiftsAliasedAlongFromProtonDimension \
 	-shiftsAliasedAlongToProtonDimension \
 	-verboseOutsideSpectrum \
 	$basePhaseFlag \
	-aveExp $aveExp
    

    #
    # correct the shift offset between the proton dimensions,
    # yielding updated peak locations
    #
    
    set needsRematch 0
    
    if { $doCorrectShiftOffset } {
	if { [catch {correctShiftOffset \
			 -pot $pot  \
			 -correctToProtonDimension \
			 -remarksVariableName peakRemarks} msg] } {
	
	    if {$msg == "Shift offset detected--peak locations changed"} {
		set needsRematch 1
		set errorInfo ""
	    } else {
		set savedInfo $errorInfo
		error $savedInfo
	    }
	}
    }
    
    #
    # clean up the peak list based on global spectral features
    #
    
    removeDiagonalPeaks \
 	-peakList [$pot peaks] \
 	-pot $pot \
 	-tolerance [expr 0.25 * ($fromProtonTightTol + $fromProtonTightTol)] \
 	-remarksVariableName peakRemarks
        
    #
    # Optional broad-tolerance rematch, if a shift offset was detected
    #
    
    if {$needsRematch} {
	
 	match4d \
 	    -peakList [$pot peaks] \
 	    -fromProtonTolerancePPM 0.075 \
 	    -fromProtonSpectralRangePPM  $fromProtonRange \
 	    -fromHeavyatomTolerancePPM 0.75 \
 	    -fromHeavyatomSpectralRangePPM $fromHeavyatomRange \
 	    -toProtonTolerancePPM 0.075 \
 	    -toProtonSpectralRangePPM $toProtonRange \
 	    -toHeavyatomTolerancePPM 0.75 \
 	    -toHeavyatomSpectralRangePPM $toHeavyatomRange \
 	    -pot $pot \
 	    -remarksVariableName peakRemarks \
	    -signChangesUponFoldingFromHeavyatomDimension\
 	    -signChangesUponFoldingToHeavyatomDimension\
 	    -shiftsAliasedAlongFromHeavyatomDimension \
 	    -shiftsAliasedAlongToHeavyatomDimension \
 	    -shiftsAliasedAlongFromProtonDimension \
 	    -shiftsAliasedAlongToProtonDimension \
 	    $basePhaseFlag \
	    -aveExp $aveExp
	
    }
    
    #
    # generate the distance bounds
    #
    
    if {$allWeakPeaks} {
 	foreach curPeak [$pot peaks] {
 	    Peak -this $curPeak
 	    $curPeak setUpBound 5.0
 	    $curPeak setLowBound 1.8
 	    rename $curPeak ""
 	}
    } else {
 	generateDistanceBounds -peakList [$pot peaks]
    }
    
    correctUpBoundsForMethyls -peakList [$pot peaks]
    
    #
    # generate inverse exceptions for statistics gathering
    #
    
    createExplicitInverseExceptions \
	-pot $pot \
	-remarksVariableName saRemarks \
	-failedFiltersCutoff -1 \
	-fromProtonSolventRange $fromProtonSolventRange \
	-toProtonSolventRange $toProtonSolventRange
    
    
    if {$writeEachStage} {
	
 	writeExplicitInverseExceptions -pot $pot -filename [format "%s.stage1" $exceptionsFileName]
	
	if {$saFileName != ""} {
	    writeShiftAssignments \
		-fileName [format "%s.stage1" $saFileName] \
		-shiftAssignments [$pot shiftAssignments]
	}
	
	if {$peakFileName != ""} {
	    writeMarvinPeaks \
		-fileName [format "%s.stage1" $peakFileName] \
		-pot $pot
	}
    }
    
    #
    # record performance for posterity
    #
    
    initialShiftAssignmentAnalysis \
 	-pot $pot \
 	-referenceStructureFile $refStructFileName \
 	-minLikelihood 0.9 \
 	-description "Broad-tolerance match" \
 	-remarksVariableName saRemarks
    
    initialPeakAnalysis \
 	-pot $pot \
 	-violCutoff 0.5 \
 	-minLikelihood 0.9 \
 	-referenceStructureFile $refStructFileName \
 	-description "Broad-tolerance match" \
 	-remarksVariableName peakRemarks
    
    #
    # remove exceptions
    #
    
    $pot removeAllExplicitInverseExceptions
    
    #
    # correct the shift assignments' chemical shifts by detecting stripes,
    # yielding updated shift assignment values
    #
    
    if {$doStripeCorrection} {
	stripeCorrection \
	    -pot $pot \
	    -remarksVariableName saRemarks
    }
    
    
#     stripeCorrection\
#   	-pot $pot \
#   	-stripeGenerationProtonTolerance 0.03 \
#   	-stripeGenerationHeavyatomTolerance 0.3 \
#   	-estimateMissingTargets \
#   	-remarksVariableName peakRemarks \
#   	-geminalFilter \
#   	-toFromFilter \
#   	-toFromShare \
#   	-useBackboneSequential 
    
    #
    # Third matching stage, done with tight tolerances
    #

    match4d \
 	-peakList [$pot peaks] \
 	-fromProtonTolerancePPM $fromProtonTightTol \
 	-fromProtonSpectralRangePPM  $fromProtonRange \
 	-fromHeavyatomTolerancePPM $fromHeavyTightTol \
 	-fromHeavyatomSpectralRangePPM $fromHeavyatomRange \
 	-toProtonTolerancePPM $toProtonTightTol \
 	-toProtonSpectralRangePPM $toProtonRange \
 	-toHeavyatomTolerancePPM $toHeavyTightTol \
 	-toHeavyatomSpectralRangePPM $toHeavyatomRange \
 	-pot $pot \
 	-remarksVariableName peakRemarks \
 	-signChangesUponFoldingFromHeavyatomDimension\
     	-signChangesUponFoldingToHeavyatomDimension\
 	-shiftsAliasedAlongFromHeavyatomDimension \
 	-shiftsAliasedAlongToHeavyatomDimension \
 	-shiftsAliasedAlongFromProtonDimension \
 	-shiftsAliasedAlongToProtonDimension \
 	$basePhaseFlag \
	-aveExp $aveExp


    #
    # generate the distance bounds
    #
    
    if {$allWeakPeaks} {
	foreach curPeak [$pot peaks] {
	    Peak -this $curPeak
	    $curPeak setUpBound 5.0
	    $curPeak setLowBound 1.8
	    rename $curPeak ""
	}
    } else {
	generateDistanceBounds -peakList [$pot peaks]
    }

    correctUpBoundsForMethyls -peakList [$pot peaks]
    
    #
    # generate exceptions for statistics gathering
    #

    createExplicitInverseExceptions \
	-pot $pot \
	-remarksVariableName saRemarks \
	-failedFiltersCutoff -1 \
	-fromProtonSolventRange $fromProtonSolventRange \
	-toProtonSolventRange $toProtonSolventRange


    if {$writeEachStage} {

	writeExplicitInverseExceptions -pot $pot -filename [format "%s.stage2" $exceptionsFileName]

	if {$saFileName != ""} {
	    writeShiftAssignments \
		-fileName [format "%s.stage2" $saFileName] \
		-shiftAssignments [$pot shiftAssignments]
	}
	
	if {$peakFileName != ""} {
	    writeMarvinPeaks \
		-fileName [format "%s.stage2" $peakFileName] \
		-pot $pot
	}
    }

    #
    # record performance of tight-tolerance match for posterity
    #

    initialShiftAssignmentAnalysis \
	-pot $pot \
	-referenceStructureFile $refStructFileName \
	-minLikelihood 0.9 \
	-description "Tight-tolerance match" \
	-remarksVariableName saRemarks

    initialPeakAnalysis \
	-pot $pot \
	-violCutoff 0.5 \
	-minLikelihood 0.9 \
	-referenceStructureFile $refStructFileName \
	-description "Tight-tolerance match" \
	-remarksVariableName peakRemarks


    #
    # report statistics.  If there is a known structure, you can 
    # enter it here to generate accuracy statistics in the starting
    # NOE file.
    #
    
    writeExplicitInverseExceptions -pot $pot -filename $exceptionsFileName

    initialShiftAssignmentAnalysis \
	-pot $pot \
	-referenceStructureFile $refStructFileName \
	-minLikelihood 0.9 \
	-description "end of initial match" \
	-remarksVariableName saRemarks
        
    if {$saFileName != ""} {
	writeShiftAssignments \
	    -fileName $saFileName \
	    -shiftAssignments [$pot shiftAssignments] \
	    -remarks $saRemarks
    }
        
    initialPeakAnalysis \
	-pot $pot \
	-violCutoff 0.5 \
	-minLikelihood 0.9 \
	-description "end of initial match" \
	-referenceStructureFile $refStructFileName \
	-remarksVariableName peakRemarks
        
    if {$peakFileName != ""} {
	writeMarvinPeaks \
	    -fileName $peakFileName \
	    -pot $pot \
	    -remarks $peakRemarks
    }
    

    return [list $peakRemarks $saRemarks]
}



proc standardJointFilter args {

    #
    # gathers together the standard joint-network-filtering script
    #
    
    set potList              [requiredFlagVal $args -potList]
    set potListForMap        [flagVal $args -potListForMap $potList]
    set refStructFileName    [flagVal $args -referenceStructureFile]
    set exceptionsFileNames  [requiredFlagVal $args -exceptionsFileNameList]
    set peakFileNames        [requiredFlagVal $args -peakFileNameList]
    set saFileNames          [requiredFlagVal $args -shiftAssignmentsFileNameList]    
    set minLikelihood        [flagVal $args -minLikelihood 0.9]
    set maxLikelihood        [flagVal $args -maxLikelihood 2.0]
    set passNetFrac          [flagVal $args -passNetFrac 0.1]
    set minExpectedNetScore  [flagVal $args -minExpectedNetScore 0.2]
    set writeEachStage       [flagExists $args -writeEachStage]
    set primarySeqFilter     [flagExists $args -primarySequenceDistanceFilter]
    set expectedContacts     [flagVal $args -expectedContactsPerResidue 8]
    set dontPreferIntra      [flagVal $args -dontPreferIntra 0]

    foreach curPot $potList {
	initializePeakAssignmentNumFiltersFailed -pot $curPot
    }

    #convert TCL list of PASDPots to Python sequence
    marvinPyth command "from pasdPot import PASDPot"
    marvinPyth command "from pyInterp import fromStringRep"

    marvinPyth command {pots=[]}
    foreach pot $potList {
	set ptr [$pot smartPtr]
	$pot -disown
	$ptr -disown
	marvinPyth command "pots.append(ptr)" [list [list ptr $ptr]]
    }

    marvinPyth command {pots = [fromStringRep(ptr) for ptr in pots]} 

    marvinPyth command {potsToFilter=[]}
    foreach pot $potList {
	set ptr [$pot smartPtr]
	$ptr -disown
	marvinPyth command "potsToFilter.append(ptr)" [list [list ptr $ptr]]
    }
    marvinPyth command {potsToFilter = [fromStringRep(ptr)
					for ptr in potsToFilter]} 
    marvinPyth command {for pot in potsToFilter: pot.thisown=False}


    marvinPyth command "from pasd.netfilter import netFilter"
    set pythCmd "remarks=netFilter(pots,potsToFilter,
                                   passFrac=$passNetFrac,
				   minPeakScore=-1,
                                   minExpectedScore=$minExpectedNetScore,
                                   printResiduePairScores=True,
				   verbose=2,
				   tclOutput=True,
                                   preferIntra=(not $dontPreferIntra))"
    set peakRemarksJoint [lindex \
			      [lindex \
				   [marvinPyth command $pythCmd {} {remarks}] \
					 0] 1]

    foreach curPot $potList {
	foreach curPeak [$curPot peaks] {
	    Peak -this $curPeak
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		rename $curPA ""
	    }
	    rename $curPeak ""
	}
    }

#	-dontPreferIntra \


#     set expContacts [list]
#     foreach curPot $potList {
# 	lappend expContacts 10
#     }


    set retVal [list]

    foreach curPot $potList curExceptionFileName $exceptionsFileNames curSAFileName $saFileNames curPeakFileName $peakFileNames {

	set curPeakRemarks $peakRemarksJoint
	set curSARemarks   [list]
	$curPot -disown
	
	
# 	initializeLikelihoodsFromFilters -pot $curPot
	
 	if {$writeEachStage} {
	    
 	    createExplicitInverseExceptions -pot $curPot -remarksVariableName curSARemarks -failedFiltersCutoff 0
 	    writeExplicitInverseExceptions -pot $curPot -filename [format "%s.stage3" $curExceptionFileName]
	    
 	    writeShiftAssignments \
 		-fileName [format "%s.stage3" $curSAFileName] \
 		-shiftAssignments [$curPot shiftAssignments]
	    
 	    writeMarvinPeaks \
 		-fileName [format "%s.stage3" $curPeakFileName] \
 		-pot $curPot 
 	}
	
	
 	initialShiftAssignmentAnalysis \
 	    -pot $curPot \
 	    -referenceStructureFile $refStructFileName \
 	    -minLikelihood 0.9 \
 	    -description "after joint network filter" \
 	    -remarksVariableName curSARemarks
	
 	initialPeakAnalysis \
 	    -pot $curPot \
 	    -violCutoff 0.5 \
 	    -minLikelihood $minLikelihood \
 	    -maxLikelihood $maxLikelihood \
 	    -referenceStructureFile $refStructFileName \
 	    -description "after joint network filter" \
 	    -remarksVariableName curPeakRemarks

	addPrimarySequenceExceptions \
	    -pot $curPot \
	    -remarksVariableName curPeakRemarks

	if {$primarySeqFilter} {
	    primarySequenceDistanceFilter \
		-pot $curPot \
		-intraresidue \
		-remarksVariableName curPeakRemarks    

	    initializeLikelihoodsFromFilters -pot $curPot
	}


	if {$writeEachStage} {

	    createExplicitInverseExceptions -pot $curPot -remarksVariableName curSARemarks -failedFiltersCutoff 0
	    writeExplicitInverseExceptions -pot $curPot -filename [format "%s.stage4" $curExceptionFileName]

	    writeShiftAssignments \
		-fileName [format "%s.stage4" $curSAFileName] \
		-shiftAssignments [$curPot shiftAssignments] 
	    
	    writeMarvinPeaks \
		-fileName [format "%s.stage4" $curPeakFileName] \
		-pot $curPot 
	}

	writePASDFiles -pot $curPot \
	    -referenceStructureFile $refStructFileName \
	    -exceptionFileName $curExceptionFileName \
            -peakFileName $curPeakFileName \
	    -shiftAssignmentsFileName $curSAFileName \
	    -minLikelihood $minLikelihood \
	    -maxLikelihood $maxLikelihood \
	    -peakRemarks $curPeakRemarks \
	    -shiftAssignmentRemarks $curSARemarks

	lappend retVal $curPeakRemarks
	lappend retVal $curSARemarks
    }

    return $retVal
}
	
proc writePASDFiles args {
    #
    # write out peak, shift assignment and exceptions files for the given term
    #

    set pot                 [requiredFlagVal $args -pot]
    set refPDB              [flagVal $args -referenceStructureFile]
    set exceptionsFilename  [requiredFlagVal $args -exceptionFileName]
    set peakFilename        [requiredFlagVal $args -peakFileName]
    set saFilename          [requiredFlagVal $args -shiftAssignmentsFileName]
    set minLikelihood       [flagVal $args -minLikelihood 0.9]
    set maxLikelihood       [flagVal $args -maxLikelihood 2.0]
    set peakRemarks         [flagVal $args -peakRemarks]
    set saRemarks           [flagVal $args -shiftAssignmentRemarks]

	      

    
    createExplicitInverseExceptions -pot $pot \
	-remarksVariableName saRemarks -failedFiltersCutoff 0
    writeExplicitInverseExceptions -pot $pot \
	-filename $exceptionsFilename

    initialShiftAssignmentAnalysis \
	-pot $pot \
	-minLikelihood 0.9 \
	-description "end of joint filter" \
	-referenceStructureFile $refPDB \
	-remarksVariableName saRemarks

    writeShiftAssignments \
	-fileName $saFilename \
	-shiftAssignments [$pot shiftAssignments] \
	-remarks $saRemarks

    initialPeakAnalysis \
	-pot $pot \
	-violCutoff 0.5 \
	-minLikelihood $minLikelihood \
	-maxLikelihood $maxLikelihood \
	-description "end of joint filter" \
	-referenceStructureFile $refPDB \
	-remarksVariableName peakRemarks
		
    writeMarvinPeaks \
	-fileName $peakFilename \
	-pot $pot \
	-remarks $peakRemarks
    return
}
