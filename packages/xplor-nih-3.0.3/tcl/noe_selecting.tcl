#
# procedures for selecting NOEs
#
# JJK 2/24/03
#

package provide marvin 1.0
package require xplorpot


proc selectNOEpeaks args {

    set pot       [flagVal $args -pot ""]
    set startList [flagVal $args -from ""]
    set aProc     [requiredFlagVal $args -proc]

    if {($pot == "") && ($startList == "")} {
	return [list]
    }

    if {$pot != ""} {
	set startList [$pot peaks]
    }

    set result [list]

    foreach peakPtr $startList {

	if {[eval $aProc $peakPtr]} {
	    lappend result $peakPtr
	}
    }

    return $result
}

proc markNOEpeaks args {

    set pot           [flagVal $args -pot ""]
    set startList     [flagVal $args -list ""]
    set commentString [flagVal $args -comment ""]
    set aProc         [flagVal $args -proc]

    if {($pot == "") && ($startList == "")} {
	set startList [list]
    }

    if {$pot != ""} {
	set startList [$pot peaks]
    }

    foreach curPeak $startList {

	#
	# generate a new comment, either from
	# the static string or dynamically from
	# a proc
	#

	if {$commentString != ""} {
	    set newComment $commentString
	} else {
	    set newComment [eval $aProc $curPeak]
	}

	#
	# append it to that peak's comments
	#

	Peak -this $curPeak
	$curPeak appendToNote $newComment
	rename $curPeak ""
	
    }
}


#
# selections based on primary sequence distance
#


#
# the primary sequence distance of a Peak is the minimum
# of that peak's peakAssignments' primary sequence distances
#

proc peakSequenceDist {tempPeak} {

    Peak -this $tempPeak

    set retVal [$tempPeak primarySeqDist]
    rename $tempPeak ""
    return $retVal
}

#
# The sequence distance for a peakAssignment is the
# maximum primary sequence distance for any pair of 
# atoms in the peakAssignment's selections.
#

proc peakAssignmentSequenceDist {tempPA} {

    PeakAssignment -this $tempPA

    set retVal [$tempPA primarySeqDist]
    rename $tempPA ""
    return $retVal
}

proc sequenceDistLessThan {target aPeakPtr} {

    set msd [peakSequenceDist $aPeakPtr]

    if {($msd != -1) && ($msd < $target)} {
	return 1
    } else  {
	return 0
    }
}

proc sequenceDistGreaterThan {target aPeakPtr} {

    set msd [minSequenceDist $aPeakPtr]
    if {($msd != -1) && ($msd > $target)} {
	return 1
    } else  {
	return 0
    }
}

proc isIntraresiduePeak {aPeakPtr} {

    set msd [peakSequenceDist $aPeakPtr]

    if {$msd == 0} {
	return 1
    } else {
	return 0
    }   
}

proc isSequentialPeak {aPeakPtr} {

    set msd [peakSequenceDist $aPeakPtr]

    if {$msd == 1} {
	return 1
    } else {
	return 0
    }
}

proc isShortRangePeak args {

    if {[llength $args] == 1} {
	set longRangeCutoff 6
	set aPeakPtr [lindex $args 0]
    } else {
	set longRangeCutoff [lindex $args 0]
	set aPeakPtr [lindex $args 1]
    }

    set msd [peakSequenceDist $aPeakPtr]

    if {($msd > 1) && ($msd < $longRangeCutoff)} {
	return 1
    } else {
	return 0
    }
}

proc isLongRangePeak args {

    if {[llength $args] == 1} {
	set aPeakPtr [lindex $args 0]
	set longRangeCutoff 6
    } else {
	set longRangeCutoff [lindex $args 0]
	set aPeakPtr [lindex $args 1]
    }

    set msd [peakSequenceDist $aPeakPtr]

    if {($msd != -1) && ($msd >= $longRangeCutoff)} {
	return 1
    } else {
	return 0
    }
}



#
# selections based on violations
#

proc bestViol {tempPeak} {

    Peak -this $tempPeak

    set bestViol [$tempPeak lowestViolation]
 
    rename $tempPeak ""

    return $bestViol
}


proc violLessThan {target aPeakPtr} {

    set v [bestViol $aPeakPtr]

    if {($v != -1) && ($v < $target)} {
	return 1
    } else {
	return 0
    }
}

proc violGreaterThan {target aPeakPtr} {

    if {[bestViol $aPeakPtr] > $target} {
	return 1
    } else {
	return 0
    }
}

#
# selections based on likelihood
#

proc likelihoodLessThan {target tempPeak} {

    Peak -this $tempPeak
    set l [$tempPeak previousLikelihood]
    rename $tempPeak ""

    if {($l != -1) && ($l < $target)} {
	return 1
    } else {
	return 0
    }
}

proc likelihoodLessThanOrEqualTo {target tempPeak} {

    Peak -this $tempPeak
    set l [$tempPeak previousLikelihood]
    rename $tempPeak ""

    if {($l != -1) && ($l <= $target)} {
	return 1
    } else {
	return 0
    }
}

proc likelihoodGreaterThan {target tempPeak} {

    Peak -this $tempPeak
    set l [$tempPeak previousLikelihood]
    rename $tempPeak ""

    if {$l > $target} {
	return 1
    } else {
	return 0
    }
}

proc likelihoodGreaterThanOrEqualTo {target tempPeak} {

    Peak -this $tempPeak
    set l [$tempPeak previousLikelihood]
    rename $tempPeak ""

    if {$l >= $target} {
	return 1
    } else {
	return 0
    }
}

#
# selections based on degeneracy
#

proc noeDegeneracy {tempPeak} {

    Peak -this $tempPeak

    set na [$tempPeak numPeakAssignments]
    rename $tempPeak ""

    return $na
}

proc degeneracyGreaterThan {maxDegen aPeakPtr} {

    if {[noeDegeneracy $aPeakPtr] > $maxDegen} {
	return 1
    } else {
	return 0
    }
}

proc degeneracyLessThan {maxDegen aPeakPtr} {

    if {[noeDegeneracy $aPeakPtr] < $maxDegen} {
	return 1
    } else {
	return 0
    }
}

proc isUnassigned {aPeakPtr} {

    if {[noeDegeneracy $aPeakPtr] == 0} {
	return 1
    } else {
	return 0
    }
}

#
# selections based on string analysis
#

proc peakIsNamed {target tempPeak} {

    Peak -this $tempPeak

    set na [$tempPeak name]
    rename $tempPeak ""

    if {$na == $target} {
	return 1
    } else {
	return 0
    }
}

proc peakNameIncludesString {target tempPeak} {

    Peak -this $tempPeak

    set na [$tempPeak name]
    rename $tempPeak ""

    if {[string first $target $na] != -1} {
	return 1
    } else {
	return 0
    }
}

proc peakNameDoesntIncludeString {target tempPeak} {

    Peak -this $tempPeak

    set na [$tempPeak name]
    rename $tempPeak ""

    if {[string first $target $na] == -1} {
	return 1
    } else {
	return 0
    }
}

proc peakNoteIncludesString {target tempPeak} {

    Peak -this $tempPeak

    set na [$tempPeak note]
    rename $tempPeak ""

    if {[string first $target $na] != -1} {
	return 1
    } else {
	return 0
    }
}

proc peakNoteDoesntIncludeString {target tempPeak} {

    Peak -this $tempPeak

    set na [$tempPeak note]
    rename $tempPeak ""

    if {[string first $target $na] == -1} {
	return 1
    } else {
	return 0
    }
}


