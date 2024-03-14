package provide primseqexcept 0.1
package require marvin

namespace eval PrimarySequenceExceptions {

namespace export addPrimarySequenceExceptions 

proc recordSAsInResidues {pot resList fSAsInResArrayName tSAsInResArrayName} {

    upvar $fSAsInResArrayName fromSAsInResidue
    upvar $tSAsInResArrayName   toSAsInResidue

    array unset fromSAsInResidue
    array unset   toSAsInResidue

    foreach curRes $resList {
	set fromSAsInResidue($curRes) [list]
	set   toSAsInResidue($curRes) [list]
    }
    
    foreach curSA [$pot shiftAssignments] {
	ShiftAssignment -this $curSA 
	
	set curSel [$curSA protonSelection]
	AtomSel -this $curSel
	set residuesInCurSel [residuesInSelection $curSel]
	
	foreach curResidue $residuesInCurSel {
	    if {[$curSA isFrom]} {
		lappend fromSAsInResidue($curResidue) [$curSA name]
	    } else {
		lappend toSAsInResidue($curResidue)   [$curSA name]
	    }
	}
	
	rename $curSel ""
	rename $curSA ""
    }
    
    return ""
}


proc recordSAsInNeighborhoods {PSDdegree resList fSAsInResArrayName tSAsInResArrayName fSAsInNeighArrayName tSAsInNeighArrayName} {

    
    upvar $fSAsInResArrayName fromSAsInResidue
    upvar $tSAsInResArrayName   toSAsInResidue
    
    upvar $fSAsInNeighArrayName fromSAsInNeighborhood
    upvar $tSAsInNeighArrayName   toSAsInNeighborhood
    
    foreach curRes $resList {
	
	set fromSAsInNeighborhood($curRes) $fromSAsInResidue($curRes)
	set   toSAsInNeighborhood($curRes)   $toSAsInResidue($curRes)
	
	set targRes $curRes
	for {set count 1} {$count <= $PSDdegree} {incr count} {

	    if {$targRes == ""} {
		break
	    }

	    set targRes [nextResidue $targRes]
	    
	    if {[info exists fromSAsInResidue($targRes)]} {
		set fromSAsInNeighborhood($curRes) [concat $fromSAsInNeighborhood($curRes) $fromSAsInResidue($targRes)]
	    }
	    
	    if {[info exists toSAsInResidue($targRes)]} {
		set   toSAsInNeighborhood($curRes) [concat   $toSAsInNeighborhood($curRes)   $toSAsInResidue($targRes)]
	    }
	}
	
	set targRes $curRes
	for {set count 1} {$count <= $PSDdegree} {incr count} {
	    
	    if {$targRes == ""} {
		break
	    }

	    set targRes [previousResidue $targRes] 
	    
	    if {[info exists fromSAsInResidue($targRes)]} {
		set fromSAsInNeighborhood($curRes) [concat $fromSAsInNeighborhood($curRes) $fromSAsInResidue($targRes)]
	    }
	    
	    if {[info exists toSAsInResidue($targRes)]} {
		set   toSAsInNeighborhood($curRes) [concat   $toSAsInNeighborhood($curRes)   $toSAsInResidue($targRes)]
	    }
	}
	
	set fromSAsInNeighborhood($curRes) [lsort -unique $fromSAsInNeighborhood($curRes)]
	set   toSAsInNeighborhood($curRes) [lsort -unique   $toSAsInNeighborhood($curRes)]
    }
    
    return ""
}


proc createPSDExceptions {pot resList fSAsInResArrayName tSAsInResArrayName fSAsInNeighArrayName tSAsInNeighArrayName} {

    upvar $fSAsInResArrayName fromSAsInResidue
    upvar $tSAsInResArrayName   toSAsInResidue

    upvar $fSAsInNeighArrayName fromSAsInNeighborhood
    upvar $tSAsInNeighArrayName   toSAsInNeighborhood


    #
    # count how many explicit exceptions we already have
    #

    set nStartingExceptions 0
    foreach {fromname toname} [$pot explicitInverseExceptions] {
	incr nStartingExceptions
    }

    foreach curRes $resList {

	foreach curFromSA $fromSAsInResidue($curRes) {
	    foreach otherToSA $toSAsInNeighborhood($curRes) {
		$pot addExplicitInverseException $curFromSA $otherToSA
	    }
	}

	foreach curToSA $toSAsInResidue($curRes) {
	    foreach otherFromSA $fromSAsInNeighborhood($curRes) {
		$pot addExplicitInverseException $otherFromSA $curToSA
	    }
	}
    }

    #
    # count how many explicit exceptions we have now
    #

    set nFinishedExceptions 0
    foreach {fromname toname} [$pot explicitInverseExceptions] {
	incr nFinishedExceptions
    }

    return [expr $nFinishedExceptions - $nStartingExceptions]
}


proc addPrimarySequenceExceptions args {

    set pot        [requiredFlagVal $args -pot]
    set PSDrange   [flagVal $args -range 4]
    set remVar     [flagVal $args -remarksVariableName ""]
    set residueSel [flagVal $args -residueSelection [AtomSel -args "(all)"]]

    updateUser "Adding primary sequence neighbor exceptions \r"

    set residueList [residuesInSelection $residueSel]
    
    recordSAsInResidues $pot $residueList fromSAsInResidue toSAsInResidue
    recordSAsInNeighborhoods $PSDrange $residueList fromSAsInResidue toSAsInResidue fromSAsInNeighborhood toSAsInNeighborhood

    set nAdded [createPSDExceptions $pot $residueList fromSAsInResidue toSAsInResidue fromSAsInNeighborhood toSAsInNeighborhood]

    if {$remVar != ""} {
	
	upvar $remVar tempRem

	lappend tempRem [format "addPrimarySequenceExceptions:  with PSD range +/- %d residues, added %d inverse exceptions" \
			     $PSDrange $nAdded]
    }
    

    return ""
}

}

namespace import PrimarySequenceExceptions::*
