#
# pdb_support.tcl
#
# Routines to handle writing PDB files
#
# JJK 5/21/02
#

package provide marvin 1.0
package require pdbtool

proc writePDB args {

    #
    # grab the flags, enforcing defaults
    #

    set outName     [requiredFlagVal $args -fileName]
    set remarksList [flagVal $args -remarks]
    set aux1        [flagVal $args -aux1]
    set aux2        [flagVal $args -aux2]

    set needToCleanSel 0
    set aSel [flagVal $args -selection]
    if {$aSel == {}} {	
	set needToCleanSel 1
	set aSel [AtomSel -args "all"]
    }

    set pt [PDBTool]
    $pt setSelection [$aSel cget -this]
    $pt setFilename $outName
    $pt clearRemarks

    foreach rem $remarksList {
	foreach l [split $rem \n] {
	    $pt addRemark $l
	}
    }

    if {[llength $aux1] == [$aSel size]} {
	for {set x 0} {$x < [$aSel size]} {incr x} {
	    $pt setAux1 [$aSel atomAtPos $x] [lindex $aux1 $x]
	}
    }

    if {[llength $aux2] == [$aSel size]} {
	for {set x 0} {$x < [$aSel size]} {incr x} {
	    $pt setAux2 [$aSel atomAtPos $x] [lindex $aux2 $x]
	}
    }

    if {[catch {$pt write}]} {

	error "Error writing PDB file $outName"
    }

    rename $pt ""

    if {$needToCleanSel} {
	rename $aSel ""
    }
}


proc readPDB args {

    #
    # grab the flags, enforcing defaults
    #

    set inName   [requiredFlagVal $args -fileName]
    set aSimPtr  [flagVal $args -simulation [PtrToCurrentSimulation]]
    set modelNum [flagVal $args -model]

    set needToCleanSel 0
    set aSel [flagVal $args -selection]
    if {$aSel == {}} {	
	set needToCleanSel 1
	set aSel [AtomSel -args "all"]
    }

    set inNames [glob $inName] 
    if {[llength $inNames] == 0} {
	error "Error reading PDB file $inName.  File does not exist."
    } elseif {[llength $inNames] > 1} {
	error "Error reading PDB file $inName.  >1 file matches pattern."
    } else {
	set inName [lindex $inNames 0]
    }

    set pt [PDBTool]
    $pt setFilename $inName
    $pt setSelection [$aSel cget -this]

    if {$modelNum != ""} {

	if {[catch {$pt read $modelNum}]} {
	    error "Error reading PDB file $inName model number $modelNum"
	}

    } else {

	if {[catch {$pt read}]} {    
	    error "Error reading PDB file $inName"
	}
    }
    
    set retVal [$pt remarks]

    rename $pt ""

    if {$needToCleanSel} {
	rename $aSel ""
    }

    return $retVal
}


#
# Given a PDB file name, return the remarks portion.
# Stops reading when it hits an "atom" line in the PDB
#

proc readPDBremarks args {

    set fileName [requiredFlagVal $args -fileName]

    global errorInfo
    
    if {[catch {set inUnit [open $fileName r]}]} {
    
	error "Error opening input file $fileName"
    }

    set retVal {}

    while {1} {

	if {[catch {set curLine [nextLine $inUnit]}]} {
	    set errorInfo ""
	    break
	}

	set firstWord [string tolower [lindex $curLine 0]]

	if {[string match "atom*" $firstWord]} {
	    set errorInfo ""
	    break
	} elseif {[string match "rema*" $firstWord]} {
	    set rest [lrange $curLine 1 end]
	    lappend retVal $rest
	}
    }

    close $inUnit

    return $retVal
}

