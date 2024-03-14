#
# TCL code for importing shift tables from 
# NMR-STAR files 
#

package require marvin
package provide nmrstar 1.0

#
# Given an input unit to an NMR-STAR file, read the 
# next "loop" header and return it in the form
# {{colTitle colPos} {colTitle colPos} ...}
#

proc readNMRSTARLoopHeader {inUnit} {

    if {[catch {eatNormalTextUpTo $inUnit "loop_"}]} {
	error "No more loops"
    }

    set retVal [list]
    set curColumnNum -1
    
    while {1} {

	set curLine [nextLine $inUnit]

	#
	# column headers end with a blank line
	#

	if {[string trim $curLine] == ""} {
	    break
	}
	
	incr curColumnNum
	
	set curColumnTitle [lindex $curLine 0]

	lappend retVal [list $curColumnTitle $curColumnNum]
    }

    return $retVal
}

#
# Search for a specific NMR-STAR loop (defined by a list of column titles) and 
# return the data in it or an error
#

proc readNMRSTARLoop {inUnit colTitleList} {

    global errorInfo

    #
    # start reading from the beginning of the file
    #

    seek $inUnit 0 start

    #
    # read loop headers until we find one with all the required column titles
    #

    while {1} {

	#
	# readNMRSTARLoopHeader returns an error when there are no more loops
	#

	if {[catch {set curHeader [readNMRSTARLoopHeader $inUnit]}]} {
	    error "Can't find loop header with required column titles $colTitleList"
	}

	#
	# see if this loop's header defines all the values 
	# we're looking for
	#

	set colPosList [list]

	foreach reqColTitle $colTitleList {

	    set reqColPos -1

	    foreach elem $curHeader {

		set curColTitle [lindex $elem 0]
		set curColPos   [lindex $elem 1]

		if {$curColTitle == $reqColTitle} {
		    set reqColPos $curColPos
		} 
	    }

	    lappend colPosList $reqColPos
	}

	set isTargetHeader 1
	foreach elem $colPosList {
	    if {$elem == -1} {
		set isTargetHeader 0
		break
	    }
	}

	if {$isTargetHeader} {
	    break
	}
    }


    #
    # we've found the proper loop.  Now read the data
    #

    set rawData [list]
    
    while {1} {
	
	if {[catch {set l [nextLine $inUnit]}]} {
	    set errorInfo ""
	    break
	}
	
	#
	# skip blank lines
	#
	
	if {[string trim $l] == ""} {
	    continue
	}
	
	#
	# data end at "stop_"
	#
	
	if {[string trim $l] == "stop_"} {
	    break
	}
	
	set curVal [list]
	
	#
	# try extracting out the columns' data
	#
	
	if {[catch {
	    
	    foreach pos $colPosList {
		lappend curVal [lindex $l $pos]
	    }
	    
	}]} {
	    
	    error "cant parse line: $l from file $fname"
	}
	
	lappend rawData $curVal
    }

    return $rawData
}


proc readNMRSTARShiftTable args {

    set fileName          [requiredFlagVal $args -fileName]
    set remVar            [flagVal $args -remarksVariableName ""]
    set useAmbiguityCodes [flagExists $args -useAmbiguityCodes]

    updateUser "Importing NMR-STAR formatted shift table \r"

    #
    # open the file
    #

    if { [catch {set inUnit [open $fileName r]}] } {
    	error "Error opening input file $fileName"
    }

    #
    # read the atom shift loop.  First try to find the "residue author seq 
    # code," which is the residue number if it didn't start from 1.  If it 
    # doesn't exist,  use the "residue seq code" as the residue number.
    #

    if {[catch {set rawShiftData \
		    [readNMRSTARLoop $inUnit \
			 [list "_Atom_shift_assign_ID" \
			      "_Residue_author_seq_code" \
			      "_Residue_label" \
			      "_Atom_name" \
			      "_Chem_shift_value" \
			      "_Chem_shift_ambiguity_code"]] }]} {
	

	set rawShiftData \
	    [readNMRSTARLoop $inUnit [list "_Atom_shift_assign_ID" \
					  "_Residue_seq_code" \
					  "_Residue_label" \
					  "_Atom_name" \
					  "_Chem_shift_value" \
					  "_Chem_shift_ambiguity_code"]]
    }
    
    #
    # The NMR-STAR ambiguity code is
    #
    # 1:  unique
    # 2:  ambiguity of geminal atoms or geminal methyl proton groups
    # 3:  ambiguity of aromatic atoms on opposite sides of the ring (eg., phe HD1 | HD2)
    # 4:  intraresidue ambiguities (eg., lys HG vs HD protons
    # 5:  interresidue ambiguities (eg., lys 12 vs lys 27)
    # 9:  other ambiguity
    #
    # see http://www.bmrb.wisc.edu/elec_dep/gen_aa.html
    #
    # If there are any entries with ambiguity codes 4 or 5, we need to find 
    # the corresponding ambiguity detail loop
    #

    set hasHighAmbiguity 0
    foreach elem $rawShiftData {
	set curAmbigCode [lindex $elem 5]
	if {($curAmbigCode == 4) || ($curAmbigCode == 5)} {
	    set hasHighAmbiguity 1
	    break
	}
    }

    
    if {$hasHighAmbiguity} {

	#
	# NMR-STAR shift ambiguity entries are a comma-delimited list of atom shift assign IDs.
	# Each is the assignment ID for a chemical shift assignment that has been given 
	# an ambiguity code of 4 or 5.  Each set indicates that the observed chemical shifts
	# are related to the defined atoms, but have not been assigned uniquely to a specific
	# atom in the set.
	# 
     
	set rawAmbiguityData [readNMRSTARLoop $inUnit [list "_Atom-shift_assign_ID_ambiguity"]]

	set ambiguitySets [list]
	foreach elem $rawAmbiguityData {
	    lappend ambiguitySets [split $elem ","]
	}
    }	
	
    #
    # convert the raw data into my shift list format:
    # { {shift sel sel ...} {shift sel sel ...} ... }
    #

    set shiftData [list]

    foreach elem $rawShiftData {

	set curAssignID  [lindex $elem 0]
	set curResNum    [lindex $elem 1]
	set curResName   [lindex $elem 2]
	set curAtomName  [lindex $elem 3]
	set curShift     [lindex $elem 4]
	set curAmbiguity [lindex $elem 5]

	if { [string length $curResName] == 1 } {
	    marvinPyth command {import selectTools}
	    set curResName [lindex \
			    [lindex \
			     [marvinPyth command \
		  {out=selectTools.renameResidues(inVar)[0]} \
			  [list [list inVar $curResName]] {out}] 0] 1]
	}

	#
	# ambiguity codes are often wrong, so default behavior is 
	# to make everything unique and let the standard 
	# non-stereo processing handle it
	#

	if {! $useAmbiguityCodes} {
	    set curAmbiguity 1
	}

	switch $curAmbiguity {

	    1 {set curSels [list [correctedSelection $curResNum $curResName $curAtomName]]}

	    2 {set curSels [geminalAmbigSels $curResNum $curResName $curAtomName]}
	    
	    3 {set curSels [aromaticAmbigSels $curResNum $curResName $curAtomName]}

	    4 {error "Can't yet deal with ambiguity type 4"}
	    
	    5 {error "Can't yet deal with ambiguity type 5" }

	    9 {error "Can't yet deal with ambiguity type 9"}
	}

	set curData [list $curShift]
	foreach sel $curSels {
	    lappend curData $sel
	}

	lappend shiftData $curData
    }

    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "%d shifts were read from NMR-STAR formatted file %s" [llength $shiftData] $fileName]
    }

	
    close $inUnit
    return $shiftData
}

proc readNMRSTAR args {

    set fileName          [requiredFlagVal $args -fileName]

    updateUser "Importing NMR-STAR information \r"

    #
    # open the file
    #

    if { [catch {set inUnit [open $fileName r]}] } {
    	error "Error opening input file $fileName"
    }

    marvinPyth command "import pasd"
    marvinPyth command "pasd.readSTAR('$fileName')"
}

proc readNMRSTARShifts args {

#    set fileName          [requiredFlagVal $args -fileName]
    set remVar            [flagVal $args -remarksVariableName ""]
    set useAmbiguityCodes [flagExists $args -useAmbiguityCodes]
    set saveSet           [flagVal $args -saveSet ""]

    updateUser "Importing NMR-STAR shift data \r"

    

    if { $useAmbiguityCodes } {
	set ret [marvinPyth command \
		     {ret=pasd.starShifts(True,saveSet=saveSet)} \
		     [list [list saveSet $saveSet]] \
		     {ret}]
    } else {
	set ret [marvinPyth command \
		     {ret=pasd.starShifts(False,saveSet=saveSet)} \
		     [list [list saveSet $saveSet]] \
		     {ret}]
    }

    set ret [lindex [lindex $ret 0] 1]

    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "%d NMR-STAR formatted shifts were read" \
			 [llength $ret]]
    }

    return $ret
}

#
# NMR-STAR formatted shift tables typically have three 
# types of atom name problems:
#
# 1.  Backbone amides are named H, rather than HN.  
#
# 2.  methylenes are typically named H*2 and H*3 instead of 
# H*1 and H*2.  Not sure how those names should be mapped to 
# each other.  For now, I'm just chaning H*3 into H*1 and leaving
# H*2 alone.  
#
# 3.  Selections involving methyls usually only select one
# of the three protons.  
#
# Returns an xplor selection with the atom name corrected/expanded 
# as necessary.
#

proc correctedSelection {resNum resName atomName} {


    # format for mappings is resType, nmrStarAtomName, xplorAtomName

    set mappings [list]

    lappend mappings [list "*" "h" "hn"]

    lappend mappings [list "gly" "ha3" "ha1"]
    lappend mappings [list "ile" "hg13" "hg11"]
    lappend mappings [list "leu" "hb3" "hb1"]
    lappend mappings [list "phe" "hb3" "hb1"]
    lappend mappings [list "trp" "hb3" "hb1"]
    lappend mappings [list "cys" "hb3" "hb1"]
    lappend mappings [list "ser" "hb3" "hb1"]
    lappend mappings [list "asn" "hb3" "hb1"]
    lappend mappings [list "tyr" "hb3" "hb1"]
    lappend mappings [list "his" "hb3" "hb1"]
    lappend mappings [list "asp" "hb3" "hb1"]
    lappend mappings [list "met" "hb3" "hb1"]
    lappend mappings [list "met" "hg3" "hg1"]
    lappend mappings [list "gln" "hb3" "hb1"]
    lappend mappings [list "gln" "hg3" "hg1"]
    lappend mappings [list "glu" "hb3" "hb1"]
    lappend mappings [list "glu" "hg3" "hg1"]
    lappend mappings [list "pro" "hb3" "hb1"]
    lappend mappings [list "pro" "hg3" "hg1"]
    lappend mappings [list "pro" "hd3" "hd1"]
    lappend mappings [list "arg" "hb3" "hb1"]
    lappend mappings [list "arg" "hg3" "hg1"]
    lappend mappings [list "arg" "hd3" "hd1"]
    lappend mappings [list "lys" "hb3" "hb1"]
    lappend mappings [list "lys" "hg3" "hg1"]
    lappend mappings [list "lys" "hd3" "hd1"]
    lappend mappings [list "lys" "he3" "he1"]

    lappend mappings [list "ala" "hb*"  "hb\#"]
    lappend mappings [list "val" "hg1*" "hg1\#"]
    lappend mappings [list "val" "hg2*" "hg2\#"]
    lappend mappings [list "ile" "hg2*" "hg2\#"]
    lappend mappings [list "ile" "hd*"  "hd*"]
    lappend mappings [list "leu" "hd1*" "hd1\#"]
    lappend mappings [list "leu" "hd2*" "hd2\#"]
    lappend mappings [list "met" "he*"  "he*"]
    lappend mappings [list "thr" "hg2*" "hg2\#"]

    foreach map $mappings {

	set mapResName [lindex $map 0]
	set mapOldAtomName [lindex $map 1]
	set mapNewAtomName [lindex $map 2]

	if {([atomNameMatches $mapResName $resName]) && 
	    ([atomNameMatches $mapOldAtomName $atomName])} {
	    set atomName $mapNewAtomName
	    break
	}
    }

    return [format "(resid %d and resn %s and name %s)" \
		$resNum [string tolower $resName] [string tolower $atomName]]
}
	



proc atomNameMatches {pattern atomName} {

    if {[string match [string tolower $pattern] [string tolower $atomName]]} {
	return 1
    } else {
	return 0
    }
}


proc aromaticAmbigSels {resNum resName atomName} {

    if {($resName != "PHE") && ($resName != "TYR")} {
	error "can't parse aromatic ambiguity resid $resNum resn $resName atomname $atomName"
    }

    if {[atomNameMatches hd* $atomName]} {
	return [list [correctedSelection $resNum $resName "hd*"]]

    } elseif {[atomNameMatches he* $atomName]} {
	return [list [correctedSelection $resNum $resName "he*"]]

    } elseif {[atomNameMatches cd* $atomName]} {
	return [list [correctedSelection $resNum $resName "cd*"]]

    } elseif {[atomNameMatches ce* $atomName]} {
	return [list [correctedSelection $resNum $resName "ce*"]]

    } else {
	error "can't parse aromatic ambiguity resid $resNum resn $resName atomname $atomName"
    }
}


proc geminalAmbigSels {resNum resName atomName} {

    switch $resName {

	"GLY" {

	    if {[atomNameMatches ha* $atomName]} {
		return [list [correctedSelection $resNum $resName "ha*"]]
	    } else {
		error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	    }
	}

	"ALA" {

	    error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	}

	"VAL" {

	    if {[atomNameMatches hg* $atomName]} {
		return [list [correctedSelection $resNum $resName "hg*"]]
	    }  elseif {[atomNameMatches cg* $atomName]} {
		return [list [correctedSelection $resNum $resName "cg*"]]
	    } else {
		error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	    }
	}

	"ILE" {

	    if {[atomNameMatches hg1* $atomName]} {
		return [list [correctedSelection $resNum $resName "hg1*"]]
	    } elseif {[atomNameMatches cg* $atomName]} {
		return [list [correctedSelection $resNum $resName "cg1"]]
	    } else {
		error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	    }
	}

	"LEU" {

	    if {[atomNameMatches hb* $atomName]} {
		return [list [correctedSelection $resNum $resName "hb*"]]
	    } elseif {[atomNameMatches hd* $atomName]} {
		return [list [correctedSelection $resNum $resName "hd*"]]
	    } elseif {[atomNameMatches cd* $atomName]} {
		return [list [correctedSelection $resNum $resName "cd*"]]
	    } else {
		error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	    }
	}


	"PRO" {

	    if {[atomNameMatches hb* $atomName]} {
		return [list [correctedSelection $resNum $resName "hb*"]]
	    } elseif {[atomNameMatches hg* $atomName]} {
		return [list [correctedSelection $resNum $resName "hg*"]]
	    } elseif {[atomNameMatches hd* $atomName]} {
		return [list [correctedSelection $resNum $resName "hd*"]]
	    } else {
		error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	    }
	}


	"THR" {

	    error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	}
	    

	"GLU" -
	"MET" {

	    if {[atomNameMatches hb* $atomName]} {
		return [list [correctedSelection $resNum $resName "hb*"]]
	    } elseif {[atomNameMatches hg* $atomName]} {
		return [list [correctedSelection $resNum $resName "hg*"]]
	    } else {
		error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	    }
	}

	"GLN" {

	    if {[atomNameMatches hb* $atomName]} {
		return [list [correctedSelection $resNum $resName "hb*"]]
	    } elseif {[atomNameMatches hg* $atomName]} {
		return [list [correctedSelection $resNum $resName "hg*"]]
	    } elseif {[atomNameMatches he2* $atomName]} {
		return [list [correctedSelection $resNum $resName "he2*"]]
	    } else {
		error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	    }
	}


	"ASN" {

	    if {[atomNameMatches hb* $atomName]} {
		return [list [correctedSelection $resNum $resName "hb*"]]
	    } elseif {[atomNameMatches hd2* $atomName]} {
		return [list [correctedSelection $resNum $resName "hd2*"]]
	    } else {
		error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	    }
	}


	"PHE" -
	"TRP" -
	"CYS" -
	"SER" -
	"TYR" -
	"HIS" -
	"ASP" {

	    if {[atomNameMatches hb* $atomName]} {
		return [list [correctedSelection $resNum $resName "hb*"]]
	    } else {
		error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	    }
	}

	"LYS" {

	    if {[atomNameMatches hb* $atomName]} {
		return [list [correctedSelection $resNum $resName "hb*"]]
	    } elseif {[atomNameMatches hg* $atomName]} {
		return [list [correctedSelection $resNum $resName "hg*"]]
	    } elseif {[atomNameMatches hd* $atomName]} {
		return [list [correctedSelection $resNum $resName "hd*"]]
	    } elseif {[atomNameMatches he* $atomName]} {
		return [list [correctedSelection $resNum $resName "he*"]]
	    } else {
		error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	    }
	}


	"ARG" {

	    if {[atomNameMatches hb* $atomName]} {
		return [list [correctedSelection $resNum $resName "hb*"]]
	    } elseif {[atomNameMatches hg* $atomName]} {
		return [list [correctedSelection $resNum $resName "hg*"]]
	    } elseif {[atomNameMatches hd* $atomName]} {
		return [list [correctedSelection $resNum $resName "hd*"]]
	    } else {
		error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	    }
	}

	default { 
	    error "can't parse geminal ambiguity resid $resNum resn $resName atomname $atomName"
	}
    }
}

#
# support for NMR-STAR 2.1 format
#

proc process3dCNMRSTARPeaks args {

    set remVar [flagVal $args -remarksVariableName ""]
    if {$remVar != ""} {
	upvar $remVar $remVar
    }

    if {! [flagExists $args -namePrefix]} {
	lappend args "-namePrefix" "3dc"
    }

    eval process3dNMRSTARPeaks $args
    return ""
}
	
proc process3dNMRSTARPeaks args {

    set pot        [requiredFlagVal $args -pot]
    set H1name     [requiredFlagVal $args -fromProtonColumnName]
    set C1name     [requiredFlagVal $args [list -fromHeavyatomColumnName -fromCarbonColumnName -fromNitrogenColumnName]]
    set H2name     [requiredFlagVal $args -toProtonColumnName]
    set saveSet    [flagVal $args -saveSet ""]
    set namePrefix [flagVal $args -namePrefix "3d"]
    set remVar     [flagVal $args -remarksVariableName ""]
     
    set peaks [marvinPyth command \
     "ret=pasd.starPeaks(h1Name,h2Name,c1Name,saveSet=saveSet,tclOutput=True)"\
		   [list [list saveSet $saveSet] \
			[list h1Name $H1name] \
			[list h2Name $H2name] \
			[list c1Name $C1name]] {ret}]
    set peaks [lindex [lindex $peaks 0] 1]
    set nPeaksAdded 0

    puts [llength $peaks]

    foreach peak $peaks {

	set curID   [lindex $peak 0]
	set curName [format "%s%s" $namePrefix $curID]

	#
	# PIPP PCK/DEG tables sometimes have multiple lines referring to a single peak
	# with multiple possible peakAssignments.  Ignore them.
	#

	if {[$pot hasPeakNamed $curName]} {
	    continue
	}

	set tempPeak [Peak -args $curName]
	$tempPeak setIntensity [lindex $peak 1]
	$tempPeak setFromProtonShift [lindex $peak 2]
	$tempPeak setToProtonShift [lindex $peak 3]
	$tempPeak setFromHeavyatomShift [lindex $peak 4]
	$tempPeak appendToNote [format "from set %s, peak %s" $saveSet $curID]
	
	$pot addPeak [$tempPeak cget -this]
	$tempPeak -disown
	rename $tempPeak ""
	incr nPeaksAdded
    }

    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "Created %d peaks from set %s" \
			     $nPeaksAdded $saveSet]
    }
    return ""
}

proc process2dNMRSTARPeaks args {

    set pot        [requiredFlagVal $args -pot]
    set H1name     [requiredFlagVal $args -fromProtonColumnName]
    set H2name     [requiredFlagVal $args -toProtonColumnName]
    set saveSet    [flagVal $args -saveSet ""]
    set namePrefix [flagVal $args -namePrefix "2d"]
    set remVar     [flagVal $args -remarksVariableName ""]
     
    set peaks [marvinPyth command \
     "ret=pasd.starPeaks(h1Name,h2Name,saveSet=saveSet,tclOutput=True)" \
		   [list [list saveSet $saveSet] \
			[list h1Name $H1name] \
			[list h2Name $H2name]] {ret}]
    set peaks [lindex [lindex $peaks 0] 1]
    set nPeaksAdded 0

    puts [llength $peaks]

    foreach peak $peaks {

	set curID   [lindex $peak 0]
	set curName [format "%s%s" $namePrefix $curID]

	#
	# PIPP PCK/DEG tables sometimes have multiple lines referring to a single peak
	# with multiple possible peakAssignments.  Ignore them.
	#

	if {[$pot hasPeakNamed $curName]} {
	    continue
	}

	set tempPeak [Peak -args $curName]
	$tempPeak setIntensity [lindex $peak 1]
	$tempPeak setFromProtonShift [lindex $peak 2]
	$tempPeak setToProtonShift [lindex $peak 3]
	$tempPeak appendToNote [format "from set %s, peak %s" $saveSet $curID]
	
	$pot addPeak [$tempPeak cget -this]
	$tempPeak -disown
	rename $tempPeak ""
	incr nPeaksAdded
    }

    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "Created %d peaks from set %s" \
			     $nPeaksAdded $saveSet]
    }
    return ""
}

# don't know how to do this
#proc readNMRSTARSpectralRanges args {
#    set saveSet    [flagVal $args -saveSet "c13noesy"]
#
#    set widths [marvinPyth command \
#     "ret=pasd.starSpectralWidths3d(saveSet,tclOutput=True)" \
#		   [list [list saveSet $saveSet]] {ret}]
#
#    set widths [lindex [lindex $widths 0] 1]
#
#    return widths
