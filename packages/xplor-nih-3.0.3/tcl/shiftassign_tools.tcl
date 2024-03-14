#
# procedures to read, write, and generate shift assignments
#
# JJK 6/10/05
#

package provide marvin 1.0
package require pasdpot

proc readShiftAssignments args {
    
    set fileName   [requiredFlagVal $args -fileName]
    set potential  [requiredFlagVal $args -pot]
    
    global errorInfo
    
    if { [catch {set inUnit [open $fileName r]}] } {
    	error "Error opening input file $fileName"
    }
    
    XplorCommand "set message off echo off end"
    
    # read the shift assignments
    
    while {![eof $inUnit]} {
	
	if { [catch {set curShiftAssign [readOneShiftAssign $inUnit]} msg]} {
	    	    
	    if {$msg == "No more shift assignments"} {
		set errorInfo ""
		break
	    } else {
		set savedInfo $errorInfo
		error [format "Shift assignment reading error in file %s:\n %s" $fileName $msg] $savedInfo
	    }
	}
	
	processOneShiftAssign $curShiftAssign $potential
	updateUser [format "read shift assignment %s \r" [lindex $curShiftAssign 0]]
    }
    
    
    close $inUnit
}


proc readOneShiftAssign {fileID} {

    # 
    # Try to find the beginning of another shift assignment,
    # ignoring text in xplor-style comments 
    #
    
    global errorInfo 

    if {[catch {eatXplorTextUpTo $fileID "shiftAssignment"}]} {

	set errorInfo ""
	error "No more shift assignments"
    }

    # defaults for optional fields 
    
    set curNote ""
    set protonShift "NONE"
    set heavyatomShift "NONE"
    set foldedProtonShift "NONE"
    set foldedHeavyatomShift "NONE"
    set prevLikelihood "NONE"
    set protonSel "NONE"
    set heavyatomSel "NONE"
    set whichKind ""
    set toFromPartnerName "NONE"
    set stereoPartnerName "NONE"
    set isMethyl 0
    set isGood 0
    set linkedAssignList [list]

    # grab the required fields

    set curName [nextWord $fileID]
    if {$curName == ""} {
	error "Missing shift assignment name"
    }
    
    # search for optional fields

    set done 0
    while {! $done} {
	
	set curWord [nextWord $fileID]
	
	switch -exact -- $curWord {
	    
	    "" {
		error [format "unexpected end of file at shift assignment %s" $curName]
	    }
	    
	    "end" {
		set done 1
	    }
	    
	    "shiftAssignment" {
		error [format "missing end statement after shift assignment %s" $curName]
	    }
	    
	    "-protonShift" {

		set protonShift [nextWord $fileID]

		if {! [isNumber $protonShift]} {
		    error [format "Missing or invalid proton shift for shift assignment %s" $curName]
		}
	    }
	   
	    "-heavyatomShift" {

		set heavyatomShift [nextWord $fileID]

		if {! [isNumber $heavyatomShift]} {
		    error [format "Missing or invalid heavyatom shift for shift assignment %s" $curName]
		}
	    }

	    "-foldedProtonShift" {

		set foldedProtonShift [nextWord $fileID]

		if {! [isNumber $foldedProtonShift]} {
		    error [format "Missing or invalid folded proton shift for shift assignment %s" $curName]
		}
	    }
	   
	    "-foldedHeavyatomShift" {

		set foldedHeavyatomShift [nextWord $fileID]

		if {! [isNumber $foldedHeavyatomShift]} {
		    error [format "Missing or invalid folded heavyatom shift for shift assignment %s" $curName]
		}
	    }
	   	   
	    "-note" {
		
		if {$curNote == ""} {
		    set curNote [nextLine $fileID]
		} else {
		    set curNote [format "%s\n%s" $curNote [nextLine $fileID]]
		}
	    }

	    "-protonSelection" {
		set protonSel [nextSel $fileID]
		if {$protonSel == ""} {
		    error [format "Missing or invalid proton selection for shift assignment %s" $curName]
		}
	    }

	    "-heavyatomSelection" {
		set heavyatomSel [nextSel $fileID]
		if {$heavyatomSel == ""} {
		    error [format "Missing or invalid heavyatom selection for shift assignment %s" $curName]
		}
	    }

	    "-previousLikelihood" {
		set prevLikelihood [nextWord $fileID]
		
		if {! [isZeroToOneNumber $prevLikelihood]} {
		    error [format "Missing or invalid previous likelihood for shift assignment %s" $curName]
		}
	    }

	    
	    "-symmetricShiftAssignment" -
	    "-toFromPartner" {
		
		set toFromPartnerName [nextWord $fileID]
	    }

	    "-stereoPartner" {
		
		set stereoPartnerName [nextWord $fileID]
	    }

	    
	    "-from" {
		set whichKind "from"
	    }

	    "-to" {
		set whichKind "to"
	    }

	    "-methyl" {
		set isMethyl 1
	    }

	    "-good" {
		set isGood 1
	    }

	    "-linkedTo" {
		set curLinkedAssName [nextWord $fileID]
		lappend linkedAssignList $curLinkedAssName
	    }
	}
    }
    
    #
    # make sure we actually read at least a proton selection and a from/to flag
    #
    
    if {$protonSel == "NONE"} {
	error [format "No proton selection found for shift assignment %s" $curName]
    }

    if {$whichKind == ""} {
	error [format "No from or to flag found for shift assignment %s" $curName]
    }

    return [list $curName $protonShift $heavyatomShift $protonSel $heavyatomSel \
		$prevLikelihood $isMethyl $isGood $whichKind $toFromPartnerName $stereoPartnerName $curNote $linkedAssignList $foldedProtonShift $foldedHeavyatomShift]
}


proc processOneShiftAssign {shiftAssignData pot} {
    
    #
    # takes the data read by readOneShiftAssign
    # and actually creates a ShiftAssignment and 
    # attaches it to a MarvinNOEPotential
    #
    
    set assignName         [lindex $shiftAssignData 0]
    set protonShift        [lindex $shiftAssignData 1]
    set heavyatomShift     [lindex $shiftAssignData 2]
    set protonSelString    [lindex $shiftAssignData 3]
    set heavyatomSelString [lindex $shiftAssignData 4]
    set previousLikelihood [lindex $shiftAssignData 5]
    set isMethyl           [lindex $shiftAssignData 6]
    set isGood             [lindex $shiftAssignData 7]
    set whichKind          [lindex $shiftAssignData 8]
    set toFromPartnerName  [lindex $shiftAssignData 9]
    set stereoPartnerName  [lindex $shiftAssignData 10]
    set note               [lindex $shiftAssignData 11]
    set linkedAssignList   [lindex $shiftAssignData 12]
    set foldedProtonShift  [lindex $shiftAssignData 13]
    set foldedHeavyShift   [lindex $shiftAssignData 14]
    
    #
    # Create the shift assignment
    #
    
    set anAssign [ShiftAssignment -args $assignName]
     
    $anAssign setNote $note
    
    if {$protonShift != "NONE"} {
	$anAssign setProtonShift $protonShift
    }
    
    if {$heavyatomShift != "NONE"} {
	$anAssign setHeavyatomShift $heavyatomShift
    }

    if {$foldedProtonShift != "NONE"} {
	$anAssign setFoldedProtonShift $foldedProtonShift
    }
    
    if {$foldedHeavyShift != "NONE"} {
	$anAssign setFoldedHeavyatomShift $foldedHeavyShift
    }

    if {$previousLikelihood != "NONE"} {
	$anAssign setPreviousLikelihood $previousLikelihood
    }

    if {$toFromPartnerName != "NONE"} {
	$anAssign setToFromPartnerName $toFromPartnerName
    }

    if {$stereoPartnerName != "NONE"} {
	$anAssign setStereoPartnerName $stereoPartnerName
    }

    if {$whichKind == "from"} {
	$anAssign setIsFrom
    } else {
	$anAssign setIsTo
    }

    if {$isMethyl} {
	$anAssign setIsMethyl
    } else {
	$anAssign setNotMethyl
    }

    if {$isGood} {
	$anAssign setGood
    } else {
	$anAssign setBad
    }

    foreach elem $linkedAssignList {
	$anAssign addLinkedShiftAssignmentName $elem
    }
    
    # note that these selections are attached to the currentSimulation
    
    if {$protonSelString != "NONE"} {
	$anAssign setProtonSelectionString $protonSelString
    }

    if {$heavyatomSelString != "NONE"} {
	$anAssign setHeavyatomSelectionString $heavyatomSelString
    }
       
    #
    # add it to the pot
    #
    
    $pot addShiftAssignment [$anAssign cget -this]
    $anAssign -disown
    rename $anAssign ""
}


proc writeShiftAssignments args {

    set fileName     [requiredFlagVal $args -fileName]
    set remarksList  [flagVal $args -remarks]
    set shiftAssigns [requiredFlagVal $args -shiftAssignments]
    
    if {[catch {set outUnit [open $fileName w]}]} {
	
	error "Error opening output file $fileName"
    }

    #
    # write remarks
    #
    
    if {[llength $remarksList] > 0} {
	puts $outUnit [format "!"]
	foreach rem $remarksList {
	    foreach l [split $rem "\n"] {
		puts $outUnit [format "! %s" $l]
	    }
	    puts $outUnit [format "!"]
	}
    }
    
    #
    # write each shift assignment
    #

    foreach assign $shiftAssigns {
	writeOneShiftAssign $outUnit $assign
    }

    close $outUnit
}


proc writeOneShiftAssign {outUnit anAssign} {
    
    ShiftAssignment -this $anAssign
    
    puts $outUnit [format "shiftAssignment %s " [$anAssign name]]

    puts $outUnit [format "   -protonSelection (%s)" [$anAssign protonSelectionString]]

    # heavyatom selection is optional, so don't write anything if it doesn't exist
    
    if {[$anAssign hasHeavyatomSelection]} {
	puts $outUnit [format "   -heavyatomSelection (%s)" [$anAssign heavyatomSelectionString]]
    }

    if {[$anAssign isFrom]} {
	puts $outUnit "   -from"
    } else {
	puts $outUnit "   -to"
    }

    if {[$anAssign isMethyl]} {
	puts $outUnit "   -methyl"
    }

    if {[$anAssign isGood]} {
	puts $outUnit "   -good"
    }

    if {[$anAssign hasProtonShift]} {
	puts $outUnit [format "   -protonShift %f" [$anAssign protonShift]]
    }

    if {[$anAssign hasHeavyatomShift]} {
	puts $outUnit [format "   -heavyatomShift %f" [$anAssign heavyatomShift]]
    }

    if {[$anAssign hasFoldedProtonShift]} {
	puts $outUnit [format "   -foldedProtonShift %f" [$anAssign foldedProtonShift]]
    }

    if {[$anAssign hasFoldedHeavyatomShift]} {
	puts $outUnit [format "   -foldedHeavyatomShift %f" [$anAssign foldedHeavyatomShift]]
    }

    if {[$anAssign hasPreviousLikelihood]} {
	puts $outUnit [format "   -previousLikelihood %f" [$anAssign previousLikelihood]]
    }

    if {[$anAssign hasToFromPartnerName]} {
	puts $outUnit [format "   -toFromPartner %s" [$anAssign toFromPartnerName]]
    }

    if {[$anAssign hasStereoPartnerName]} {
	puts $outUnit [format "   -stereoPartner %s" [$anAssign stereoPartnerName]]
    }

    foreach elem [$anAssign linkedShiftAssignmentNames] {
	puts $outUnit [format "   -linkedTo %s" $elem]
    }

    if {[$anAssign note] != ""} {
	foreach l [split [$anAssign note] "\n"] {
	    puts $outUnit [format "   -note %s" $l]
	}
    }

    puts $outUnit "end"
    puts $outUnit ""

    rename $anAssign ""
}

proc consolidateIdenticalStereopartners args {

    set pot    [requiredFlagVal $args -pot]
    set remVar [flagVal $args -remarksVariableName ""]
    global errorInfo

    #
    # Lots of shift tables have entries for an averaged stereopair that 
    # look like two stereospecifically assigned protons that happen to have
    # the same chemical shift.  After non-stereo processing, I won't be able to 
    # tell them apart, and will end up with multiple to-from partners.  
    #
    # Eliminate these extra copies of shift assignments prior to non-stereo
    # processing
    #
    # This is run directly after createShiftAssignments, before to-from
    # partners are set or non-stereo processing is done.
    #

    set rpt [list]

    set shiftAssigns [$pot shiftAssignments]
    set saNames [list]
    foreach curSA $shiftAssigns {
	ShiftAssignment -this $curSA
	lappend saNames [$curSA name]
	rename $curSA ""
    }

    set assignCount 0

    foreach curSAname $saNames {

	if {! [$pot hasShiftAssignmentNamed $curSAname]} {
	    continue
	}
	
	set curSA [$pot shiftAssignmentNamed $curSAname]

	updateUser [format \
			"Consolidating identical stereopartners %d of %d \r" \
			[incr assignCount] [llength $saNames]]
	
	if {[catch {set partner [findStereopartnerForShiftAssignment \
			     -shiftAssignmentList [$pot shiftAssignments] \
			     -targetShiftAssignment $curSA] } msg]} {

	    if {$msg == "No possible partner"} {
		set errorInfo ""
		continue
	    } else {
		set savedInfo $errorInfo
		error $msg $savedInfo
	    }
	}

	if {$partner == ""} {
	    continue
	}

	ShiftAssignment -this $curSA
	ShiftAssignment -this $partner

	if {[$curSA hasProtonShift] && [$partner hasProtonShift]} {
	    if {[$curSA protonShift] == [$partner protonShift]} {

		$pot removeShiftAssignmentNamed $curSAname
		$curSA -acquire

		lappend rpt [list [$curSA name] [$partner name]]
	    }
	}

	rename $curSA ""
	rename $partner ""
    }

    if {$remVar != ""} {
	
	upvar $remVar tempRem

	set newRem [format "Consolidated %d pairs of identical stereopartners:" [llength $rpt]]
	foreach elem $rpt {
	    appendLineToString newRem [format "   %s  --> %s" [lindex $elem 0] [lindex $elem 1]]
	}

	lappend tempRem $newRem
    }
    return ""
}




proc recordStereoPartnersForShiftAssignments args {

    global errorInfo
    
    set pot    [requiredFlagVal $args -pot]
    set remVar [flagVal $args -remarksVariableName ""]
    
    set shiftAssignments [$pot shiftAssignments]
    
    set numWith 0
    set count 0
    
    foreach curAssign $shiftAssignments {
	ShiftAssignment -this $curAssign
	$curAssign resetStereoPartnerName
	rename $curAssign ""
    }
    
    foreach curSA [$pot shiftAssignments] {
		
	updateUser [format \
		"Determining stereo partner for shiftAssign %d of %d \r" \
			[incr count] [llength $shiftAssignments]]
	
	if {[catch {set partner [findStereopartnerForShiftAssignment \
			     -shiftAssignmentList [$pot shiftAssignments] \
				     -targetShiftAssignment $curSA] } msg]} {

	    if {$msg == "No possible partner"} {
		set errorInfo ""
		continue
	    } else {
		set savedInfo $errorInfo
		error $msg $savedInfo
	    }
	}

	if {$partner != ""} {
	    
	    ShiftAssignment -this $curSA 
	    ShiftAssignment -this $partner
	    $curSA setStereoPartnerName [$partner name]
	    rename $partner ""
	    rename $curSA ""
	    incr numWith
	}
	
    }
    
    if {$remVar != ""} {
	
	upvar $remVar tempRem
	lappend tempRem "recordStereoPartnersForShiftAssignments: $numWith shiftAssignments had partners"
    }
    return ""
}



proc recordToFromPartnersForShiftAssignments args {

    set pot    [requiredFlagVal $args -pot]
    set remVar [flagVal $args -remarksVariableName ""]

    set shiftAssignments [$pot shiftAssignments]

    set numWith 0
    set numWithout 0
    set count 0

    foreach curAssign $shiftAssignments {
	ShiftAssignment -this $curAssign
	$curAssign resetToFromPartnerName
	rename $curAssign ""
    }
    
    foreach curAssign $shiftAssignments {
	
	updateUser [format "Determining to-from partner for shiftAssign %d of %d \r" \
			[incr count] [llength $shiftAssignments]]
	
	set partnerName [findToFromPartnerForShiftAssignment \
			     -shiftAssignmentList $shiftAssignments \
			     -targetShiftAssignment $curAssign]
	
	ShiftAssignment -this $curAssign
	
	if {$partnerName == ""} {
	    incr numWithout
	    $curAssign appendToNote "No to-from partner found"

	} else {

	    incr numWith
	    $curAssign setToFromPartnerName $partnerName
	}

	rename $curAssign ""
    }

    if {$remVar != ""} {
	
	upvar $remVar tempRem
	lappend tempRem "recordToFromPartnersForShiftAssignments: $numWith shiftAssignments had partners"
	lappend tempRem "recordToFromPartnersForShiftAssignments: $numWithout shiftAssignments had no partners"
    }
    return ""
}




proc findToFromPartnerForShiftAssignment args {

    set shiftAssigns [requiredFlagVal $args -shiftAssignmentList]
    set target       [requiredFlagVal $args -targetShiftAssignment]

    ShiftAssignment -this $target
    set targetProtSelString [$target protonSelectionString]
    set targetIsFrom        [$target isFrom]

    if {[$target hasProtonShift]} {
	set targetProtShift [$target protonShift]
    } else {
	set targetProtShift "NONE"
    }

    rename $target ""

    set targetProtSel [AtomSel -args $targetProtSelString]

    #
    # scan the list of shiftAssignments for ones with opposite 
    # from/to sense, equal proton shift, and an equal proton selection
    #

    set matchList [list]

    foreach otherAssign $shiftAssigns {
	
	ShiftAssignment -this $otherAssign
	set otherProtSel [$otherAssign protonSelection]
	set otherIsFrom  [$otherAssign isFrom]
	set otherName    [$otherAssign name]

	if {[$otherAssign hasProtonShift]} {
	    set otherProtShift [$otherAssign protonShift]
	} else {
	    set otherProtShift "NONE"
	}

	rename $otherAssign ""
	
	AtomSel -this $otherProtSel
	    
	if {($targetIsFrom != $otherIsFrom) && 
	    ($targetProtShift == $otherProtShift) &&
	    [$targetProtSel isEqualTo [$otherProtSel cget -this]]} {

	    lappend matchList $otherName
	    
	}
	
	rename $otherProtSel ""
    }
    
    rename $targetProtSel ""

    if {[llength $matchList] == 0} {
	return ""
    } elseif {[llength $matchList] == 1} {
	return [lindex $matchList 0]
    } else {
	ShiftAssignment -this $target
	set targetName [$target name]
	rename $target ""
	puts [format "%d to-from partners found for shift assignment %s: %s" [llength $matchList] $targetName $matchList]
	return [lindex $matchList 0]
    }
}


#
# Given a list of shiftAssigns and a target, 
# find a shiftAssign with the same from/to sense
# that's the target's stereopartner.
# Throws an error if the target is either not part
# of a stereopair or already covers an entire stereopair.
#

proc findStereopartnerForShiftAssignment args {

    set shiftAssigns [requiredFlagVal $args -shiftAssignmentList]
    set target       [requiredFlagVal $args -targetShiftAssignment]

    ShiftAssignment -this $target
    set targetName          [$target name]
    set targetProtSelString [$target protonSelectionString]
    set targetIsFrom        [$target isFrom]
    rename $target ""
	       
    # construct the stereopartner's proton selection
    
    set partnerProtSelString  [stereoPartnerForSelection $targetProtSelString]

    # If there's no change to the selection strings, this 
    # shift assignment is either not a part of a stereopair, or 
    # already covers the entire stereopair.  
    # Either way, throw an error
	
    if {$partnerProtSelString == $targetProtSelString} {
	error "No possible partner"
    }
    
    set partnerProtSel  [AtomSel -args $partnerProtSelString]

    #
    # see if this stereopartner selection already exists
    # in another shift assignment of the same from/to kind
    #

    set partnerShiftAssign ""
	
    foreach otherAssign $shiftAssigns {

	ShiftAssignment -this $otherAssign
	    
	if {$targetIsFrom == [$otherAssign isFrom]} {
   
	    if {[$partnerProtSel isEqualTo [$otherAssign protonSelection]]} {
		    
		set partnerShiftAssign $otherAssign
	    }
	}

	rename $otherAssign "" 
	    
	if {$partnerShiftAssign != ""} {
	    break
	}
    }
    
    return $partnerShiftAssign
}



proc expandHalfStereopairShiftAssignments args {

    #
    # Given a set of shiftAssignments, 
    #
    # Look for shiftAssignments which represent only one member of a
    # stereopair.  Make sure there isn't a second shiftAssignment
    # defining the other member.  Expand that single-member
    # shiftAssignment's selection to cover the entire stereopair.
    #
    # This is very conservative, in that it doesn't erase any existing
    # stereoassigns.  It's intended to handle situations like "phe
    # hd1" instead of "phe hd#"
    #

    set shiftAssigns [requiredFlagVal $args -shiftAssignments]
    set remVar       [flagVal $args -remarksVariableName ""]


    global errorInfo 

    set assignCount 0
    set numExpanded 0
    
    foreach curAssign $shiftAssigns {
	
	updateUser [format "Looking for stereopartner of shift assign %d of %d \r" \
			[incr assignCount] [llength $shiftAssigns]]
       
	if {[catch {set partner [findStereopartnerForShiftAssignment \
				     -shiftAssignmentList $shiftAssigns \
				     -targetShiftAssignment $curAssign] } msg]} {

	    if {$msg == "No possible partner"} {
		set errorInfo ""
		continue
	    } else {
		set savedInfo $errorInfo
		error $msg $savedInfo
	    }
	}

	#
	# if I didn't find a shiftAssignment with the stereopartner's selections, 
	# expand the original shiftAssignment's selections to cover both stereopartners
	#

	if {$partner == ""} {

	    ShiftAssignment -this $curAssign
	    
	    $curAssign setProtonSelectionString    [makeSelectionNonStereoSpecific [$curAssign protonSelectionString]]
	    if {[$curAssign hasHeavyatomSelection]} {
		$curAssign setHeavyatomSelectionString \
		    [makeSelectionNonStereoSpecific [$curAssign heavyatomSelectionString]]
	    }

	    $curAssign appendToNote "selections expanded by expandHalfStereopairShiftAssignments"
	    
	    rename $curAssign ""
	    
	    incr numExpanded
	}
    }

    if {$remVar != ""} {
	
	upvar $remVar tempRem
	lappend tempRem [format "%d of %d shiftAssignments were solitary members of stereopairs and were expanded" \
			     $numExpanded [llength $shiftAssigns]]
    }

    return ""
}


proc expandStereospecificShiftAssignments args {

    #
    # Given a set of shiftAssignments, 
    #
    # Expand the selections of each shiftAssignment that
    # corresponds to half a stereopair to cover the entire
    # stereopair.
    #
    # Note that this includes all functionality of
    # expandHalfStereopairShiftAssignments
    #

    set shiftAssigns [requiredFlagVal $args -shiftAssignments]
    set remVar       [flagVal $args -remarksVariableName ""]

    global errorInfo

    set assignCount 0
    set numExpanded 0
    
    foreach curAssign $shiftAssigns {
	
	updateUser [format \
		"Looking for stereopartner of shift assign %d of %d \r" \
			[incr assignCount] [llength $shiftAssigns]]
	
	if {[catch {set partner [findStereopartnerForShiftAssignment \
				   -shiftAssignmentList $shiftAssigns \
				   -targetShiftAssignment $curAssign] } msg]} {

	    if {$msg == "No possible partner"} {
		set errorInfo ""
		continue
	    } else {
		set savedInfo $errorInfo
		error $msg $savedInfo
	    }
	}

	ShiftAssignment -this $curAssign
	
	$curAssign setProtonSelectionString    [makeSelectionNonStereoSpecific [$curAssign protonSelectionString]]
	if {[$curAssign hasHeavyatomSelection]} {
	    $curAssign setHeavyatomSelectionString \
		[makeSelectionNonStereoSpecific [$curAssign heavyatomSelectionString]]
	}

	$curAssign appendToNote "selections expanded by expandStereospecificShiftAssignments"
	rename $curAssign ""
	incr numExpanded

	if {$partner != ""} {

	    ShiftAssignment -this $partner
	    
	    $partner setProtonSelectionString    [makeSelectionNonStereoSpecific [$partner protonSelectionString]]
	    if {[$partner hasHeavyatomSelection]} {
		$partner setHeavyatomSelectionString [makeSelectionNonStereoSpecific [$partner heavyatomSelectionString]]
	    }
	    
	    $partner appendToNote "selections expanded by expandStereospecificShiftAssignments"

	    rename $partner ""
	    incr numExpanded
	}
    }
    
    if {$remVar != ""} {
	
	upvar $remVar tempRem
	lappend tempRem \
	    [format "%d of %d stereospecific shiftAssignments were expanded by expandStereospecificShiftAssignments" \
		 $numExpanded [llength $shiftAssigns]]
    }
    
    return ""
}
 

proc cloneStereospecificShiftAssignments args {

    #
    # Given a set of shiftAssignments, 
    #
    # Find each shiftAssignment that corresponds to half a stereopair.
    #
    # If a stereopartner shiftAssignment exists, then create a new shiftAssignment
    # with the original shiftAssignment's name (appended by _stereopartner) and chemical shifts, 
    # but with the partner shiftAssignment's selections.
    #
    # If no stereopartner shiftAssignment exists, expand the original shiftAssignment's 
    # selections to cover the entire stereopair.
    #
    # Note that this includes all functionality of expandHalfStereopairShiftAssignments
    #

    set pot          [requiredFlagVal $args -pot]
    set remVar       [flagVal $args -remarksVariableName ""]

    global errorInfo

    set assignCount    0
    set numExpanded    0
    set numClonesAdded 0

    set shiftAssigns [$pot shiftAssignments]

    foreach curAssign $shiftAssigns {
	
	updateUser [format "Looking for stereopartner of shift assign %d of %d \r" \
			[incr assignCount] [llength $shiftAssigns]]
	
	if {[catch {set partner [findStereopartnerForShiftAssignment \
				     -shiftAssignmentList $shiftAssigns \
				     -targetShiftAssignment $curAssign] } msg]} {

	    if {$msg == "No possible partner"} {
		set errorInfo ""
		continue
	    } else {
		set savedInfo $errorInfo
		error $msg $savedInfo
	    }
	}

	#
	# if I didn't find a shiftAssignment with the stereopartner's selections, 
	# expand the original shiftAssignment's selections to cover both stereopartners
	#

	if {$partner == ""} {

	    ShiftAssignment -this $curAssign
	    
	    $curAssign setProtonSelectionString    [makeSelectionNonStereoSpecific [$curAssign protonSelectionString]]
	    if {[$curAssign hasHeavyatomSelection]} {
		$curAssign setHeavyatomSelectionString \
		    [makeSelectionNonStereoSpecific [$curAssign heavyatomSelectionString]]
	    }

	    $curAssign appendToNote "selections expanded by cloneStereospecificShiftAssignments"
	    rename $curAssign ""
	    incr numExpanded

	} else {

	    ShiftAssignment -this $curAssign
	    ShiftAssignment -this $partner
	    set clone [ShiftAssignment -args [format "%s_stereopartner" [$curAssign name]]]
	    $pot addShiftAssignment [$clone cget -this]

	    $clone setProtonSelectionString    [$partner protonSelectionString]

	    if {[$partner hasHeavyatomSelection]} {
		$clone setHeavyatomSelectionString [$partner heavyatomSelectionString]
	    }

	    if {[$curAssign isFrom]} {
		$clone setIsFrom
	    } else {
		$clone setIsTo
	    }
	    
	    if {[$curAssign hasProtonShift]} {
		$clone setProtonShift [$curAssign protonShift]
	    }

	    if {[$curAssign hasHeavyatomShift]} {
		$clone setHeavyatomShift [$curAssign heavyatomShift]
	    }

	    if {[$curAssign isGood]} {
		$clone setIsGood
	    } 

	    if {[$curAssign isMethyl]} {
		$clone setIsMethyl
	    }

	    if {[$curAssign hasPreviousLikelihood]} {
		$clone setPreviousLikelihood [$curAssign previousLikelihood]
	    }

	    $clone setNote [format "stereopartner clone of shiftAssign %s, selections from shiftAssign %s" \
				[$curAssign name] [$partner name]]

	    $clone -disown
	    rename $clone ""

	    incr numClonesAdded
	}
    }
    
    if {$remVar != ""} {

	set rpt "$numExpanded unpartnered stereospecific shiftAssignments were expanded by cloneStereospecificShiftAssignments"
	appendLineToString rpt "$numClonesAdded stereopartner clone shiftAssignments were added by cloneStereospecificAshiftAssignments"
	
	upvar $remVar tempRem
	lappend tempRem $rpt
    }

    return ""
}



proc markMethylShiftAssignments args {

    set shiftAssignments [requiredFlagVal $args -shiftAssignments]
    set remVar           [flagVal $args -remarksVariableName ""]

    set methylSel [AtomSel -args [methylSelectionString]]

    set count 0
    set nMethyls 0

    foreach curAssign $shiftAssignments {

	ShiftAssignment -this $curAssign

	updateUser [format "Looking for methyls in shift assignment %s (%d of %d) \r" \
			[$curAssign name] [incr count] [llength $shiftAssignments]]

	
	if {[$methylSel intersects [$curAssign protonSelection]]} {
	    incr nMethyls
	    $curAssign setIsMethyl
	} else {
	    $curAssign setNotMethyl
	}
	
	rename $curAssign ""
    }

    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "Found %d methyl shift assignments." $nMethyls]
    }
    return ""
}


proc swapShiftAssignments {a b} {

    ShiftAssignment -this $a
    ShiftAssignment -this $b

    set tempProt  [$a protonShift]
    set tempHeavy [$a heavyatomShift]

    $a setProtonShift    [$b protonShift]
    $a setHeavyatomShift [$b heavyatomShift]

    $b setProtonShift    $tempProt
    $b setHeavyatomShift $tempHeavy

    rename $a ""
    rename $b ""
}


proc determineSALikelihoodsFromConvergedStructs args {

    set files      [requiredFlagVal $args -files]
    set pot        [requiredFlagVal $args -pot]
    set violCutoff [requiredFlagVal $args -violCutoff]
    set ncCutoff   [requiredFlagVal $args -noeCompletenessCutoff]
    set psCutoff   [requiredFlagVal $args -scatterCutoff]
    set invBound   [requiredFlagVal $args -inverseBound]
    set invMeth    [requiredFlagVal $args -inverseMethylCorrection]

    #
    # init array to store times each shiftAssignment has completeness > cutoff
    # and peak position scatter < cutoff
    #
    
    foreach curSA [$pot shiftAssignments] {

	ShiftAssignment -this $curSA
	set timesOK([$curSA name]) 0
	rename $curSA ""
    }
    
    # for each converged struct,
    
    $pot setInverseBound $invBound
    $pot setInverseMethylCorrection $invMeth
    
    foreach curFile $files {
		
	xplorSim setAtomPosArr [lindex $curFile 1]

	#
	# calculate NOE completeness for every shiftAssignment, 
	# using the non-violated NOE assignments as the active set
	#

	activateNonviolatedPeakAssigns \
	    -pot $pot \
	    -violCutoff $violCutoff

	$pot updateNoeCompleteness
	$pot updatePeakPositionScatter

	foreach curSA [$pot shiftAssignments] {

	    ShiftAssignment -this $curSA

	    if {([$curSA noeCompleteness] > $ncCutoff) && ([$curSA protonPeakPositionScatter] < $psCutoff)} {	
		incr timesOK([$curSA name]) 
	    }

	    rename $curSA ""
	}
    }

    #
    # now calc prevLikelihoods as frac of converged structs in which
    # noeCompleteness is > cutoff & scatter < cutoff
    #

    set nConverged [llength $files]

    foreach curSA [$pot shiftAssignments] {
	
	ShiftAssignment -this $curSA
	$curSA setPreviousLikelihood [expr $timesOK([$curSA name]) / double($nConverged)]
	rename $curSA ""
    }

    return ""
}


proc reportShiftAssignLikelihoods args {

    set pot [requiredFlagVal $args -pot]

    set temp [list]
    foreach curSA [$pot shiftAssignments] {
	ShiftAssignment -this $curSA
	lappend temp [$curSA previousLikelihood]
	rename $curSA ""
    }

    return [histogram -data $temp -min 0 -max 1 -title "ShiftAssignment likelihoods"]
}
	       

proc shiftAssignmentsInvolvingSelectedAtoms args {

    set shiftAssignments [requiredFlagVal $args -shiftAssignments]
    set sel              [requiredFlagVal $args -sel]

    set retVal [list]

    foreach curSA $shiftAssignments {
	ShiftAssignment -this $curSA
	
	set curMatches 0

	if {[$curSA hasProtonSelection]} {
	    set tempSel [$curSA protonSelection]
	    AtomSel -this $tempSel
	    if {[$tempSel intersects [$sel cget -this]]} {
		set curMatches 1
	    }
	    rename $tempSel ""
	}

	if {[$curSA hasHeavyatomSelection]} {
	    set tempSel [$curSA heavyatomSelection]
	    AtomSel -this $tempSel
	    if {[$tempSel intersects [$sel cget -this]]} {
		set curMatches 1
	    }
	    rename $tempSel ""
	}

	if {$curMatches == 1} {
	    lappend retVal [$curSA name]
	}

	rename $curSA ""
    }

    return $retVal
}

proc shiftAssignmentsWithOverlappingSelections args {
    
    set shiftAssignments [requiredFlagVal $args -shiftAssignments]
    set target           [requiredFlagVal $args -target]

    set retVal [list]

    foreach curSA $shiftAssignments {
	ShiftAssignment -this $curSA
	
	if {[$curSA name] == [$target name]} {
	    continue
	}
	
	set curMatches 0

	if {[$curSA hasProtonSelection] && [$target hasProtonSelection]} {
	    set tempSel [$curSA protonSelection]
	    AtomSel -this $tempSel
	    if {[$tempSel intersects [$target protonSelection]]} {
		set curMatches 1
	    }
	    rename $tempSel ""
	}
	
	if {[$curSA hasHeavyatomSelection] && [$target hasHeavyatomSelection]} {
	    set tempSel [$curSA heavyatomSelection]
	    AtomSel -this $tempSel
	    if {[$tempSel intersects [$target heavyatomSelection]]} {
		set curMatches 1
	    }
	    rename $tempSel ""
	}
	
	if {$curMatches == 1} {
	    lappend retVal [$curSA name]
	}
	
	rename $curSA ""
    }
    
    return $retVal
}

proc neighboringShiftAssignments args {

    set shiftAssignments [requiredFlagVal $args -shiftAssignments]
    set nStructs         [flagVal $args -numStructures 50]
    set distCutoff       [flagVal $args -distanceCutoff 3.0]
    set fracCutoff       [flagVal $args -fractionCutoff 0.9]

    for {set x 0} {$x < $nStructs} {incr x} {

	updateUser [format "Calculating neighboring shift assignments, cycle %d of %d\r" $x $nStructs]

	if {$x == 0} {
	    makeRandomStructure -setupIVM
	} else {
	    makeRandomStructure
	}

	set curDists [recordInterShiftAssignDistances -shiftAssignments $shiftAssignments]

	if {$x == 0} {
	    foreach elem $curDists {
		set timesNear([join [list [lindex $elem 0] [lindex $elem 1]]]) 0
	    }
	} 

	foreach elem $curDists {
	    if {[lindex $elem 2] <= $distCutoff} {
		incr timesNear([join [list [lindex $elem 0] [lindex $elem 1]]]) 
	    }
	}
    }

    set retVal [list]

    foreach elem [array names timesNear] {
	set curFrac [expr $timesNear($elem) / double($nStructs)]
	if {$curFrac >= $fracCutoff} {
	    lappend retVal [split $elem]
	}
    }

    return $retVal
}

proc recordInterShiftAssignDistances args {

    set shiftAssignments [requiredFlagVal $args -shiftAssignments]

    set fromSAs [list]
    set toSAs [list]

    foreach SA $shiftAssignments {

	ShiftAssignment -this $SA
	if {[$SA isFrom]} {
	    lappend fromSAs [$SA cget -this]
	} else {
	    lappend toSAs [$SA cget -this]
	}
	rename $SA ""
    }

    set retVal [list]

    foreach fromSA $fromSAs {
	ShiftAssignment -this $fromSA
	foreach toSA $toSAs {
	    ShiftAssignment -this $toSA

	    set curDist [$fromSA distanceToShiftAssignment [$toSA cget -this]]
	    
	    lappend retVal [list [$fromSA name] [$toSA name] $curDist]

	    rename $toSA ""
	}
	rename $fromSA ""
    }

    return $retVal
}



proc initialShiftAssignmentAnalysis args {

    set pot                 [requiredFlagVal $args -pot]
    set referenceStructFile [flagVal $args -referenceStructureFile ""]
    set cutoff              [flagVal $args -completenessCutoff 4.0]
    set methylCorrect       [flagVal $args -methylCorrection 0.0]
    set violCutoff          [flagVal $args -violCutoff 0.5]
    set remVar              [flagVal $args -remarksVariableName ""]
    set longRangeCutoff     [flagVal $args -longRangeCutoff 5]
    set minLikelihoodCutoff [flagVal $args -minLikelihood -1]
    set maxLikelihoodCutoff [flagVal $args -maxLikelihood 2]
    set description         [flagVal $args -description "initialMatch"]
    set isVerbose           [flagExists $args -verbose]

    set remList [list]
    
    lappend remList [format "initial shiftAssignment analysis at stage %s, using peak assignments with likelihoods between %f and %f" $description $minLikelihoodCutoff $maxLikelihoodCutoff]

    if {$referenceStructFile != ""} {
		
 	lappend remList "reference structure is $referenceStructFile"
	
 	readPDB -fileName $referenceStructFile
	
	markGoodPeakAssignments \
	    -peakList [$pot peaks] \
	    -violCutoff $violCutoff
    }


    #
    # report on SA use
    #
    
    foreach sa [$pot shiftAssignments] {
 	ShiftAssignment -this $sa
	
	$sa setInRightNeighborhood
	
 	set LRPAs([$sa name])        [list]
 	set SRPAs([$sa name])        [list]
 	set SeqPAs([$sa name])       [list]
 	set IntraPAs([$sa name])     [list]
	
 	set GoodLRPAs([$sa name])    [list]
 	set GoodSRPAs([$sa name])    [list]
 	set GoodSeqPAs([$sa name])   [list]
 	set GoodIntraPAs([$sa name]) [list]
	
 	rename $sa ""
    }

    #
    # gather up statistics of peakAssigns used by each SA.  
    #
    # Also activiate/inactivate them based on previousLikelihood selection 
    # for use during completeness calculation below
    #
    # Also build an index of peakAssigns by ShiftAssign names for use in 
    # completeness calculation below
    #
    
    foreach curPeak [$pot peaks] {
 	Peak -this $curPeak
 	foreach curPA [$curPeak peakAssignments] {
 	    PeakAssignment -this $curPA 
	    
	    if {([$curPA previousLikelihood] > $maxLikelihoodCutoff) ||
		([$curPA previousLikelihood] < $minLikelihoodCutoff)} {

		$curPA inactivate

		rename $curPA ""
		continue
	    }

	    #
	    # activate this selected PA, and record both its SA pair and their toFromPartners 
	    # as corresponding to an active peakAssign contact
	    #

	    $curPA activate

	    set peakAssignContacts([join [list [$curPA fromAssignmentName] [$curPA toAssignmentName]]]) 1
	    
	    set fromSA [$curPA fromAssignment] 
	    set toSA   [$curPA   toAssignment]
	    ShiftAssignment -this $fromSA
	    ShiftAssignment -this $toSA
	    if {[$fromSA hasToFromPartnerName] && [$toSA hasToFromPartnerName]} {
		set peakAssignContacts([join [list [$toSA toFromPartnerName] [$fromSA toFromPartnerName]]]) 1
	    }
	    rename $fromSA ""
	    rename $toSA ""

 	    if {[$curPA isIntraresidue]} {
		
 		lappend IntraPAs([$curPA fromAssignmentName]) \
		    [list [$curPeak name] [$curPA name]]
 		lappend IntraPAs([$curPA toAssignmentName])   \
		    [list [$curPeak name] [$curPA name]]
		
 		if {[$curPA isGood]} {
 		    lappend GoodIntraPAs([$curPA fromAssignmentName]) [list [$curPeak name] [$curPA name]]
 		    lappend GoodIntraPAs([$curPA toAssignmentName])   [list [$curPeak name] [$curPA name]]
 		}
		
 	    } elseif {[$curPA isSequential]} {
		
 		lappend SeqPAs([$curPA fromAssignmentName]) [list [$curPeak name] [$curPA name]]
 		lappend SeqPAs([$curPA toAssignmentName])   [list [$curPeak name] [$curPA name]]
		
 		if {[$curPA isGood]} {
 		    lappend GoodSeqPAs([$curPA fromAssignmentName]) [list [$curPeak name] [$curPA name]]
 		    lappend GoodSeqPAs([$curPA toAssignmentName])   [list [$curPeak name] [$curPA name]]
 		}
		
 	    } elseif {[$curPA isShortRange $longRangeCutoff]} {
		
 		lappend SRPAs([$curPA fromAssignmentName]) [list [$curPeak name] [$curPA name]]
 		lappend SRPAs([$curPA toAssignmentName])   [list [$curPeak name] [$curPA name]]
		
 		if {[$curPA isGood]} {
 		    lappend GoodSRPAs([$curPA fromAssignmentName]) [list [$curPeak name] [$curPA name]]
 		    lappend GoodSRPAs([$curPA toAssignmentName])   [list [$curPeak name] [$curPA name]]
 		}
		
 	    } else {
		
 		lappend LRPAs([$curPA fromAssignmentName]) [list [$curPeak name] [$curPA name]]
 		lappend LRPAs([$curPA toAssignmentName])   [list [$curPeak name] [$curPA name]]
		
 		if {[$curPA isGood]} {
 		    lappend GoodLRPAs([$curPA fromAssignmentName]) [list [$curPeak name] [$curPA name]]
 		    lappend GoodLRPAs([$curPA toAssignmentName])   [list [$curPeak name] [$curPA name]]
 		}
 	    }
	    
 	    rename $curPA ""
 	}
 	rename $curPeak ""
    }
    
    set nIntraPAs [list]
    set nSeqPAs   [list]
    set nSRPAs    [list]
    set nLRPAs    [list]
    set nTotPAs   [list]
    
    set nGoodIntraPAs [list]
    set nGoodSeqPAs   [list]
    set nGoodSRPAs    [list]
    set nGoodLRPAs    [list]
    set nGoodTotPAs   [list]
    
    foreach curSAname [array names LRPAs] {
 	lappend nIntraPAs     [llength $IntraPAs($curSAname)]
 	lappend nSeqPAs       [llength $SeqPAs($curSAname)]
 	lappend nSRPAs        [llength $SRPAs($curSAname)]
 	lappend nLRPAs        [llength $LRPAs($curSAname)]
 	lappend nTotPAs       [expr [llength $IntraPAs($curSAname)] + [llength $SeqPAs($curSAname)] + [llength $SRPAs($curSAname)] + [llength $LRPAs($curSAname)]]
 	lappend nGoodIntraPAs [llength $GoodIntraPAs($curSAname)]
 	lappend nGoodSeqPAs   [llength $GoodSeqPAs($curSAname)]
 	lappend nGoodSRPAs    [llength $GoodSRPAs($curSAname)]
 	lappend nGoodLRPAs    [llength $GoodLRPAs($curSAname)]
 	lappend nGoodTotPAs   [expr [llength $GoodIntraPAs($curSAname)] + [llength $GoodSeqPAs($curSAname)] + [llength $GoodSRPAs($curSAname)] + [llength $GoodLRPAs($curSAname)]]
    }
    
    lappend remList [histogram -data $nTotPAs   -min 0 -integerBins -title "Number of peak assignments per shiftAssignment"]
    lappend remList [histogram -data $nIntraPAs -min 0 -integerBins -title "Number of intraresidue peak assignments per shiftAssignment"]
    lappend remList [histogram -data $nSeqPAs   -min 0 -integerBins -title "Number of sequential peak assignments per shiftAssignment"]
    lappend remList [histogram -data $nSRPAs    -min 0 -integerBins -title "Number of short range peak assignments per shiftAssignment"]
    lappend remList [histogram -data $nLRPAs    -min 0 -integerBins -title "Number of long range peak assignments per shiftAssignment"]
    
    if {[min $nTotPAs] == 0} {
 	set newRem "Shift assignments that make no peak assignments:"
 	foreach curSAname [array names LRPAs] {
 	    set curTot [expr [llength $IntraPAs($curSAname)] + [llength $SeqPAs($curSAname)] + [llength $SRPAs($curSAname)] + [llength $LRPAs($curSAname)]]
 	    if {$curTot == 0} {
 		appendLineToString newRem [format "   %s" $curSAname]
 	    }
 	}
 	lappend remList $newRem
    }
    
    
    if {[max $nGoodTotPAs] > 0} {
 	lappend remList [histogram -data $nGoodTotPAs   -min 0 -integerBins -title "Number of good peak assignments per shiftAssignment"]
 	lappend remList [histogram -data $nGoodIntraPAs -min 0 -integerBins -title "Number of good intraresidue peak assignments per shiftAssignment"]
 	lappend remList [histogram -data $nGoodSeqPAs   -min 0 -integerBins -title "Number of good sequential peak assignments per shiftAssignment"]
 	lappend remList [histogram -data $nGoodSRPAs    -min 0 -integerBins -title "Number of good short range peak assignments per shiftAssignment"]
 	lappend remList [histogram -data $nGoodLRPAs    -min 0 -integerBins -title "Number of good long range peak assignments per shiftAssignment"]
	
 	if {[min $nGoodTotPAs] == 0} {
 	    set newRem "Shift assignments that make no good peak assignments:"
 	    foreach curSAname [array names LRPAs] {
 		set curGoodTot [expr [llength $GoodIntraPAs($curSAname)] + [llength $GoodSeqPAs($curSAname)] + [llength $GoodSRPAs($curSAname)] + [llength $GoodLRPAs($curSAname)]]
 		if {$curGoodTot == 0} {
 		    appendLineToString newRem [format "   %s" $curSAname]
 		}
 	    }
 	    lappend remList $newRem
 	}
    }
    
    if {$referenceStructFile != ""} {
	
 	#
 	# Report on completeness
 	#
		
 	$pot setInverseBound $cutoff
 	$pot setInverseMethylCorrection $methylCorrect
	
 	$pot updateNoeCompleteness
	
 	set NCs [list]
 	foreach sa [$pot shiftAssignments] {
 	    ShiftAssignment -this $sa
 	    lappend NCs [expr 100.0 * [$sa noeCompleteness]]
 	    $sa appendToNote [format "ref struct completeness %f%%" \
				  [expr 100.0 * [$sa noeCompleteness]]]
 	    rename $sa ""
 	}
	
 	lappend remList [histogram -data $NCs -title "Distribution of NOE completeness in all shiftAssignments"]
 	lappend remList [format "average completeness %f%% +/- %f%% \noverall completeness %f%%" [mean $NCs] [standardDeviation $NCs] [expr 100.0 * [$pot noeCompletenessScore]]]
	
	
	
 	#
 	# report on explicit exception accuracy
 	#
	
	
 	set fromSAs [$pot fromShiftAssignments]
 	set toSAs   [$pot toShiftAssignments]
	
	#
	# record every explicit exception, and every exception's toFromPartners, as exceptions
	#

 	foreach {fromname toname} [$pot explicitInverseExceptions] {
 	    set exceptions([join [list $fromname $toname]]) 1

  	    set fromSA [$pot shiftAssignmentNamed $fromname]
  	    set   toSA [$pot shiftAssignmentNamed $toname]
  	    ShiftAssignment -this $fromSA
  	    ShiftAssignment -this $toSA
	    
  	    if {[$fromSA hasToFromPartnerName] && [$toSA hasToFromPartnerName]} {
  		set exceptions([join [list [$toSA toFromPartnerName] [$fromSA toFromPartnerName]]]) 1
  	    }
  	    rename $fromSA ""
  	    rename $toSA ""
 	}
	

	set nExceptions [llength [array names exceptions]]
	set nPeakAssignContacts [llength [array names peakAssignContacts]]
	
 	set nCloseContacts 0
 	set nGoodExceptions 0
	set nGoodPeakAssignContacts 0
	set nMissingExceptions 0
	set nNonPAExceptions 0
	set nGoodNonPAExceptions 0
	
 	foreach fromSA $fromSAs {
 	    ShiftAssignment -this $fromSA
 	    foreach toSA $toSAs {
 		ShiftAssignment -this $toSA

		set curPairHasException [info exists exceptions([join [list [$fromSA name] [$toSA name]]])]
		set curPairHasPA        [info exists peakAssignContacts([join [list [$fromSA name] [$toSA name]]])]
		
		if {$curPairHasException && ! $curPairHasPA} {
		    incr nNonPAExceptions
		}


 		set curDist [$fromSA distanceToShiftAssignment [$toSA cget -this]]
		
 		if {$curDist < $cutoff} {

 		    incr nCloseContacts
		    
 		    if {$curPairHasException} {
 			incr nGoodExceptions
 		    }

		    if {$curPairHasPA} {
			incr nGoodPeakAssignContacts
		    }

		    if {$curPairHasException && ! $curPairHasPA} {
			incr nGoodNonPAExceptions
		    }

		    if {! ($curPairHasException || $curPairHasPA)} {
			incr nMissingExceptions
		    } 
 		}
		
 		rename $toSA ""
 	    }
 	    rename $fromSA ""
 	}
	
 	set nBadExceptions [expr $nExceptions - $nGoodExceptions]
	set nBadPeakAssignContacts [expr $nPeakAssignContacts - $nGoodPeakAssignContacts]
 	set nPossibleContacts [expr [llength $fromSAs] * [llength $toSAs]]

	set pctGoodExceptions 0
	if {$nExceptions != 0} {
	    set pctGoodExceptions [expr ($nGoodExceptions / double($nExceptions)) * 100]
	}

	set pctGoodPeakAssignContacts 0
	if {$nPeakAssignContacts != 0} {
	    set pctGoodPeakAssignContacts [expr ($nGoodPeakAssignContacts / double($nPeakAssignContacts)) * 100]
	}

	set pctGoodNonPAExceptions 0
	if {$nNonPAExceptions != 0} {
	    set pctGoodNonPAExceptions [expr ($nGoodNonPAExceptions / double($nNonPAExceptions)) * 100]
	}


	set pctMissingExceptions 0
	if {$nExceptions != 0} {
	    set pctMissingExceptions [expr ($nMissingExceptions / double($nCloseContacts)) * 100]
	}

 	set newRem [format "Inverse exception accuracy:\n"]

 	appendLineToString newRem [format "Of %d shiftAssignment pairs corresponding to peakAssignments currently selected, %d (%f%%) are good--ie., correspond to contacts closer than %f A in the reference structure" \
 				       $nPeakAssignContacts $nGoodPeakAssignContacts $pctGoodPeakAssignContacts $cutoff ]

 	appendLineToString newRem [format "Of %d explicit inverse exceptions, %d (%f%%) are good" \
 				       $nExceptions $nGoodExceptions $pctGoodExceptions ]

	appendLineToString newRem [format "Of %d explicit inverse exceptions NOT corresponding to currently selected peakAssignments, %d (%f%%) are good" \
 				       $nNonPAExceptions $nGoodNonPAExceptions $pctGoodNonPAExceptions ]


 	appendLineToString newRem [format "%d pairs of shiftAssignments (%f%% of the total) made close contacts" $nCloseContacts [expr ($nCloseContacts / double($nPossibleContacts)) * 100]]

 	appendLineToString newRem [format "Of them, %d close contacts (%f%% of all close contacts) were missing explicit exceptions or currently selected peakAssignments" \
 				       $nMissingExceptions $pctMissingExceptions]
	
 	lappend remList $newRem
    }
    
    if {$remVar != ""} {
 	upvar $remVar r
 	foreach elem $remList {
 	    lappend r $elem
 	}
    }
    
    return ""
}


proc covalentNeighborhood args {

    set shiftAssignList [requiredFlagVal $args -shiftAssignList]
    set target          [requiredFlagVal $args -target]
    set bondArrayName   [requiredFlagVal $args -bondArrayName]

    upvar $bondArrayName bondArray

    #
    # Find SAs in the given list whose heavyatom selection is covalently bonded to 
    # the target's heavyatom selection.
    #
    # Missing heavyatom selections are generated automatically from the proton selections
    #

    if {[$target hasHeavyatomSelection]} {
	set targetHeavySel [$target heavyatomSelection]
	AtomSel -this $targetHeavySel
    } else {
	set targetProtonSel [$target protonSelection]
	AtomSel -this $targetProtonSel
	set targetHeavySel [AtomSel -args [format "(bondedto (%s))" [$targetProtonSel string]]]
	rename $targetProtonSel ""
    }

    set retVal [list]
    
    foreach curSA $shiftAssignList {
	ShiftAssignment -this $curSA

	if {[$curSA name] == [$target name]} {
	    lappend retVal [$curSA name]
	    continue
	}

	if {[$curSA hasHeavyatomSelection]} {
	    set curHeavySel [$curSA heavyatomSelection]
	    AtomSel -this $curHeavySel
	} else {
	    set curProtonSel [$curSA protonSelection]
	    AtomSel -this $curProtonSel
	    set curHeavySel [AtomSel -args [format "(bondedto (%s))" [$curProtonSel string]]]
	    rename $curProtonSel ""
	}

	
	if {[selectionsAreBonded $targetHeavySel $curHeavySel bondArray] || 
	    [$targetHeavySel intersects [$curHeavySel cget -this]]} {

	    lappend retVal [$curSA name]
	}
	
	rename $curHeavySel ""
	rename $curSA ""
    }
    
    return $retVal
}


proc groupShiftAssignmentsByResidue args {

    set shiftAssignList [requiredFlagVal $args -shiftAssignList]

    foreach curSA $shiftAssignList {

	ShiftAssignment -this $curSA
	set curSel [$curSA protonSelection]
	AtomSel -this $curSel

	set curResidues [residuesInSelection $curSel]
	
	foreach curRes $curResidues {
	    lappend SAsInRes($curRes) [$curSA name]
	}

	rename $curSel ""
	rename $curSA ""
    }

    set retVal [list]

    foreach curRes [array names SAsInRes] {
	
	lappend retVal $SAsInRes($curRes)
    }

    return $retVal
}



proc findEquivShiftAssignments {potA potB} {

    #
    # Given two pots, find any shiftAssigns in the second pot that 
    # are equivaelent to the SAs in the first pot and vice versa.
    #
    # Equivalence is based on proton selections.  If there are >1 
    # equivalent SAs found (eg., non-stereo-processed pairs), the 
    # single best match is chosen on the basis of minimal proton 
    # chemical shift difference.
    #
    # Note that I can't require exact proton shift matches, because 
    # of stripe based correction in each spectrum
    #

    #
    # record atoms selected by proton selection for each SA, as a string
    # of indexes for fast lookup
    #

    foreach curPot [list $potA $potB] {
	foreach curSA [$curPot shiftAssignments] {
	    ShiftAssignment -this $curSA 
	
	    set curSel [$curSA protonSelection]
	    set curIndexCode [join [$curSel indices]]
	    set curData [list [$curPot instanceName] [$curSA name] [$curSA protonShift]]

	    set indexesForSA([join [list [$curPot instanceName] [$curSA name]]]) $curIndexCode
	    lappend SAsByIndexes($curIndexCode) $curData

	    rename $curSA ""
	    rename $curSel ""
	}
    }

    set retVal [list]
    set noEquivs [list]
    set oneEquivs [list]
    set multiEquivs [list]

    foreach curPot [list $potA $potB] otherPot [list $potB $potA] {

	foreach curSA [$curPot shiftAssignments] {
	    ShiftAssignment -this $curSA 
	
	    set curIndexCode $indexesForSA([join [list [$curPot instanceName] [$curSA name]]])
	
	    set curEquivs [list]
	    set minDeltaShift 9999999
	    
	    foreach equivSA $SAsByIndexes($curIndexCode) {
		
		set equivPot    [lindex $equivSA 0]
		set equivSAname [lindex $equivSA 1]
		set equivShift  [lindex $equivSA 2] 
		
		if {$equivPot == [$curPot instanceName]} {
		    continue
		} 
		
		set curDeltaShift [expr abs($equivShift - [$curSA protonShift])]

		if {$curDeltaShift < $minDeltaShift} {
		    set curEquivs [list]
		    set minDeltaShift $curDeltaShift
		}
		
		if {$curDeltaShift == $minDeltaShift} {
		    lappend curEquivs [list [$otherPot instanceName] $equivSAname]
		}
	    }
	    
	    switch [llength $curEquivs] {
		
		0 { 
		    lappend retVal [list [$curPot instanceName] [$curSA name] "NONE"] 
		    lappend noEquivs [format "%s %s" [$curPot instanceName] [$curSA name]]
		}
		
		1 { 
		    lappend retVal [list [$curPot instanceName] [$curSA name] [lindex $curEquivs 0]] 
		    lappend oneEquivs [format "%s %s --> %s %s" [$curPot instanceName] [$curSA name] [lindex [lindex $curEquivs 0] 0] [lindex [lindex $curEquivs 0] 1]]
		} 
		
		default { 
		    lappend retVal [list [$curPot instanceName] [$curSA name] [lindex $curEquivs 0]] 
		    lappend multiEquivs [format "%s %s  -->  %s" [$curPot instanceName] [$curSA name] $curEquivs]
		}
	    }

	    rename $curSA ""
	}
    }


	
    puts "No equivalent SAs:"
    foreach elem $noEquivs {
	puts $elem
    }

    puts "\n\nOne equivalent SA:" 
    foreach elem $oneEquivs {
	puts $elem
    }

    puts "\n\n multiple equivalent SAs:"
    foreach elem $multiEquivs {
	puts $elem
    }

    return $retVal
}



proc dropShiftAssignsInShiftRange args {

    set pot    [requiredFlagVal $args -pot]
    set dim    [requiredFlagVal $args -dimension]
    set minVal [requiredFlagVal $args -minVal]
    set maxVal [requiredFlagVal $args -maxVal]
    set remVar [flagVal $args -remarksVariableName ""]


    set legalDimensions [list "fromProton" "fromHeavyatom" "toProton" "toHeavyatom"]
    if {[lsearch $legalDimensions $dim] == -1} {
	error [format "Dimension type %s not in %s" $dim $legalDimensions]
    }

    if {$dim == "fromProton"} {
	set tf "isFrom"
	set cmd "protonShift"
    } elseif {$dim == "fromHeavyatom"} {
	set tf "isFrom"
	set cmd "heavyatomShift"
    } elseif {$dim == "toProton"} {
	set tf "isTo"
	set cmd "protonShift"
    } elseif {$dim == "toHeavyatom"} {
	set tf "isTo"
	set cmd "heavyatomShift"
    }

    set droppedSAs [list]

    foreach curSA [$pot shiftAssignments] {
	ShiftAssignment -this $curSA

	if {[$curSA $tf]} {
	    set curShift [$curSA $cmd] 
	    if {($curShift >= $minVal) && ($curShift <= $maxVal)} {
		$pot removeShiftAssignmentNamed [$curSA name]
		$curSA -acquire
		lappend droppedSAs [$curSA name]
	    }
	}

	rename $curSA ""
    }

    if {$remVar != ""} {
	upvar $remVar tempRem

	set newRem [format "DropShiftAssignmentsInShiftRange: Dropped %d shiftAssigns with %s in range %f .. %f ppm" \
			[llength $droppedSAs] $dim $minVal $maxVal]

	lappend tempRem $newRem
    }

    return ""
}

