#
# procedures to read and write NOE peaks,
# along with procs to choose good/bad peaks and peakAssigns
#
# JJK 4/24/02
#

package provide marvin 1.0
package require pasdpot




proc readMarvinPeaks args {

    global errorInfo
    
    set fileName   [requiredFlagVal $args -fileName]
    set potential  [requiredFlagVal $args -pot]
    
    if { [catch {set inUnit [open $fileName r]}] } {
    	error "Error opening input file $fileName"
    }
    
    XplorCommand "set message off echo off end"
    
    while {![eof $inUnit]} {
	
	if { [catch {set curNOE [readOneMarvinPeak $inUnit]} msg]} {

	    if {$msg == "No more peaks"} {
		set errorInfo ""
		break
	    } else {
		error [format "Peak reading error in file %s:\n %s" $fileName $msg]
	    }
	}

	processOnePeak $curNOE $potential
	updateUser [format "read peak %s \r" [lindex $curNOE 0]]
    }

    $potential updateLinkedPeakAssignments
    $potential updatePrimarySeqDists 

    close $inUnit
}


proc readOneMarvinPeak {fileID} {

    # 
    # Try to find the beginning of another peak,
    # ignoring text in xplor-style comments 
    #
    
    if {[catch {eatXplorTextUpTo $fileID [list "peak" "restraint"]}]} {

	error "No more peaks"
    }

    # defaults for optional fields 
    
    set curNote ""
    set curUpBound 0.0
    set curLowBound 0.0
    set intensity "NONE"
    set fps "NONE"
    set tps "NONE"
    set fhs "NONE"
    set ths "NONE"
    set peakAssignList {}

    # grab the required fields

    set curName [nextWord $fileID]
    if {$curName == ""} {
	error "Missing peak name"
    }
    
    # search for optional fields

    set done 0
    while {! $done} {
	
	set curWord [nextWord $fileID]
	
	switch -exact -- $curWord {
	    
	    "" {
		error [format "unexpected end of file at peak %s" $curName]
	    }
	    
	    "end" {
		set done 1
	    }
	    
	    "peak" -
	    "restraint" {
		error [format "missing end statement after peak %s" $curName]
	    }
	    
	    "-bounds" {
		
		set boundA [nextWord $fileID]
		set boundB [nextWord $fileID]
		
		if {$boundA < $boundB} {
		    
		    set curLowBound $boundA
		    set curUpBound $boundB
		    
		} else {
		    
		    set curUpBound $boundA
		    set curLowBound $boundB
		    
		}

		if {! [isPosNumber $curUpBound]} {
		    error [format "Missing or invalid upBound for peak %s" $curName]
		}

		if {! [isPosNumber $curLowBound]} {
		    error [format "Missing of invalid lowBound for peak %s" $curName]
		}
	    }

	    "-intensity" {

		set intensity [nextWord $fileID]

		if {! [isNumber $intensity]} {
		    error [format "Missing or invalid intensity for peak %s" $curName]
		}
	    }


	    "-fromProtonShift" {

		set fps [nextWord $fileID]

		if {! [isNumber $fps]} {
		    error [format "Missing or invalid from proton shift for peak %s" $curName]
		}
	    }
	   
	    "-toProtonShift" {

		set tps [nextWord $fileID]

		if {! [isNumber $tps]} {
		    error [format "Missing or invalid to proton shift for peak %s" $curName]
		}
	    }
 
	    "-fromHeavyatomShift" {

		set fhs [nextWord $fileID]

		if {! [isNumber $fhs]} {
		    error [format "Missing or invalid from heavyatom shift for peak %s" $curName]
		}
	    }
	   
	    "-toHeavyatomShift" {

		set ths [nextWord $fileID]

		if {! [isNumber $ths]} {
		    error [format "Missing or invalid to heavyatom shift for peak %s" $curName]
		}
	    }
	   
	    "-note" {
		
		if {$curNote == ""} {
		    set curNote [nextLine $fileID]
		} else {
		    set curNote [format "%s\n%s" $curNote [nextLine $fileID]]
		}
	    }

	    "peakAssign" -
	    "assign" {
		
 		if {[catch {set curPeakAssign [readOnePeakAssign $fileID]} msg]} {
		    error [format "Peak Assignment reading error at peak %s:\n %s" $curName $msg]
		}
		
		lappend peakAssignList $curPeakAssign
	    }
	}
    }
    
    return [list $curName $curUpBound $curLowBound $intensity $fps $tps \
		$fhs $ths $curNote $peakAssignList]
}


proc readOnePeakAssign {fileID} {

    # when this is called, readOneMarvinPeak has already eaten the 
    # leading "assign" or "peakAssign"

    # defaults for optional fields

    set curLikelihood "NONE"
    set curNote ""
    set curLowCorrect 0
    set curUpCorrect 0
    set curGood 0
    set curUnfoldedFromProtonPP "NONE"
    set curUnfoldedFromHeavyPP "NONE"
    set curUnfoldedToProtonPP "NONE"
    set curUnfoldedToHeavyPP "NONE"
    set aveExp "NONE"
    set curNumFiltersFailed 0
    set linkedPeakAssignList [list]

    # grab the required fields

    set curName    [nextWord $fileID]
    if {$curName == ""} {
	error "Missing peak assignment name"
    }

    set curFromAssName [nextWord $fileID]
    if {$curFromAssName == ""} {
	error [format "Missing or invalid from shift assignment name for peak assignment %s" $curName]
    }

    set curToAssName   [nextWord $fileID]
    if {$curToAssName == ""} {
	error [format "Missing or invalid to shift assignment name for peak assignment %s" $curName]
    }

    # search for optional fields

    set done 0	

    while {! $done} {
	
	set curWord [nextWord $fileID]

	switch -exact -- $curWord {

	    "" {
		error [format "unexpected end of file at peak assignment %s" $curName] 
	    }

	    "end" { 
		set done 1 
	    }

	    "restraint" -
	    "peak" -
	    "assign" -
	    "peakAssign" {
		error [format "missing end statement after peak assignment %s" $curName]
	    }

	    "-likelihood" {
	    
		set curLikelihood [nextWord $fileID]

		if {![isZeroToOneNumber $curLikelihood]} {
		    error [format "Missing or invalid likelihood for peak assignment %s" $curName]
		}
	    }
	    
	    "-upBoundCorrection" {
	    
		set curUpCorrect [nextWord $fileID]
	    
		if {![isNumber $curUpCorrect]} {
		    error [format "Missing or invalid upBoundCorrection for peak assignment %s" $curName]
		}
	    }

	    "-lowBoundCorrection" {
	    
		set curLowCorrect [nextWord $fileID]
	    
		if {![isNumber $curLowCorrect]} {
		    error [format "Missing or invalid lowBoundCorrection for peak assignment %s" $curName]
		}
	    }

	    "-unfoldedFromProtonPeakPosition" {
	    
		set curUnfoldedFromProtonPP [nextWord $fileID]
	    
		if {![isNumber $curUnfoldedFromProtonPP]} {
		    error [format "Missing or invalid unfoldedFromProtonPeakPosition for peak assignment %s" $curName]
		}
	    }

	    "-unfoldedFromHeavyatomPeakPosition" {
	    
		set curUnfoldedFromHeavyPP [nextWord $fileID]
	    
		if {![isNumber $curUnfoldedFromHeavyPP]} {
		    error [format "Missing or invalid unfoldedFromHeavyatomPeakPosition for peak assignment %s" $curName]
		}
	    }

	    "-unfoldedToProtonPeakPosition" {
	    
		set curUnfoldedToProtonPP [nextWord $fileID]
	    
		if {![isNumber $curUnfoldedToProtonPP]} {
		    error [format "Missing or invalid unfoldedToProtonPeakPosition for peak assignment %s" $curName]
		}
	    }

	    "-unfoldedToHeavyatomPeakPosition" {
	    
		set curUnfoldedToHeavyPP [nextWord $fileID]
	    
		if {![isNumber $curUnfoldedToHeavyPP]} {
		    error [format "Missing or invalid unfoldedToHeavyatomPeakPosition for peak assignment %s" $curName]
		}
	    }


	    "-numFiltersFailed" {
	    
		set curNumFiltersFailed [nextWord $fileID]
	    
		if {![isInteger $curNumFiltersFailed]} {
		    error [format "Missing or invalid numFiltersFailed for peak assignment %s" $curName]
		}

		if {$curNumFiltersFailed < 0} {
		    error [format "NumFiltersFailed is < 0 for peak assignment %s" $curName]
		}
	    }


	    "-aveExp" {

		set aveExp [nextWord $fileID]
	    }

	    "-good" {

		set curGood 1
	    }

	    "-linkedTo" {

		set curLinkedPeakName    [nextWord $fileID]
		set curLinkedPeakAssName [nextWord $fileID]

		lappend linkedPeakAssignList $curLinkedPeakName
		lappend linkedPeakAssignList $curLinkedPeakAssName
	    }

	    "-note" {

		if {$curNote == ""} {
		    set curNote [nextLine $fileID]
		} else {
		    set curNote [format "%s\n%s" $curNote [nextLine $fileID]]
		}
	    }
	}
    }
	    
    return [list $curName $curFromAssName $curToAssName $curUpCorrect \
		$curLowCorrect $curLikelihood $curGood $curNote \
		$curUnfoldedFromProtonPP $curUnfoldedFromHeavyPP \
		$curUnfoldedToProtonPP $curUnfoldedToHeavyPP \
		$curNumFiltersFailed $linkedPeakAssignList $aveExp]
}




proc readXplorNOEs args {

    global errorInfo 

    set fileName              [requiredFlagVal $args -fileName]
    set pot                   [requiredFlagVal $args -pot]
    set peakNamePrefix        [flagVal $args -peakNamePrefix "xplor"]
    set shiftAssignNamePrefix [flagVal $args -shiftAssignNamePrefix "xplor"]
    set minLowBound           [flagVal $args -minLowBound 1.8]

    if {[catch {set inUnit [open $fileName r]}]} {
    
	error "Error opening input file $fileName"
    }
    
    XplorCommand "set message off echo off end"

    #
    # read the xplor-format restraints
    #

    set shiftAssignCount 0
    set peakCount 0

    while {![eof $inUnit]} {
	
	set peakName [format "%s_%d" $peakNamePrefix [incr peakCount]]
	set note [format "Created by readXplorNOEs, file %s" $fileName]

	if {[catch {set curNOE [readOneXplorNOE $inUnit $peakName $minLowBound $note]} msg]} {

	    if {$msg == "No more peaks"} {
		set errorInfo ""
		break
	    } else {
		error [format "Error reading restraint %s:\n %s" $peakName $msg]
	    }
	} 

	#
	# convert the raw selections in the xplor NOE into 
	# shiftAssignments, record them, and convert the 
	# current NOE data so it uses them
	#

	set note [format "Created by readXplorNOEs, file %s, restraint %s" $fileName $peakName]

	set newNOE [convertSelsToShiftAssigns $curNOE $pot $shiftAssignNamePrefix shiftAssignCount $note]

	processOnePeak $newNOE $pot
    }
    
    close $inUnit
}
    

proc readOneXplorNOE {fileID peakName minLowBound note} {

    #
    # Note that this routine returns a structure containing selection strings for each PeakAssignment, 
    # rather than ShiftAssignment names.  
    #
    # ShiftAssignments are generated from these selections later on, and these selection strings
    # are replaced by shiftAssignment names at that point
    #

    #
    # FIX:  This routine doesn't parse Michael Nilges's extended format, ie., ASSIgn (sel)(sel) f f f OR (sel)(sel) OR...
    #

    global errorInfo

    #
    # eat leading text
    #

    if {[catch {eatXplorTextUpTo $fileID [list "assign" "assi" "ASSI" "ASSIGN" "ASSIgn"]}]} {

	error "No more peaks"
    }
	
    #
    # the format is sel sel real real real
    #
    
    # this needs error checking ala the Marvin version

    set sel1   [nextSel $fileID]
    set sel2   [nextSel $fileID]
    set d      [nextXplorWord $fileID]
    set dminus [nextXplorWord $fileID]
    set dplus  [nextXplorWord $fileID]
	
    #
    # translate into up and low bounds
    #
	
    set upBound  [expr $d + $dplus]
    set lowBound [expr $d - $dminus]
	
    if {$lowBound < $minLowBound} {
	set lowBound $minLowBound
    }

    #
    # try to break up the selections based on top-level ORs
    #

    if {[catch {set sel1List [splitSelection -string $sel1]}]} {

	set errorInfo ""
	set sel1List [list $sel1]
    }

    if {[catch {set sel2List [splitSelection -string $sel2]}]} {

	set errorInfo ""
	set sel2List [list $sel2]
    }
	
    set paList [list]
    set paCount -1

    foreach s1 $sel1List {
	foreach s2 $sel2List {
	    set curPAName [format "%s_%d" $peakName [incr paCount]]
	    set curPA [list $curPAName $s1 $s2 0.0 0.0 "NONE" 0 "" "NONE" "NONE" "NONE" "NONE" 0 [list]]
	    lappend paList $curPA
	}
    }

    return [list $peakName $upBound $lowBound "NONE" "NONE" "NONE" "NONE" "NONE" $note $paList]
}


#
# Xplor-format NOE restraints don't have explicit shift assignments, so
# this convert the raw selections in the xplor NOE into 
# shiftAssignments, records them, and converts the 
# NOE data so it uses them
#

proc convertSelsToShiftAssigns {peakData potential prefix countVarName note} {

    upvar $countVarName shiftAssignCount

    #
    # grab the current peak's peakAssignment list from the 
    # output of readOneXplorNOE
    #

    set paList        [lindex $peakData end]

    set newPAList [list]
    
    foreach curPA $paList {
	
	set fromSelString        [lindex $curPA 1]
	set toSelString          [lindex $curPA 2]
	
	set fromShiftAssignName "NONE"
	set toShiftAssignName   "NONE"

	#
	# search for existing shiftAssigns whose proton 
	# selections and from/to flags match the from and to
	# selections read from the xplor table
	#
	
	foreach shiftAss [$potential shiftAssignments] {
	    
	    ShiftAssignment -this $shiftAss
	    if {[$shiftAss isFrom]} {
		
		if {[$shiftAss protonSelectionString] == $fromSelString} {
		    set fromShiftAssignName [$shiftAss name]
		}
	    } else {
		if {[$shiftAss protonSelectionString] == $toSelString} {
		    set toShiftAssignName [$shiftAss name]
		}
	    }
	    rename $shiftAss ""
	}
	
	# 
	# create new shiftAssigns if needed
	#
	
	if {$fromShiftAssignName == "NONE"} {
	    
	    set fromShiftAssignName [format "%s_%d" $prefix [incr shiftAssignCount]]

	    set newShiftAssData [list $fromShiftAssignName "NONE" "NONE" $fromSelString  "NONE" \
				     "NONE" 0 0 "from" "NONE" "NONE" $note [list] "NONE" "NONE"] 
	    
	    processOneShiftAssign $newShiftAssData $potential
	}
	
	if {$toShiftAssignName == "NONE"} {
	    
	    set toShiftAssignName [format "%s_%d" $prefix [incr shiftAssignCount]]
	    
	    set newShiftAssData [list $toShiftAssignName "NONE" "NONE" $toSelString  "NONE" \
				     "NONE" 0 0 "to" "NONE" "NONE" $note [list] "NONE" "NONE"] 
	    
	    processOneShiftAssign $newShiftAssData $potential
	}


	#
	# replace the selection strings with the shiftAssign names in the peakAssignment's data
	#

	set newPAData [lreplace $curPA 1 2 $fromShiftAssignName $toShiftAssignName]

	lappend newPAList $newPAData
    }

    #
    # re-record the peak's data with the new peakAssignment data
    #

    set newPeakData [lreplace $peakData end end $newPAList]

    return $newPeakData
}


proc processOnePeak {peakData anNOEpotential} {
    
    #
    # takes the data read by readOneMarvinPeak or readOneXplorNOE
    # and actually creates a Marvin peak
    # and attaches it to the given MarvinNOEPotential
    #

    set curPeakName    [lindex $peakData 0]
    set curUpBound     [lindex $peakData 1]
    set curLowBound    [lindex $peakData 2]
    set intensity      [lindex $peakData 3]
    set fps            [lindex $peakData 4]
    set tps            [lindex $peakData 5]
    set fhs            [lindex $peakData 6]
    set ths            [lindex $peakData 7]
    set curNote        [lindex $peakData 8]
    set peakAssignList [lindex $peakData 9]
    
    #
    # Create the peak and add it to the pot
    #
    
    set aPeak [Peak -args $curPeakName]
    
    $aPeak setUpBound            $curUpBound
    $aPeak setLowBound           $curLowBound
    $aPeak setNote               $curNote

    if {$intensity != "NONE"} {
	$aPeak setIntensity          $intensity
    }

    if {$fps != "NONE"} {
	$aPeak setFromProtonShift    $fps
    }
    
    if {$tps != "NONE"} {
	$aPeak setToProtonShift      $tps
    }

    if {$fhs != "NONE"} {
	$aPeak setFromHeavyatomShift $fhs
    }

    if {$ths != "NONE"} {
	$aPeak setToHeavyatomShift   $ths
    }
    
    #
    # create the peakAssignments and add them to the peak
    #
    
    marvinPyth command {import pasd}
    set nMono [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.nMono} \
			   {} {out}] \
		      0] 1]
    set defAveExp [lindex \
		       [lindex \
			    [marvinPyth command {out=pasd.aveExp} \
				 {} {out}] \
			    0] 1]
    
    foreach curPA $peakAssignList {
	
	set curPAName              [lindex $curPA 0]
	set fromShiftAssignName    [lindex $curPA 1]
	set toShiftAssignName      [lindex $curPA 2]
	set upCorr                 [lindex $curPA 3]
	set lowCorr                [lindex $curPA 4]
	set likelihood             [lindex $curPA 5]
	set curGood                [lindex $curPA 6]
	set curPeakAssignNote      [lindex $curPA 7]
	set unfoldedFromProtonPP   [lindex $curPA 8]
	set unfoldedFromHeavyPP    [lindex $curPA 9]
	set unfoldedToProtonPP     [lindex $curPA 10]
	set unfoldedToHeavyPP      [lindex $curPA 11]
	set numFiltersFailed       [lindex $curPA 12]
	set linkedPeakAssignList   [lindex $curPA 13]
	set aveExp                 [lindex $curPA 14]
	
	#
	# create the peakAssignment and add it to the peak
	#

	set aPA [PeakAssignment -args $curPAName]
	$aPA setNMono $nMono

	if {$aveExp == "NONE"} {
	    set aveExp $defAveExp
	}

	$aPA setAveExp $aveExp
	
	set fromAss [$anNOEpotential shiftAssignmentNamed $fromShiftAssignName]
	set toAss   [$anNOEpotential shiftAssignmentNamed $toShiftAssignName]
 
	ShiftAssignment -this $fromAss
	ShiftAssignment -this $toAss

	if {! [$fromAss isFrom]} {
	    error [format "shiftAssign %s is used as from in peak assign %s but is flagged as to" \
		       $fromShiftAssignName $curPAName]
	}

	if {! [$toAss isTo]} {
	    error [format "shiftAssign %s is used as to in peak assign %s but is flagged as from" \
		       $toShiftAssignName $curPAName]
	}

	$aPA setFromAssignment      [$fromAss cget -this]
	$aPA setToAssignment        [$toAss cget -this]
	$aPA setUpBoundCorrection   $upCorr
	$aPA setLowBoundCorrection  $lowCorr
	$aPA setNote                $curPeakAssignNote
	
	if {$likelihood != "NONE"} {
	    $aPA setPreviousLikelihood  $likelihood
	}

	if {$curGood == 1} {
	    $aPA setGood
	} else {
	    $aPA setBad
	}

	if {$unfoldedFromProtonPP != "NONE"} {
	    $aPA setUnfoldedFromProtonPeakPosition $unfoldedFromProtonPP
	}

	if {$unfoldedFromHeavyPP != "NONE"} {
	    $aPA setUnfoldedFromHeavyatomPeakPosition $unfoldedFromHeavyPP
	}

	if {$unfoldedToProtonPP != "NONE"} {
	    $aPA setUnfoldedToProtonPeakPosition $unfoldedToProtonPP
	}

	if {$unfoldedToHeavyPP != "NONE"} {
	    $aPA setUnfoldedToHeavyatomPeakPosition $unfoldedToHeavyPP
	}

	for {set count 0} {$count < $numFiltersFailed} {incr count} {
	    $aPA incrementNumFiltersFailed
	}
	
	foreach {pName paName} $linkedPeakAssignList {
	    $aPA addLinkedPeakAssignmentName $pName $paName
	}

	#
	# add the peak assignment to the peak at the very end, so that
	# we're sure that it's valid
	#

	$aPeak addPeakAssignment [$aPA cget -this]

	# clean up TCL objects without deleting them on the C++ side

	rename $fromAss ""
	rename $toAss ""
	$aPA -disown
	rename $aPA ""
    }

    #
    # add the peak to the potential at the very end, so that 
    # we're sure it's valid
    #

    $anNOEpotential addPeak $aPeak
    
    $aPeak -disown
    rename $aPeak ""
}


proc writeMarvinPeaks args {

    set fileName        [requiredFlagVal $args -fileName]
    set remarksList     [flagVal $args -remarks]
    set anNOEpotential  [flagVal $args -pot]
    set anNOElist       [flagVal $args -peakList]

    if {($anNOEpotential == "") && ($anNOElist == "")} {
	error "Neither -pot or -peakList defined in call to writeMarvinPeaks"
    }
    
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
    # grab peak list from pot if it's defined
    #

    if {$anNOEpotential != ""} {
	set anNOElist [$anNOEpotential peaks]
    }

    #
    # write each NOE
    #

    foreach peakPtr $anNOElist {
	writeOneMarvinPeak $outUnit $peakPtr
    }

    close $outUnit
}

proc writeOneMarvinPeak {outUnit aPeak} {
    
    # convert pointer to an object of the same name
    
    Peak -this $aPeak
    
    puts $outUnit [format "peak %s " [$aPeak name]]
    puts $outUnit [format "   -bounds %5.2f %5.2f" [$aPeak upBound] [$aPeak lowBound]]

    if {[$aPeak hasIntensity]} {
	puts $outUnit [format "   -intensity %f" [$aPeak intensity]]
    }

    if {[$aPeak hasFromProtonShift]} {
	puts $outUnit [format "   -fromProtonShift %f" [$aPeak fromProtonShift]]
    }

    if {[$aPeak hasToProtonShift]} {
	puts $outUnit [format "   -toProtonShift %f" [$aPeak toProtonShift]]
    }

    if {[$aPeak hasFromHeavyatomShift]} {
	puts $outUnit [format "   -fromHeavyatomShift %f" [$aPeak fromHeavyatomShift]]
    }

    if {[$aPeak hasToHeavyatomShift]} {
	puts $outUnit [format "   -toHeavyatomShift %f" [$aPeak toHeavyatomShift]]
    }
    
    if {[$aPeak note] != ""} {
	foreach l [split [$aPeak note] "\n"] {
	    puts $outUnit [format "   -note %s" $l]
	}
    }
    
    foreach aPA [$aPeak peakAssignments] {
		
	PeakAssignment -this $aPA
	
	set fromAss [$aPA fromAssignment]
	set toAss   [$aPA toAssignment]

	ShiftAssignment -this $fromAss
	set fromAssName [$fromAss name]
	rename $fromAss ""

	ShiftAssignment -this $toAss
	set toAssName [$toAss name]
	rename $toAss ""

	puts $outUnit [format "   peakAssign %s %s %s" [$aPA name] $fromAssName $toAssName]
	
	if {[$aPA aveExp] != $marvinassignment::PeakAssignment_defaultAveExp} {
	    puts $outUnit [format "      -aveExp  %5.2f" [$aPA aveExp]]
	}
	
	if {[$aPA upBoundCorrection] != 0.0} {
	    puts $outUnit [format "      -upBoundCorrection  %5.2f" [$aPA upBoundCorrection]]
	}
	
	if {[$aPA lowBoundCorrection] != 0.0} {
	    puts $outUnit [format "      -lowBoundCorrection %5.2f" [$aPA lowBoundCorrection]]
	}

	if {[$aPA hasPreviousLikelihood]} {
	    puts $outUnit [format "      -likelihood %f"  \
			       [$aPA previousLikelihood]]
	}

	if {[$aPA isGood]} {
	    puts $outUnit [format "      -good"]
	}

	if {[$aPA hasUnfoldedFromProtonPeakPosition]} {
	    puts $outUnit [format "      -unfoldedFromProtonPeakPosition %f" \
			       [$aPA unfoldedFromProtonPeakPosition]]
	}

	if {[$aPA hasUnfoldedFromHeavyatomPeakPosition]} {
	    puts $outUnit [format "      -unfoldedFromHeavyatomPeakPosition %f" \
			       [$aPA unfoldedFromHeavyatomPeakPosition]]
	}

	if {[$aPA hasUnfoldedToProtonPeakPosition]} {
	    puts $outUnit [format "      -unfoldedToProtonPeakPosition %f" \
			       [$aPA unfoldedToProtonPeakPosition]]
	}

	if {[$aPA hasUnfoldedToHeavyatomPeakPosition]} {
	    puts $outUnit [format "      -unfoldedToHeavyatomPeakPosition %f" \
			       [$aPA unfoldedToHeavyatomPeakPosition]]
	}

	puts $outUnit [format "      -numFiltersFailed %d" [$aPA numFiltersFailed]]

	if {[$aPA note] != ""} {
	    foreach l [split [$aPA note] "\n"] {
		puts $outUnit [format "      -note %s" $l]
	    }
	}

	foreach {pName paName} [$aPA linkedPeakAssignmentNames] {
	    puts $outUnit [format "      -linkedTo %s %s" $pName $paName]
	}

	puts $outUnit "   end"
	rename $aPA ""
    }
    
    puts $outUnit "end"
    puts $outUnit ""
    
    rename $aPeak ""
}


proc writeXplorNOEs args {
    
    set fileName            [requiredFlagVal $args -fileName]
    set remarksList         [flagVal $args -remarks]
    set anNOEpotential      [flagVal $args -pot]
    set aPeakList           [flagVal $args -peakList]
    set noRemoveIntersected [flagExists $args -noRemoveIntersected]


    if {($anNOEpotential == "") && ($aPeakList == "")} {
	error "Neither -pot or -peakList defined in call to writeXplorNOEs"
    }

    if {$anNOEpotential != ""} {
	set aPeakList [$anNOEpotential peaks]
    }

    #create a version we can pass to Python
    set peakList [CDSList_PeakPtr]
    $peakList fromList $aPeakList
    marvinPyth command {from pyInterp import fromStringRep}
    set ptr [$peakList smartPtr]
    $peakList -disown
    marvinPyth command "tmpStrPtr='$ptr'"
    #Exceptionally ugly hack to get Python SWIG side to understand
    #what this string pointer representation is...
    marvinPyth command {tmpStrPtr=tmpStrPtr[:-1]+'10_t'}

    #there may be a (no-module) object named pasdPot- don't want to
    #destroy that.
    marvinPyth command {if not "pasdPot" in dir(): import pasdPot}

    marvinPyth command "peakList=fromStringRep(tmpStrPtr)" 

    set remarksCDSList [CDSList_String]
    $remarksCDSList fromList $remarksList
    set ptr [$remarksCDSList smartPtr]
    $ptr -disown
    $remarksCDSList -disown
    marvinPyth command "tmpStrPtr='$ptr'"
    #Exceptionally ugly hack to get Python SWIG side to understand
    #what this string pointer representation is...
    marvinPyth command {tmpStrPtr=tmpStrPtr[:-1]+'10_t'}
    # the next line for CDSList_String
    marvinPyth command "import atomSelAction"
    marvinPyth command "remarksList=fromStringRep(tmpStrPtr)"

#    # this logic would mean that nested braces will be messed up
#       that's why we're not using it...
#    marvinPyth command \
#	"remarksList=tclList.translate(None,'\{').split('\}')" \
#	[list [list tclList $remarksList]]
    

    marvinPyth command {from pasd import noeTools}
    marvinPyth command "noeTools.writeXplorAssignments(filename,
			       remarksList,
			       peaks=peakList,
	               removeIntersected=True if noRemove==\"0\" else False)" \
	[list [list filename $fileName] \
	     [list noRemove $noRemoveIntersected]] {}
    
}





proc convertXplorNOEs {inXplorName outPeakName outShiftAssignmentName peakNamePrefix} {

    set temp [create_MarvinNOEPotential tempPot]

    readXplorNOEs -fileName $inXplorName -pot $temp -namePrefix $namePrefix

    writeShiftAssignments -fileName $outShiftAssignmentName -pot $temp -remarks [list [format "autoconverted from %s" $inFileName]]

    writeMarvinPeaks -fileName $outFileName -pot $temp -remarks [list [format "autoconverted from %s" $inFileName]]
    
    rename $temp ""
}

#
# Use Monte Carlo to pick a new set of active peakAssignments.
# Hand in the number of MC steps to take, and the characteristicDeltaFactor,
# which is used to control the MC temperatures.  Each characteristicDelta* value
# in the MarvinNOEPot is scaled from its current value to (curVal * characteristicDeltaFactor)
# linearly over the course of the MC annealing
#

proc updateActivationByMonteCarlo args {
    
    set numMCsteps [flagVal $args -numMCsteps 1] 
    set CDF        [flagVal $args -characteristicDeltaFactor 1.0]
    set pot        [requiredFlagVal $args -pot]
    set isVerb     [flagExists $args -verbose]

    global errorInfo

    set origCharDeltaDV [$pot characteristicDeltaDV]
    set origCharDeltaPL [$pot characteristicDeltaPL]
    set origCharDeltaNC [$pot characteristicDeltaNoeCompleteness]
    set origCharDeltaPS [$pot characteristicDeltaScatter]
    
    set finalCharDeltaDV [expr $origCharDeltaDV * $CDF]
    set finalCharDeltaPL [expr $origCharDeltaPL * $CDF]
    set finalCharDeltaNC [expr $origCharDeltaNC * $CDF]
    set finalCharDeltaPS [expr $origCharDeltaPS * $CDF]

    if {$isVerb} {
	puts [format "Updating NOE activation in MarvinNOEPotential %s by Monte Carlo" [$pot instanceName]]
	puts ""
    }

    
    for {set stepCount 0} {$stepCount < $numMCsteps} {incr stepCount} {
	
	if {$isVerb} {
	    puts [format "MC step %d of %d:" $stepCount $numMCsteps]
	}

	$pot setCharacteristicDeltaDV   [linearScale $origCharDeltaDV $finalCharDeltaDV $numMCsteps $stepCount]
	$pot setCharacteristicDeltaPL   [linearScale $origCharDeltaPL $finalCharDeltaPL $numMCsteps $stepCount]
	$pot setCharacteristicDeltaNoeCompleteness \
	    [linearScale $origCharDeltaNC $finalCharDeltaNC $numMCsteps $stepCount]
	$pot setCharacteristicDeltaScatter \
	    [linearScale $origCharDeltaPS $finalCharDeltaPS $numMCsteps $stepCount]

	#
	# updateActivation returns a list of numbers that describe the 
	# move it just took. It throws an error if the number of 
	# unaccepted MC tries is too large, which usually means we're
	# at a local minimum
	#

	if { [catch {set curSummary [$pot updateActivation]}] } {
	    if {$isVerb} {
		puts "Gave up after too many unaccepted MC tries.  Reverting to previous step's activations"
	    }
	    set errorInfo ""
	    break
	}

	if {$isVerb} {
	    
	    puts [format "   violation score                      %.5f --> %.5f" \
		      [lindex $curSummary 0] [lindex $curSummary 1]]

	    puts [format "   previous likelihood score            %.5f --> %.5f" \
		      [lindex $curSummary 2] [lindex $curSummary 3]]

	    puts [format "   NOE completeness score               %.5f --> %.5f" \
		      [lindex $curSummary 4] [lindex $curSummary 5]]

	    puts [format "   peak posn scatter score              %.5f --> %.5f" \
		      [lindex $curSummary 6] [lindex $curSummary 7]]

	    puts [format "   normal energy                        %.5f --> %.5f" \
		      [lindex $curSummary 8] [lindex $curSummary 9]]

	    puts [format "   inverse energy                       %.5f --> %.5f" \
		      [lindex $curSummary 10] [lindex $curSummary 11]]

	    puts [format "   number of inactive peaks        %5d --> %5d" \
		      [expr round([lindex $curSummary 12])] [expr round([lindex $curSummary 13])]]

	    puts [format "   number of shift assignments in wrong neighborhoods %5d --> %5d" \
		      [expr round([lindex $curSummary 14])] [expr round([lindex $curSummary 15])]]

	    puts [format "   fraction of peak assignments inactive %.5f --> %.5f" \
		      [lindex $curSummary 16] [lindex $curSummary 17]]

	    if {([lindex $curSummary 18] > 0) || ([lindex $curSummary 19] > 0)} {
		puts [format "   frac active longrange forces that are good  %.5f --> %.5f" \
			  [lindex $curSummary 18] [lindex $curSummary 19]]
		
		puts [format "   frac good longrange peak assignments active %.5f --> %.5f" \
			  [lindex $curSummary 20] [lindex $curSummary 21]]
	    }
	    
	    puts [format "   Took %d tries, overall likelihood = %.5f" \
		      [expr round([lindex $curSummary 22])] [lindex $curSummary 23]]
	    puts ""
	    puts ""
	}
    }

    $pot setCharacteristicDeltaDV              $origCharDeltaDV
    $pot setCharacteristicDeltaPL              $origCharDeltaPL
    $pot setCharacteristicDeltaNoeCompleteness $origCharDeltaNC
    $pot setCharacteristicDeltaScatter         $origCharDeltaPS

}



proc writeActivePeakAssigns args {

    set fileName   [requiredFlagVal $args -fileName]
    set potential  [requiredFlagVal $args -pot]

    if { [catch {set outUnit [open $fileName w]}] } {
    	error "Error opening output file $fileName"
    }
    
    foreach curPeak [$potential peaks] {

	Peak -this $curPeak

	foreach curPA [$curPeak peakAssignments] {

	    PeakAssignment -this $curPA

	    if {[$curPA isActive]} {

		puts $outUnit [format "%s %s" [$curPeak name] [$curPA name]]
	    }
	    
	    rename $curPA ""
	}
	rename $curPeak ""
    }

    close $outUnit
}


proc readAndActivatePeakAssignsFromFile args {

    
    set fileName   [requiredFlagVal $args -fileName]
    set pot        [requiredFlagVal $args -pot]

    global errorInfo

    if { [catch {set inUnit [open $fileName r]}] } {
    	error "Error opening input file $fileName"
    }

    $pot inactivateAllAssigns

    while {![eof $inUnit]} {

	if {[catch {set temp [nextLine $inUnit]}]} {
	    set errorInfo ""
	    break
	}

	set curPeakName   [lindex $temp 0]
	set curPAName     [lindex $temp 1]

	if {($curPeakName == "") || ($curPAName == "")} {
	    continue
	}

	set curPeak [$pot peakNamed $curPeakName]
	Peak -this $curPeak

	set curPA [$curPeak assignmentNamed $curPAName]
	PeakAssignment -this $curPA

	$curPA activate

	rename $curPA ""
	rename $curPeak ""
    }

    close $inUnit
}


proc printMatchingPeaks args {

    set pot                     [requiredFlagVal $args -pot]
    set namesOnly               [flagExists $args -namesOnly]
    set unassignedOnly          [flagExists $args -unassignedOnly]
    set peakName                [flagVal $args -name]
    set shiftAssignName         [flagVal $args -shiftAssignName]
    set fromProtShiftRange      [flagVal $args -fromProtonShiftRange]
    set fromHeavyShiftRange     [flagVal $args -fromHeavyatomShiftRange]
    set toProtShiftRange        [flagVal $args -toProtonShiftRange]

    set matchingPeaks [$pot peaks]

    if {$unassignedOnly} {

	set temp [list]

	foreach r $matchingPeaks {

	    Peak -this $r

	    if {[$r numPeakAssignments] == 0} {
		lappend temp $r
	    }

	    rename $r ""
	}

	set matchingPeaks $temp
    }


    if {$peakName != ""} {

	set temp [list]

	foreach p $matchingPeaks {

	    Peak -this $p

	    if {[string match $peakName [$p name]]} {
		lappend temp $p 
	    }

	    rename $p ""
	}

	set matchingPeaks $temp
    }


    if {$shiftAssignName != ""} {

	set temp [list]
	
	foreach p $matchingPeaks {
	    
	    Peak -this $p
	    set matched 0

	    foreach pa [$p peakAssignments] {
		PeakAssignment -this $pa
		set fsa [$pa fromAssignment]
		set tsa [$pa toAssignment]
		ShiftAssignment -this $fsa
		ShiftAssignment -this $tsa

		if {[string match $shiftAssignName [$fsa name]] || [string match $shiftAssignName [$tsa name]]} {
		    lappend temp $p
		    set matched 1
		}

		rename $fsa ""
		rename $tsa ""
		rename $pa ""

		if {$matched == 1} {
		    break
		}
	    }
	    rename $p ""
	}

	set matchingPeaks $temp
    }


    if {$fromProtShiftRange != ""} {

	set minShift [min $fromProtShiftRange]
	set maxShift [max $fromProtShiftRange]

	set temp [list]

	foreach p $matchingPeaks {

	    Peak -this $p
	    
	    if {([$p fromProtonShift] >= $minShift) && ([$p fromProtonShift] <= $maxShift)} {
		lappend temp $p
	    }

	    rename $p ""
	}

	set matchingPeaks $temp
    }

    
    if {$fromHeavyShiftRange != ""} {

	set minShift [min $fromHeavyShiftRange]
	set maxShift [max $fromHeavyShiftRange]

	set temp [list]

	foreach p $matchingPeaks {

	    Peak -this $p
	    
	    if {([$p fromHeavyatomShift] >= $minShift) && ([$p fromHeavyatomShift] <= $maxShift)} {
		lappend temp $p
	    }

	    rename $p ""
	}

	set matchingPeaks $temp
    }

    
    if {$toProtShiftRange != ""} {

	set minShift [min $toProtShiftRange]
	set maxShift [max $toProtShiftRange]

	set temp [list]

	foreach p $matchingPeaks {

	    Peak -this $p
	    
	    if {([$p toProtonShift] >= $minShift) && ([$p toProtonShift] <= $maxShift)} {
		lappend temp $p
	    }

	    rename $p ""
	}

	set matchingPeaks $temp
    }

    
    
    if {$namesOnly} {

	foreach p $matchingPeaks {
	    Peak -this $p
	    puts [$p name]
	    rename $p ""
	}

    } else {
	
	foreach p $matchingPeaks {
	    Peak -this $p

	    puts [format "Peak %s" [$p name]]
	    puts [format "  peak position %f %f %f" [$p fromProtonShift] [$p fromHeavyatomShift] [$p toProtonShift]]
	    
	    if {[$p note] != ""} {
		foreach l [split [$p note] "\n"] {
		    puts [format "      note %s" $l]
		}
	    }
	    
	    foreach pa [$p peakAssignments] {
		PeakAssignment -this $pa
		set fsa [$pa fromAssignment]
		set tsa [$pa toAssignment]
		ShiftAssignment -this $fsa
		ShiftAssignment -this $tsa

		puts -nonewline [format "   peakAssign %s    %s %s   " [$pa name] [$fsa name] [$tsa name]]
		
		if {[$pa isGood]} {
		    puts "GOOD"
		} else {
		    puts ""
		}

		rename $fsa ""
		rename $tsa ""
		rename $pa ""
	    }

	    puts ""
	    rename $p ""
	}
    }
}
	  
	    
proc writeExplicitInverseExceptions args {
    
    set pot   [requiredFlagVal $args -pot]
    set fname [requiredFlagVal $args -filename]
    
    if { [catch {set outUnit [open $fname w]}] } {
    	error "Error opening output file $fname"
    }

    foreach {fromname toname} [$pot explicitInverseExceptions] {
	puts $outUnit [format "except %s %s" $fromname $toname]
    }
    
    close $outUnit
}


proc readExplicitInverseExceptions args {

    global errorInfo

    set pot   [requiredFlagVal $args -pot]
    set fname [requiredFlagVal $args -fileName]
    
    if { [catch {set inUnit [open $fname r]}] } {
    	error "Error opening input file $fname"
    }

    $pot removeAllExplicitInverseExceptions

    set count 0

    while {![eof $inUnit]} {

	if {[catch {eatXplorTextUpTo $inUnit "except"}]} {
	    set errorInfo ""
	    break
	}
	
	set saName1 [nextWord $inUnit]
	set saName2 [nextWord $inUnit]

	$pot addExplicitInverseException $saName1 $saName2

	updateUser [format "Reading explicit inverse exception %d \r" [incr count]]
    }

    close $inUnit
}


proc consolidateClosePeaks args {

    #
    # combines peaks that are very close to each other--used in cases where I have 
    # a single spectrum picked in several ways (eg., two different programs, or 
    # by hand and automatically)
    #

    set pot                [requiredFlagVal $args -pot]
    set protonTolerance    [flagVal $args -protonTolerance 0.01]
    set heavyatomTolerance [flagVal $args -heavyatomTolerance 0.05]
    set remVar             [flagVal $args -remarksVariableName ""]

    set peakData [list]
    set consolidated [list]

    foreach curPeak [$pot peaks] {
	Peak -this $curPeak
	lappend peakData [list [$curPeak name] [$curPeak fromProtonShift] [$curPeak fromHeavyatomShift] [$curPeak toProtonShift]]
	set deleted([$curPeak name]) 0
	rename $curPeak ""
    }

    for {set curPos 0} {$curPos < [llength $peakData]} {incr curPos} {
	set curElem       [lindex $peakData $curPos]
	set curPeakName   [lindex $curElem 0]
	set curFromProton [lindex $curElem 1]
	set curFromHeavy  [lindex $curElem 2]
	set curToProton   [lindex $curElem 3]

	updateUser [format "Looking for peaks too close to %s (%d of %d)\r" $curPeakName $curPos [llength $peakData]]

	if {$deleted($curPeakName) == 1} {
	    continue
	}
	
	for {set otherPos [expr $curPos + 1]} {$otherPos < [llength $peakData]} {incr otherPos} {
	    set otherElem       [lindex $peakData $otherPos]
	    set otherPeakName   [lindex $otherElem 0]
	    set otherFromProton [lindex $otherElem 1]
	    set otherFromHeavy  [lindex $otherElem 2]
	    set otherToProton   [lindex $otherElem 3]

	    if {$deleted($otherPeakName) == 1} {
		continue
	    }

	    if {([expr abs($curFromProton - $otherFromProton)] < $protonTolerance) &&
		([expr abs($curFromHeavy - $otherFromHeavy)] < $heavyatomTolerance) &&
		([expr abs($curToProton - $otherToProton)] < $protonTolerance)} {

		set curPeak [$pot peakNamed $curPeakName]
		Peak -this $curPeak
		set otherPeak [$pot peakNamed $otherPeakName]
		Peak -this $otherPeak
		$curPeak appendToNote [format "Consolidating peak %s into this one--here are its notes" $otherPeakName]
		$curPeak appendToNote [$otherPeak note]

		lappend consolidated [list [$curPeak name] [$otherPeak name]]

		rename $curPeak ""
		rename $otherPeak ""
		set deleted($otherPeakName) 1

	    }
	}
    }

    foreach curPeakName [array names deleted] {
	if {$deleted($curPeakName) == 1} {
	    set curPeak [$pot peakNamed $curPeakName]
	    Peak -this $curPeak
	    $pot removePeakNamed $curPeakName
	    $curPeak -acquire
	    rename $curPeak ""
	}
    }

    if {$remVar != ""} {
	upvar $remVar tempRem
	
	set nConsolidated [llength $consolidated]
	
	set rpt [format "ConsolidatePeaks: %d peaks were removed for being too close to existing peaks" $nConsolidated]
	appendLineToString rpt [format "   proton shift tolerance %f" $protonTolerance]
	appendLineToString rpt [format "   heavyatom shift tolerance %f" $heavyatomTolerance]

	if {$nConsolidated > 0} {
	    appendLineToString rpt "They are: "
	    foreach elem $consolidated {
		appendLineToString rpt [format "   %s --> %s" [lindex $elem 1] [lindex $elem 0]]
	    }
	}

	lappend tempRem $rpt
    }
    return ""
}



proc removeFilteredPeakAssignments args {

    set pot    [requiredFlagVal $args -pot]
    set cutoff [flagVal $args -cutoff 0]
    set remVar [flagVal $args -remarksVariableName ""]
    set preventUnassignment [flagExists $args -preventUnassignment]

    set nPeaksAffected 0
    set nPAsRemoved 0
    set nGoodPAsRemoved 0
    set nLRPAsRemoved 0
    set nGoodLRPAsRemoved 0


    foreach curPeak [$pot peaks] {
	Peak -this $curPeak

	set curPeakAffected 0

	if {$preventUnassignment} {

	    # find the lowest number of failed filters of this peak's peakAssignments

	    set minFailedFilters 999
	    
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		
		set minFailedFilters [min $minFailedFilters [$curPA numFiltersFailed]]
		rename $curPA ""
	    }
	    
	    set curCutoff [max $cutoff $minFailedFilters]

	} else {
	    set curCutoff $cutoff
	}

	foreach curPA [$curPeak peakAssignments] {
	    PeakAssignment -this $curPA

	    if {[$curPA numFiltersFailed] > $curCutoff} {

		incr nPAsRemoved
		
		if {[$curPA isLongRange 6]} {
		    incr nLRPAsRemoved 
		    if {[$curPA isGood]} {
			incr nGoodLRPAsRemoved 
		    } 
		} 

		if {[$curPA isGood]} {
		    incr nGoodPAsRemoved
		}
		
		set curPeakAffected 1

		$curPeak removePeakAssignmentNamed [$curPA name]
		$curPA -acquire
	    }
	    
	    rename $curPA ""
	}

	if {$curPeakAffected} {
	    incr nPeaksAffected
	}

	rename $curPeak ""
    }

    # 
    # eliminate any remaining links to PAs I just deleted
    #

    $pot clearMissingLinksToPeakAssignments

    
    if {$remVar != ""} {
	upvar $remVar tempRem
	
	set newRem [format "Removing peakAssignments with numFailedFilters > %d" $cutoff]
	appendLineToString newRem [format "Removed %d peakAssignments of %d peaks" $nPAsRemoved $nPeaksAffected]
	appendLineToString newRem [format "Of them, %d peakAssignments were long range" $nLRPAsRemoved]
	
	if {$nGoodPAsRemoved > 0} {
	    appendLineToString newRem [format "%d of the removed peakAssigns were marked as good" $nGoodPAsRemoved]
	    appendLineToString newRem [format "Of them, %d were long range" $nGoodLRPAsRemoved]
	}

	lappend tempRem $newRem
    }

    return ""
}


proc removeLowLikelihoodPeakAssignments args {

    set pot    [requiredFlagVal $args -pot]
    set cutoff [flagVal $args -cutoff 0.9]
    set remVar [flagVal $args -remarksVariableName ""]
    set preventUnassignment [flagExists $args -preventUnassignment]
    if [flagExists $args -likelihoodCutoff] {
	error "bad argument name: -likelihoodCutoff should be -cutoff"
    }

    set nPeaksAffected 0
    set nPAsRemoved 0
    set nGoodPAsRemoved 0
    set nLRPAsRemoved 0
    set nGoodLRPAsRemoved 0


    foreach curPeak [$pot peaks] {
	Peak -this $curPeak

	set curPeakAffected 0

	if {$preventUnassignment} {

	    set curCutoff [max $cutoff [$curPeak previousLikelihood]]

	} else {
	    set curCutoff $cutoff
	}

	foreach curPA [$curPeak peakAssignments] {
	    PeakAssignment -this $curPA

	    if {[$curPA previousLikelihood] < $curCutoff} {

		incr nPAsRemoved
		
		if {[$curPA isLongRange 6]} {
		    incr nLRPAsRemoved 
		    if {[$curPA isGood]} {
			incr nGoodLRPAsRemoved 
		    } 
		} 

		if {[$curPA isGood]} {
		    incr nGoodPAsRemoved
		}
		
		set curPeakAffected 1

		$curPeak removePeakAssignmentNamed [$curPA name]
		$curPA -acquire
	    }
	    
	    rename $curPA ""
	}

	if {$curPeakAffected} {
	    incr nPeaksAffected
	}

	rename $curPeak ""
    }

    # 
    # eliminate any remaining links to PAs I just deleted
    #

    $pot clearMissingLinksToPeakAssignments

    
    if {$remVar != ""} {
	upvar $remVar tempRem
	
	set newRem [format "Removing peakAssignments with previousLikelihood < %d" $cutoff]
	appendLineToString newRem [format "Removed %d peakAssignments of %d peaks" $nPAsRemoved $nPeaksAffected]
	appendLineToString newRem [format "Of them, %d peakAssignments were long range" $nLRPAsRemoved]
	
	if {$nGoodPAsRemoved > 0} {
	    appendLineToString newRem [format "%d of the removed peakAssigns were marked as good" $nGoodPAsRemoved]
	    appendLineToString newRem [format "Of them, %d were long range" $nGoodLRPAsRemoved]
	}

	lappend tempRem $newRem
    }

    return ""
}

