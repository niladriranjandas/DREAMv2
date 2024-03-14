#
# procedures to read and write NOE stripes,
# along with code formerly in aeneas_tools.tcl for 
# detecting stripes to begin with
#
# JJK 1/7/05
#

package provide marvin 1.0
package require pasdpot

proc readMarvinNOEStripes args {
    
    set fileName   [requiredFlagVal $args -fileName]
    set potential  [requiredFlagVal $args -pot]
    
    if { [catch {set inUnit [open $fileName r]}] } {
    	error "Error opening input file $fileName"
    }
    
    XplorCommand "set message off echo off end"
    
    while {![eof $inUnit]} {
	
	if { [catch {set curStripe [readOneNOEStripe $inUnit]} msg]} {

	    if {$msg == "No more stripes"} {
		break
	    } else {
		error [format "NOE stripe reading error in file %s:\n %s" $fileName $msg]
	    }
	}

	processOneNOEStripe $curStripe $potential
	updateUser [format "read %s \r" [lindex $curStripe 0]]
    }
    
    close $inUnit
}


proc readOneNOEStripe {fileID} {

    # 
    # Try to find the beginning of another stripe, ignoring 
    # text within xplor-style comments 
    #
    
    if {[catch {eatXplorTextUpTo $fileID "stripe"}]} {

	error "No more stripes"
    }

    # defaults for optional fields 
    
    set curNote ""
    set fps "NONE"
    set fhs "NONE"
    set memberList {}

    # grab the required fields

    set curName [nextWord $fileID]
    if {$curName == ""} {
	error "Missing stripe name"
    }
    
    # search for optional fields
    
    set done 0
    while {! $done} {
	
	set curWord [nextWord $fileID]
	
	switch -exact -- $curWord {
	    
	    "" {
		error [format "unexpected end of file at stripe %s" $curName]
	    }
	    
	    "end" {
		set done 1
	    }
	    
	    "stripe" {
		error [format "missing end statement after stripe %s" $curName]
	    }
	    
	    "-member" {

		lappend memberList [nextWord $fileID] 
	    }

	    "-fromProtonShift" {
		
		set fps [nextWord $fileID]

		if {! [isNumber $fps]} {
		    error [format "Missing or invalid from proton shift for stripe %s" $curName]
		}
	    }

	    "-fromHeavyatomShift" {
		
		set fhs [nextWord $fileID]

		if {! [isNumber $fhs]} {
		    error [format "Missing or invalid from heavyatom shift for stripe %s" $curName]
		}
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
    
    #
    # make sure we actually read some members
    #
    
    if {[llength $memberList] == 0} {
	error [format "No members found for NOE stripe %s" $curName]
    }
    
    return [list $curName $curNote $fps $fhs $memberList]
}


proc processOneNOEStripe {stripeData anNOEpotential} {
    
    #
    # takes the data read by readOneMarvinNOEStripe
    # and actually creates a Marvin NOEStripe
    # and attaches it to the given MarvinNOEPotential
    #

    set stripeName [lindex $stripeData 0]
    set curNote    [lindex $stripeData 1]
    set fps        [lindex $stripeData 2]
    set fhs        [lindex $stripeData 3]
    set memberList [lindex $stripeData 4]
    
    #
    # Create the stripe and add it to the pot
    #
    
    set aStripe [NOEStripe -args $stripeName]
    
    if {[catch {$anNOEpotential addToNoeStripes [$aStripe cget -this]}]} {
	error [format "Attempt to add an NOEStripe with an existing name: %s to the Marvin NOE potential %s" \
		   $stripeName [$anNOEpotential instanceName]]
    }
    
    #
    # set the stripe's note and shifts
    #
    
    $aStripe setNote $curNote

    if {$fps != "NONE"} {
	$aStripe setFromProtonShift $fps
    }

    if {$fhs != "NONE"} {
	$aStripe setFromHeavyatomShift $fhs
    }
    
    #
    # grab the members from the pot and add them to the stripe
    #
    
    foreach memberName $memberList {
	
	if {[catch {set aPeak [$anNOEpotential peakNamed $memberName]}]} {
	    error [format "Can't find peak %s in Marvin NOE potential %s to add to NOEStripe %s" \
		       $memberName [$anNOEpotential instanceName] $stripeName]
	}
	
	Peak -this $aPeak
	$aStripe addToMembers [$aPeak cget -this]
	rename $aPeak ""
    }

    $aStripe -disown
    rename $aStripe ""
}


proc writeMarvinNOEStripes args {

    set fileName        [requiredFlagVal $args -fileName]
    set remarksList     [flagVal $args -remarks]
    set anNOEpotential  [flagVal $args -pot]
    set stripeList      [flagVal $args -stripeList]

    if {($anNOEpotential == "") && ($stripeList == "")} {
	error "Neither -pot or -stripeList defined in call to writeMarvinNOEStripes"
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
    # grab NOEStripe list from pot if it's defined
    #

    if {$anNOEpotential != ""} {
	set stripeList [$anNOEpotential noeStripes]
    }

    #
    # write each NOEStripe
    #

    foreach noeStripePtr $stripeList {
	writeOneMarvinNOEStripe $outUnit $noeStripePtr
    }

    close $outUnit
}


proc writeOneMarvinNOEStripe {outUnit aStripe} {

    NOEStripe -this $aStripe

    puts $outUnit [format "stripe %s" [$aStripe name]]

    if {[$aStripe hasFromProtonShift]} {
	puts $outUnit [format "   -fromProtonShift %f" [$aStripe fromProtonShift]]
    }

    if {[$aStripe hasFromHeavyatomShift]} {
	puts $outUnit [format "   -fromHeavyatomShift %f" [$aStripe fromHeavyatomShift]]
    }

    if {[$aStripe note] != ""} {
	foreach l [split [$aStripe note] "\n"] {
	    puts $outUnit [format "   -note %s" $l]
	}
    }

    foreach peak [$aStripe members] {
	Peak -this $peak
	puts $outUnit [format "   -member %s" [$peak name]]
	rename $peak ""
    }
    
    puts $outUnit "end"
    puts $outUnit ""
    
    rename $aStripe ""
}



proc initialStripeAnalysis args {

    set anNOEpotential  [flagVal $args -pot]
    set stripeList      [flagVal $args -stripeList]
    set remVar          [flagVal $args -remarksVariableName ""]
 
    if {($anNOEpotential == "") && ($stripeList == "")} {
	error "Neither -pot or -stripeList defined in call to initialStripeAnalysis"
    }
    
    #
    # grab NOEStripe list from pot if it's defined
    #

    if {$anNOEpotential != ""} {
	set stripeList [$anNOEpotential noeStripes]
    }

    $anNOEpotential updatePrimarySeqDists

    set remList [list]

    #
    # overall numbers
    #
    
    set peakNames [list]

    foreach s $stripeList {
	NOEStripe -this $s
	foreach p [$s members] {
	    Peak -this $p
	    lappend peakNames [$p name]
	    rename $p ""
	}
	rename $s ""
    }

    set peakNames [lsort -unique $peakNames]

    lappend remList [format "There are %d stripes, covering %d different peaks." \
			 [llength $stripeList] [llength $peakNames]]

    #
    # histogram of the number of members
    #

    set memCounts [list]
    foreach s $stripeList {
	NOEStripe -this $s
	lappend memCounts [$s numMembers]
	rename $s ""
    }

    lappend remList [histogram -data $memCounts -title "Distribution of number of members in all NOE stripes"]

    #
    # overlapped stripe report
    #

    lappend remList [checkForOverlappedStripes -stripeList $stripeList]

    #
    # unassigned member report
    #

    set assignedFracs [list]
    set assignedStripeList [list]

    foreach s $stripeList {
	NOEStripe -this $s
	set curAssignedMembers 0 
	foreach m [$s members] {
	    Peak -this $m
	    if {[$m numAssignments] > 0} {
		incr curAssignedMembers
	    }
	    rename $m ""
	}
	set fracAssigned [expr $curAssignedMembers / double([$s numMembers])]
	lappend assignedFracs $fracAssigned 

	if {$fracAssigned > 0} {
	    lappend assignedStripeList [$s cget -this]
	}

	rename $s ""
    }

    set numUnassignedStripes 0
    set numPartiallyAssignedStripes 0
    set numCompletelyAssignedStripes 0

    foreach d $assignedFracs {
	if {$d == 0} {
	    incr numUnassignedStripes
	} elseif {$d == 1} {
	    incr numCompletelyAssignedStripes
	} else {
	    incr numPartiallyAssignedStripes
	}
    }

    lappend remList [format "%d stripes are completely unassigned" $numUnassignedStripes]
    lappend remList [format "%d stripes are partially assigned" $numPartiallyAssignedStripes]
    lappend remList [format "%d stripes are completely assigned" $numCompletelyAssignedStripes]

    lappend remList [histogram -data $assignedFracs \
			 -title "Distribution of fraction of assigned members in all NOE stripes"]

    #
    # primary sequence distance analysis
    #

    set fracsLR [list]
    set assignedLRStripeList [list]

    foreach s $assignedStripeList {
	NOEStripe -this $s
	set nIntra 0
	set nSeq   0
	set nSR    0
	set nLR    0
	foreach m [$s members] {
	    if {[peakNoteIncludesString "intraresidue" $m]} {
		incr nIntra
	    } elseif {[peakNoteIncludesString "sequential" $m]} {
		incr nSeq
	    } elseif {[peakNoteIncludesString "short range" $m]} {
		incr nSR
	    } elseif {[peakNoteIncludesString "long range" $m]} {
		incr nLR
	    }
	}
	set curFracLR [expr $nLR / double($nIntra + $nSeq + $nSR + $nLR)]
	lappend fracsLR $curFracLR
	if {$curFracLR > 0} {
	    lappend assignedLRStripeList [$s cget -this]
	}
	
	rename $s ""
    }

    set nNoLR 0
    set nSomeLR 0
    set nAllLR 0

    foreach f $fracsLR {
	if {$f == 0} {
	    incr nNoLR
	} elseif {$f == 1} {
	    incr nAllLR
	} else {
	    incr nSomeLR
	}
    }

    lappend remList [format "Among %s stripes with assignments, " [llength $assignedStripeList]]
    lappend remList [format "   %d have no long range peaks" $nNoLR]
    lappend remList [format "   %d have some long range peaks" $nSomeLR]
    lappend remList [format "   %d have all long range peaks" $nAllLR]
	
    lappend remList [histogram -data $fracsLR -title "Distribution of frac long range peaks per assigned stripe"]

    #
    # good/bad peak analysis
    #

    set fracs [list]
    
    foreach s $assignedLRStripeList {
	NOEStripe -this $s
	set nGood 0
	set nBad 0
	foreach m [$s members] {
	    Peak -this $m
	    set curNAssigns [$m numPeakAssignments]
	    rename $m ""
	    if {$curNAssigns > 0} {
		if {[peakNoteIncludesString "Good" $m]} {
		    incr nGood
		} else {
		    incr nBad
		}
	    }
	}

	set nAssigned [expr $nGood + $nBad]
	if {$nAssigned > 0} {
	    set fracGood [expr $nGood / double($nGood + $nBad)]
	    lappend fracs $fracGood
	}

	rename $s ""
    }

    set nCompletelyGood 0
    set nCompletelyBad  0
    set nMixed          0

    foreach f $fracs {
	if {$f == 1} {
	    incr nCompletelyGood
	} elseif {$f == 0} {
	    incr nCompletelyBad
	} else {
	    incr nMixed
	}
    }

    lappend remList [format "Among %d stripes with assigned long range members," [llength $assignedLRStripeList]]
    lappend remList [format "   %f%% are all-good" [expr 100 * $nCompletelyGood / double([llength $fracs])]]
    lappend remList [format "   %f%% are all-bad"  [expr 100 * $nCompletelyBad / double([llength $fracs])]]
    lappend remList [format "   %f%% are mixed good-bad"  [expr 100 * $nMixed / double([llength $fracs])]]

    lappend remList [format "Frac good peaks per stripe (ignoring unassigned) %f +/- %f" \
			 [mean $fracs] [standardDeviation $fracs]]
    
    lappend remList [histogram -data $fracs -title \
			 "Distribution of frac good peaks per stripe with long range peaks"]
    
    #
    # copy output
    #

    if {$remVar != ""} {
	upvar $remVar r
	foreach elem $remList {
	    lappend r $elem
	}
    }
}
    

proc checkForOverlappedStripes args {

    set anNOEpotential  [flagVal $args -pot]
    set stripeList      [flagVal $args -stripeList]

    if {($anNOEpotential == "") && ($stripeList == "")} {
	error "Neither -pot or -stripeList defined in call to checkForOverlappedStripes"
    }
    
    #
    # grab NOEStripe list from pot if it's defined
    #

    if {$anNOEpotential != ""} {
	set stripeList [$anNOEpotential noeStripes]
    }

    #
    # create an array indexed by peak name, each 
    # element of which is a list of the names of all the stripes which 
    # include it
    #

    foreach curStripe $stripeList {

	NOEStripe -this $curStripe
	set curStripeName [$curStripe name]
	set curStripeMembers [$curStripe members]

	foreach curPeak $curStripeMembers {
	    Peak -this $curPeak
	    set curPeakName [$curPeak name]

	    lappend stripeIndex($curPeakName) $curStripeName
	    
	    rename $curPeak ""
	}
	rename $curStripe ""
    }

    #
    # translate that array into a list of pairs of overlapped stripes,
    # with each pair having its stripe names sorted (to speed up the next step)
    #
    
    set overlappedStripePairs [list]
    
    foreach peakName [array names stripeIndex] {
	set curStripeList $stripeIndex($peakName)
	if {[llength $curStripeList] > 1} {
	    
	    for {set fromCount 0} {$fromCount < [expr [llength $curStripeList] - 1]} {incr fromCount} {
		for {set toCount [expr $fromCount + 1]} {$toCount < [llength $curStripeList]} {incr toCount} {
		    lappend overlappedStripePairs [lsort -dictionary [list [lindex $curStripeList $fromCount] \
									  [lindex $curStripeList $toCount]]]
		}
	    }
	}
    }
    
    #
    # eliminate duplicates
    #

    set cleanOverlaps [lsort -dictionary -unique $overlappedStripePairs]

    #
    # return a long comment with the no-duplicate list of overlapped stripe pairs
    #

    set retVal ""

    appendLineToString retVal [format "There are %d pairs of overlapped stripes: " [llength $cleanOverlaps]]

    foreach overlap $cleanOverlaps {

	appendLineToString retVal [format "   Stripes %s and %s overlap" [lindex $overlap 0] [lindex $overlap 1]]
    }

    return $retVal
}



#
# Only identify peaks as related if two covalently connected shifts are the 
# same (within tolerance) in each peak.  
#
# Note that this can imply having to check symmetric positions in 4DCC spectra
#
# Note that the shift tolerance in higher dimensions can differ from the tolerance
# along the symmetric axis because of limited digital resolution
#


proc findNOEStripes  args {

    set anNOEpotential  [requiredFlagVal $args -pot]
    set peakList        [flagVal $args -peakList]
    set FPtol           [flagVal $args -fromProtonShiftTolerance 0.02]
    set FHtol           [flagVal $args -fromHeavyatomShiftTolerance 0.02]
    set remVar          [flagVal $args -remarksVariableName ""]
    set ignoreSign      [flagExists $args -ignoreSign]
    set minStripeLength [flagVal $args -minStripeLength 1]

    if {[llength $peakList] == 0} {
	set peakList [$anNOEpotential peaks]
    }

    set remList [list]
    lappend remList [format "Running findNOEStripes with from proton shift tolerance of %f ppm" $FPtol]
    lappend remList [format "and from heavyatom shift tolerance of %f ppm" $FHtol]

    if {$ignoreSign} {
	lappend remList "Allowing peaks to have different signs within a stripe."
    } else {
	lappend remList "Requiring peaks to have the same sign within a stripe."
    }
    
    #
    # sort the peak list by fromProtonShift to speed up search 
    #

    proc fromProtonComparison {peakA peakB} {

	if {$peakA == $peakB} {
	    return 0
	}

	Peak -this $peakA
	Peak -this $peakB

	set fpA [$peakA fromProtonShift]
	set fpB [$peakB fromProtonShift]
	
	rename $peakA ""
	rename $peakB ""
	
	if {$fpA < $fpB} {
	    return -1
	} else {
	    return 1
	}
    }
    
    set sortedPeakList [lsort -increasing -command fromProtonComparison $peakList]

    #
    # for each peak, find all the other peaks that are in an NOE stripe with it
    #
    
    set allStripes [list]
    set curCount 0
    set totCount [llength $peakList]

    foreach peakA $peakList {
	
	Peak -this $peakA
	
	set curFromProt      [$peakA fromProtonShift]
	set curFromHeavyatom [$peakA fromHeavyatomShift]
	set curSign          [sign [$peakA intensity]]
	set curName          [$peakA name]
	
	rename $peakA ""
	updateUser [format "Searching for peaks related to %s (%d of %d) \r" \
			$curName [incr curCount] $totCount]
	
	set curStripe [findNOEStripeAround $curFromProt $curFromHeavyatom $curSign $sortedPeakList \
			   $FPtol $FHtol $ignoreSign]

	set curStripe [lsort -unique $curStripe]
	
	if {[llength $curStripe] >= $minStripeLength} {
	    lappend allStripes [list [list $curFromProt $curFromHeavyatom] $curStripe]
	}
	
    }

    lappend remList [format "\n%d raw stripes found." [llength $allStripes]]

    #
    # Each stripe should have been detected multiple times, with each 
    # member peak as the starting point.  First, eliminate exact duplicates 
    # quickly using the lsort -unique function
    #

    set allStripes [lsort -unique -index 1 $allStripes]

    lappend remList [format "%d unique stripes found." [llength $allStripes]]

    #
    # now recenter each stripe:  Re-pick it, using the mean of the shifts as
    # the target instead of any particular peak's position
    #

    set cleanStripes [list]
    set recenterCount -1

    foreach curStripe $allStripes {
	
	updateUser [format "recentering stripe %d of %d \r" \
			[incr recenterCount] [llength $allStripes]]

	lappend cleanStripes [recenterStripe $anNOEpotential $curStripe $sortedPeakList $FHtol $FPtol $ignoreSign]
    }

    set cleanStripes [lsort -unique -index 1 $cleanStripes]
    
    lappend remList [format "%d unique re-centered stripes found." [llength $cleanStripes]]

    #
    # Sort the stripes by length, with longest first.
    #

    proc listLenComparison {a b} {
	set la [llength $a]
	set lb [llength $b]
	if {$la < $lb} {
	    return -1
	} else {
	    return +1
	}
    }

    set cleanStripes [lsort -command listLenComparison -decreasing -index 1 $cleanStripes]

    #
    # now check whether any existing stripe is a subset of any other existing stripe.
    # If so, eliminate it.
    #

    proc isSubset {listA listB} {

	set retVal 1

	foreach elem $listB {
	    if {[lsearch -exact $listA $elem] == -1} {
		set retVal 0
		break
	    }
	}
	
	return $retVal
    }


    for {set posA 0} {$posA < [llength $cleanStripes]} {incr posA} {
	set stripeA [lindex $cleanStripes $posA]

	updateUser [format "Checking stripe %d of %d for subset stripes \r" \
			$posA [llength $cleanStripes]]
	
	#
	# Only compare it to stripes that are shorter, since a subset has to be smaller 
	# than a superset.  Since the list is sorted by stripe length, just look at stripes
	# that come after stripeA.  Do the removal of subset stripes in a separate step
	# in order to avoid confusing the loop.
	#

	set delePosList [list]

	for {set posB [expr $posA + 1]} {$posB < [llength $cleanStripes]} {incr posB} {
	    set stripeB [lindex $cleanStripes $posB]

	    if {[isSubset [lindex $stripeA 1] [lindex $stripeB 1]]} {
		lappend delePosList $posB
	    }
	}

	set delePosList [lsort -decreasing $delePosList]

	foreach elem $delePosList {
	    set cleanStripes [lreplace $cleanStripes $elem $elem]
	}
    }

    lappend remList [format "%d superset stripes found." [llength $cleanStripes]]
    
    return $cleanStripes

    #
    # look for related stripes?
    #
    
    #
    # make sure that every peak in every stripe has compatible from selections - edge cases ?
    #
    
    #
    # create the NOEStripe, and add it to the pot
    #
    
#     set stripeCount -1

#     foreach s $cleanStripes {

# 	set curName [format "%s_%d" $namePrefix [incr stripeCount]]
# 	set curStripe [NOEStripe -args $curName]

# 	$curStripe setFromProtonShift    [lindex [lindex $s 0] 0]
# 	$curStripe setFromHeavyatomShift [lindex [lindex $s 0] 1]

# 	foreach r [lindex $s 1] {
# 	    $curStripe addToMembers [$anNOEpotential peakNamed $r]
# 	}

# 	$anNOEpotential addNoeStripe [$curStripe cget -this]

# 	$curStripe -disown
# 	rename $curStripe ""
#     }
    
#     if {$remVar != ""} {
# 	upvar $remVar rtemp
# 	foreach elem $remList {
# 	    lappend rtemp $elem
# 	}
#     }
}


#
# Given a from proton and from heavyatom shift, find all the peaks with 
# from proton and heavyatom shifts within tolerance of the target.  
#
# Returns a list of their names
#

proc findNOEStripeAround {targetProtShift targetHeavyatomShift targetSign sortedPeakList FPtol FHtol ignoreSign} {

    set retVal [list]

    foreach peakB $sortedPeakList {

	#
	# read out data
	#

	Peak -this $peakB 
	
	set curProtShift      [$peakB fromProtonShift]
	set curHeavyatomShift [$peakB fromHeavyatomShift]
	set curSign           [sign [$peakB intensity]]
	set curName           [$peakB name]

	rename $peakB ""

	#
	# a quick check to see if we're in range on the proton dimension,
	# made possible by the fact that the peak list we're given is sorted
	# in order of fromProtonShift
	#
	
	set deltaProtShift [expr $targetProtShift - $curProtShift]
	
	if {$deltaProtShift > $FPtol} {
	    continue
	} elseif {$deltaProtShift <  [expr $FPtol * -1.0]} {
	    break
	}
	
	
	#
	# see if they're in an NOE stripe, using the ellipsoid cutoff 
	# and making sure they have the same sign
	#
	
	if {[ellipseCheck \
		 [list $targetProtShift $targetHeavyatomShift] \
		 [list $curProtShift    $curHeavyatomShift] \
		 [list $FPtol $FHtol]]} {
	    
	    if {($targetSign == $curSign) || $ignoreSign} {
		
		lappend retVal $curName
	    }
	}
    }

    return $retVal
}

#
# Given a stripe of peaks, calculate a new target from proton and heavyatom
# shift from the mean values in the stripe, and then re-pick the stripe
#

proc recenterStripe {pot oldStripe sortedPeakList FHtol FPtol ignoreSign} {

    set protonShifts [list]
    set heavyatomShifts [list]
    set signs [list]
    
    foreach pName [lindex $oldStripe 1] {
	
	set p [$pot peakNamed $pName]
	Peak -this $p
	
	lappend protonShifts [$p fromProtonShift]
	lappend heavyatomShifts [$p fromHeavyatomShift]
	lappend signs [sign [$p intensity]]
	
	rename $p ""
    }

    set targetProtonShift [mean $protonShifts]
    set targetHeavyatomShift [mean $heavyatomShifts]
    set targetSign [mean $signs]
    
    set newStripe [findNOEStripeAround $targetProtonShift $targetHeavyatomShift $targetSign $sortedPeakList \
		       $FPtol $FHtol $ignoreSign]

    set newStripe [lsort -unique $newStripe]

    return [list [list $targetProtonShift $targetHeavyatomShift] $newStripe]
}



proc stripesMatchingShiftAssignment args {

    set stripeList              [requiredFlagVal $args -stripeList]
    set pot                     [requiredFlagVal $args -pot]
    set target                  [requiredFlagVal $args -target]
    set invBound                [flagVal $args -inverseBound 4.0]
    set invMethyl               [flagVal $args -inverseMethylCorrection 0.0]
    set minFracNeighborsMatched [flagVal $args -minFracNeighborsMatched 0.25]
    set minFracPeaksMatched     [flagVal $args -minFracPeaksMatched 0.25]

    #
    # if the target shiftAssignment is a toShiftAssignment, return nothing
    #

    ShiftAssignment -this $target 
	
    if {! [$target isFrom]} {
	rename $target ""
	return [list]
    }

    #
    # set up the pot for inverse NOE neighborhood calculations
    #

    $pot activateAllAssigns
    $pot setInverseBound $invBound
    $pot setInverseMethylCorrection $invMethyl
    $pot updateNoeCompleteness
    $pot updatePeakPositionScatter

    #
    # record the names of the toShiftAssignments in the target's
    # expectedNeighborhood
    #
	
    set targetsNeighbors [list]
    foreach neigh [$pot expectedNeighbors [$target cget -this]] {
	ShiftAssignment -this $neigh
	lappend targetsNeighbors [$neigh name]
	rename $neigh ""
    }

    #
    # for each stripe,
    #

    set matchingStripes [list]

    foreach stripe $stripeList {

	set curStripeShifts    [lindex $stripe 0]
	set curStripePeakNames [lindex $stripe 1]

	#
	# get sign of peaks in current stripe
	#
	
	set tempPeak [$pot peakNamed [lindex $curStripePeakNames 0]]
	Peak -this $tempPeak
	set curStripeSign [sign [$tempPeak intensity]]
	rename $tempPeak ""
	
	#
	# get frac of peaks in current stripe that are assigned
	#

	set nPeaksAssignedInCurStripe 0

	foreach peakName $curStripePeakNames {
	    set curPeak [$pot peakNamed $peakName]
	    Peak -this $curPeak
	    if {[$curPeak isAssigned]} {
		incr nPeaksAssignedInCurStripe
	    }
	    rename $curPeak ""
	}

	set curStripeFracAssigned \
	    [expr $nPeaksAssignedInCurStripe / double([llength $curStripePeakNames])]

	#
	# get from shiftAssignments that have been applied to the current stripe
	#
    
	set curStripeFromShiftAssigns [list]

	foreach peakName $curStripePeakNames {
	    set curPeak [$pot peakNamed $peakName]
	    Peak -this $curPeak
	    foreach curPA [$curPeak peakAssignments] {
		PeakAssignment -this $curPA
		set curFSA [$curPA fromAssignment]
		ShiftAssignment -this $curFSA
		lappend curStripeFromShiftAssigns [$curFSA name]
		rename $curFSA ""
		rename $curPA ""
	    }
	    rename $curPeak ""
	}

	set curStripeFromShiftAssigns [lsort -dictionary -unique $curStripeFromShiftAssigns]

	#
	# collect a set of all the possible toShiftAssigns that match the peaks
	# in this stripe
	#
        
	set curStripeToMatches [list]

	foreach peakName $curStripePeakNames {
	    
	    set p [$pot peakNamed $peakName]
	    Peak -this $p

	    set temp [list]

	    foreach toSA [$p toShiftAssignments] {
		ShiftAssignment -this $toSA
		lappend temp [$toSA name]
		rename $toSA ""
	    }

	    lappend curStripeToMatches $temp
	
	    rename $p ""
	}
	
	set uniqueToMatches [lsort -dictionary -unique [flattenList $curStripeToMatches]]


	#
	# calc fraction of the target's expectedNeighbors that
	# is explained by the stripe's peaks' toShiftAssignments
	# 

	set targetNeighborsMatched 0
	foreach neighName $targetsNeighbors {
	    if {[lsearch -exact $uniqueToMatches $neighName] != -1} {
		incr targetNeighborsMatched
	    }
	}

	if {[llength $targetsNeighbors] == 0} {
	    set fracNeighborsMatched 0
	} else {
	    set fracNeighborsMatched [expr $targetNeighborsMatched / double([llength $targetsNeighbors])]
	}

	#
	# calc fraction of this stripe's peaks that have toShiftAssigns that match 
	# the target's expectedNeighbors
	#

	set curPeaksMatched 0

	foreach peakToMatches $curStripeToMatches {

	    set curPeakMatches 0

	    foreach toSA $peakToMatches {
		if {[lsearch -exact $targetsNeighbors $toSA] != -1} {
		    set curPeakMatches 1
		    break
		}
	    }

	    if {$curPeakMatches == 1} {
		incr curPeaksMatched
	    }
	}

	if {[llength $curStripePeakNames] == 0} {
	    set fracPeaksMatched 0
	} else {
	    set fracPeaksMatched [expr $curPeaksMatched / double([llength $curStripePeakNames])]
	}
	
	#
	# if a good-enough fraction of the target's expectedNeighborhood matches the stripe,
	# and a good-enough fraction of the stripe's peaks match the target's expectedNeighborhood, 
	# then declare the stripe to match the target
	#

	if {($fracNeighborsMatched >= $minFracNeighborsMatched) &&
	    ($fracPeaksMatched >= $minFracPeaksMatched)} {

	    lappend matchingStripes [list $curStripeShifts $curStripePeakNames \
					 $curStripeSign $curStripeFracAssigned \
					 $curStripeFromShiftAssigns \
					 $fracNeighborsMatched $fracPeaksMatched \
					 [expr $fracNeighborsMatched + $fracPeaksMatched]]
	}
    }
    
    # sort list of matches by overall score, to put best ones first
    
    set matchingStripes [lsort -decreasing -real -index end $matchingStripes]
 
    return $matchingStripes
}




proc shiftAssignmentsMatchingStripe args {

    set shiftAssignList         [requiredFlagVal $args -shiftAssignmentList]
    set pot                     [requiredFlagVal $args -pot]
    set target                  [requiredFlagVal $args -target]
    set invBound                [flagVal $args -inverseBound 4.0]
    set invMethyl               [flagVal $args -inverseMethylCorrection 0.0]
    set minFracNeighborsMatched [flagVal $args -minFracNeighborsMatched 0.25]
    set minFracPeaksMatched     [flagVal $args -minFracPeaksMatched 0.25]

    #
    # set up the pot for inverse NOE neighborhood calculations
    #

    $pot activateAllAssigns
    $pot setInverseBound $invBound
    $pot setInverseMethylCorrection $invMethyl
    $pot updateNoeCompleteness
    $pot updatePeakPositionScatter

    set targetShifts    [lindex $target 0]
    set targetPeakNames [lindex $target 1]

    #
    # collect a set of all the possible toShiftAssigns that match the peaks
    # in this stripe
    #
        
    set targetToMatches [list]

    foreach peakName $targetPeakNames {

	set p [$pot peakNamed $peakName]
	Peak -this $p

	set temp [list]

	foreach toSA [$p toShiftAssignments] {
	    ShiftAssignment -this $toSA
	    lappend temp [$toSA name]
	    rename $toSA ""
	}

	lappend targetToMatches $temp
	
	rename $p ""
    }

    set uniqueToMatches [lsort -dictionary -unique [flattenList $targetToMatches]]

    #
    # look for fromShiftAssigns whose expectedNeighbors are all among the 
    # current stripe's to matches
    #

    set matchingFromShiftAssigns [list]

    foreach fromSA $shiftAssignList {
	ShiftAssignment -this $fromSA 
	
	if {! [$fromSA isFrom]} {
	    rename $fromSA ""
	    continue
	}

	#
	# record the names of the toShiftAssignments in this fromShiftAssignment's
	# expectedNeighborhood
	#

	set curNeighbors [list]
	foreach neigh [$pot expectedNeighbors [$fromSA cget -this]] {
	    ShiftAssignment -this $neigh
	    lappend curNeighbors [$neigh name]
	    rename $neigh ""
	}

	#
	# calc fraction of this fromShiftAssign's expectedNeighbors that
	# is explained by the target's peaks' toShiftAssignments
	# 

	set curNeighborsMatched 0
	foreach neighName $curNeighbors {
	    if {[lsearch -exact $uniqueToMatches $neighName] != -1} {
		incr curNeighborsMatched
	    }
	}

	if {[llength $curNeighbors] == 0} {
	    set fracNeighborsMatched 0
	} else {
	    set fracNeighborsMatched [expr $curNeighborsMatched / double([llength $curNeighbors])]
	}

	#
	# calc fraction of the target's peaks that have toShiftAssigns that match 
	# this fromShiftAssign's expectedNeighbors
	#

	set targetPeaksMatched 0

	foreach peakToMatches $targetToMatches {

	    set curPeakMatches 0

	    foreach toSA $peakToMatches {
		if {[lsearch -exact $curNeighbors $toSA] != -1} {
		    set curPeakMatches 1
		    break
		}
	    }

	    if {$curPeakMatches == 1} {
		incr targetPeaksMatched
	    }
	}

	if {[llength $targetPeakNames] == 0} {
	    set fracPeaksMatched 0
	} else {
	    set fracPeaksMatched [expr $targetPeaksMatched / double([llength $targetPeakNames])]
	}

	
	#
	# if a good-enough fraction of the fromShiftAssign's expectedNeighborhood matches the target,
	# and a good-enough fraction of the target's peaks match the fromShiftAssign's 
	# expectedNeighborhood, then declare the fromShiftAssign to match the target
	#

	if {($fracNeighborsMatched >= $minFracNeighborsMatched) &&
	    ($fracPeaksMatched >= $minFracPeaksMatched)} {

	    lappend matchingFromShiftAssigns [list [$fromSA name] \
						  [$fromSA protonShift] [$fromSA heavyatomShift] \
						  [$fromSA foldedProtonShift] \
						  [$fromSA foldedHeavyatomShift] [$fromSA foldedSign] \
						  [$fromSA noeCompleteness] \
						  $fracNeighborsMatched $fracPeaksMatched \
						  [expr $fracNeighborsMatched + $fracPeaksMatched]]
	}
    }
    
    # sort list of matches by overall score, to put best ones first
    
    set matchingFromShiftAssigns [lsort -decreasing -real -index end $matchingFromShiftAssigns]
 
    return $matchingFromShiftAssigns
}
