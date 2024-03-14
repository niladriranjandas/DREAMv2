package provide newstripe 0.1
package require marvin

namespace eval StripeFilter {

namespace export stripeCorrection

proc recordShiftAssignmentNames {pot salVarName} {

    upvar $salVarName saList

    foreach curSA [$pot shiftAssignments] {
	
	ShiftAssignment -this $curSA
	lappend saList [$curSA name]
	rename $curSA ""
    }

    set saList [lsort -dictionary $saList]

    return ""
}


proc recordToFromPartners {pot arrayName} {

    upvar $arrayName toFromPartner

    foreach curSA [$pot shiftAssignments] {
	
	ShiftAssignment -this $curSA
	set curName [$curSA name]
	
	if {[$curSA hasToFromPartnerName]} {
	    set toFromPartner($curName) [$curSA toFromPartnerName]
	} else {
	    set toFromPartner($curName) "NONE"
	}
	
	rename $curSA ""
    }
    
    return ""
}


proc recordToFromSense {pot arrayName} {

    upvar $arrayName saIsTo

    foreach curSA [$pot shiftAssignments] {
	
	ShiftAssignment -this $curSA
	set curName [$curSA name]
	
	if {[$curSA isTo]} {
	    set saIsTo($curName) 1
	} else {
	    set saIsTo($curName) 0
	}

	rename $curSA ""
    }

    return ""
}


proc recordIntraresPartners {pot arrayName} {

    upvar $arrayName intraPartnerSAs

    array unset fromSAsInResidue
    array unset   toSAsInResidue

    foreach curSA [$pot shiftAssignments] {
	ShiftAssignment -this $curSA 
	
	set intraPartnerSAs([$curSA name]) [list]
	
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
    
    foreach curResidue [array names fromSAsInResidue] {
	foreach saName $fromSAsInResidue($curResidue) {
	    if {[info exists toSAsInResidue($curResidue)]} {
		set intraPartnerSAs($saName) [concat $intraPartnerSAs($saName) $toSAsInResidue($curResidue)]
	    }
	}
    }
   
    foreach curResidue [array names toSAsInResidue] {
	foreach saName $toSAsInResidue($curResidue) {
	    if {[info exists fromSAsInResidue($curResidue)]} {
		set intraPartnerSAs($saName) [concat $intraPartnerSAs($saName) $fromSAsInResidue($curResidue)]
	    }
	}
    }

    return ""
}


proc recordBbnSeqPartners {pot arrayName} {

    upvar $arrayName intraPartnerSAs

    #
    # add shift assignments to a given SA's list of intraresidue partners if
    # both SAs are (hb* or ha* or hn) and they're in sequential residues
    #
    # Do not add bbn-seq SAs to the list of intraresidue partners for every SA in 
    # a residue, since things well out in the sidechains are unlikely to see 
    # the neighbors' bbn SAs
    #



    array unset fromBbnSAsInResidue
    array unset   toBbnSAsInResidue

    set bbnSeqSel [AtomSel -args "(name hn or name ha or name ha* or name hb or name hb*)"]

    #
    # record names of all backbone SAs, along with the residues they belong to
    #

    set fromBbnSAs [list]
    set toBbnSAs [list]

    foreach curSA [$pot shiftAssignments] {
	ShiftAssignment -this $curSA 
	
	set curSel [$curSA protonSelection]
	AtomSel -this $curSel
	
	if {[$curSel intersects $bbnSeqSel]} {
	
	    if {[$curSA isFrom]} {
		lappend fromBbnSAs [$curSA name] [residuesInSelection $curSel]
	    } else {
		lappend toBbnSAs   [$curSA name] [residuesInSelection $curSel]
	    }
	}

	rename $curSel ""
	rename $curSA ""
    }

    #
    # find pairs of from and to backbone SAs that are sequential 
    # and add them to each other's lists of intraresidue partner SAs
    #

    foreach {fromSA fromResList} $fromBbnSAs {
	foreach fromRes $fromResList {
	    set prevTarg [previousResidue $fromRes]
	    set nextTarg [nextResidue $fromRes]

	    foreach {toSA toResList} $toBbnSAs {
		foreach toRes $toResList {

		    if {($toRes == $prevTarg) || ($toRes == $nextTarg)} {
			lappend intraPartnerSAs($fromSA) $toSA
			lappend intraPartnerSAs($toSA)   $fromSA
		    }
		}
	    }
	}
    }



    return ""
}


proc removeToFromPartnersFromIntraPartners {ipArrayName tfpArrayName salVarName} {

    upvar $ipArrayName  intraPartnerSAs 
    upvar $tfpArrayName toFromPartner
    upvar $salVarName   saList

    foreach curSA $saList {

	set cleanIntraPartnerSAs [list]
	foreach partnerSA $intraPartnerSAs($curSA) {
	    if {$partnerSA != $toFromPartner($curSA)} {
		lappend cleanIntraPartnerSAs $partnerSA
	    }
	}
	
	set intraPartnerSAs($curSA) [lsort -unique $cleanIntraPartnerSAs]
    }
}


proc recordGeminalPartners {pot arrayName} {

    upvar $arrayName geminalPartner

    foreach curSA [$pot shiftAssignments] {
	ShiftAssignment -this $curSA 
	set curSAname [$curSA name]
	set curSAisFrom [$curSA isFrom]

	set geminalPartner($curSAname) "NONE"

	if {(! [$curSA hasHeavyatomSelection]) ||
	    (! [$curSA hasHeavyatomShift])} {
	    rename $curSA ""
	    continue
	}

	set curSAheavySel [$curSA heavyatomSelection]
	set curSAheavyShift [$curSA heavyatomShift]
	AtomSel -this $curSAheavySel
	rename $curSA ""

	foreach otherSA [$pot shiftAssignments] {
	    ShiftAssignment -this $otherSA
	    set otherSAisFrom [$otherSA isFrom]
	    
	    if {(! [$otherSA hasHeavyatomSelection]) ||
		(! [$otherSA hasHeavyatomShift])} {
		rename $otherSA ""
		continue
	    }

	    if {([$otherSA name] != $curSAname) && ($otherSAisFrom == $curSAisFrom)} {
		
		if {[$curSAheavySel isEqualTo [$otherSA heavyatomSelection]] &&
		    ($curSAheavyShift == [$otherSA heavyatomShift])} {

		    set geminalPartner($curSAname) [$otherSA name]
		    rename $otherSA ""
		    break
		}
	    }
	    rename $otherSA ""
	}
	
	rename $curSAheavySel ""
    }

    return ""
}


proc recordStereoPartners {pot arrayName} {

    upvar $arrayName stereoPartner

    array unset stereoPartner

    foreach curSA [$pot shiftAssignments] {
	ShiftAssignment -this $curSA 

	if {[$curSA hasStereoPartnerName]} {
	    set stereoPartner([$curSA name]) [$curSA stereoPartnerName]
	} else {
	    set stereoPartner([$curSA name]) "NONE"
	}

       	rename $curSA ""
    }

    return ""
}


proc recordAllPossibleTargetPAs {pot ptArrayName ipArrayName plVarName salVarName} {

    upvar $ptArrayName possibleTargetsForSA
    upvar $ipArrayName intraPartnerSAs
    upvar $plVarName   peakList
    upvar $salVarName  saList

    array unset possibleTargetsForSA

    foreach curSA $saList {
	set possibleTargetsForSA($curSA) [list]
    }

    set peakList [list]

    foreach p [$pot peaks] {

	Peak -this $p

	lappend peakList [$p name]
		
	foreach pa [$p peakAssignments] {
	    
	    PeakAssignment -this $pa

	    set curUnfoldedFromProtonPP "NONE"
	    set curUnfoldedFromHeavyPP "NONE"
	    set curUnfoldedToProtonPP "NONE"
	    set curUnfoldedToHeavyPP "NONE"
	    
	    if {[$pa hasUnfoldedFromProtonPeakPosition]} {
		set curUnfoldedFromProtonPP [$pa unfoldedFromProtonPeakPosition]
	    } 
	    
	    if {[$pa hasUnfoldedFromHeavyatomPeakPosition]} {
		set curUnfoldedFromHeavyPP [$pa unfoldedFromHeavyatomPeakPosition]
	    } 
	    
	    if {[$pa hasUnfoldedToProtonPeakPosition]} {
		set curUnfoldedToProtonPP [$pa unfoldedToProtonPeakPosition]
	    } 
	    
	    if {[$pa hasUnfoldedToHeavyatomPeakPosition]} {
		set curUnfoldedToHeavyPP [$pa unfoldedToHeavyatomPeakPosition]
	    } 
	   
	    
	    # Record it if this is a PA of interest

	    set curFromSA [$pa fromAssignmentName]
	    set curToSA   [$pa   toAssignmentName]


	    if {[lsearch -exact $intraPartnerSAs($curFromSA) $curToSA] != -1} {

		set curSApair [join [list $curFromSA $curToSA]]
		lappend possibleTargetsForSA($curFromSA) [list [$p name] [$pa name] $curUnfoldedFromProtonPP $curUnfoldedFromHeavyPP $curToSA $curUnfoldedToProtonPP $curUnfoldedToHeavyPP $curSApair]
	    }

	    if {[lsearch -exact $intraPartnerSAs($curToSA) $curFromSA] != -1} {

		set curSApair [join [list $curToSA $curFromSA]]
		lappend possibleTargetsForSA($curToSA) [list [$p name] [$pa name] $curUnfoldedToProtonPP $curUnfoldedToHeavyPP $curFromSA $curUnfoldedFromProtonPP $curUnfoldedFromHeavyPP $curSApair]
	    }
	    
	    rename $pa ""
	}
	
	rename $p ""
    }

    return ""
}


proc recordSAsWithNoTargets {ptArrayName ztArrayName salVarName} {

    upvar $ptArrayName  possibleTargetsForSA 
    upvar $ztArrayName  saHasZeroTargets
    upvar $salVarName   saList

    array unset saHasZeroTargets 

    foreach curSA $saList {

	if {[llength $possibleTargetsForSA($curSA)] == 0} {
	    set saHasZeroTargets($curSA) 1
	} else {
	    set saHasZeroTargets($curSA) 0
	}
    }

    return ""
}


proc targetIsReasonable {curPossTarget curSAName ptArrayName tfArrayName gpArrayName spArrayName ztArrayName salVarName protonTol heavyTol} {
 
    upvar $ptArrayName  possibleTargetsForSA 
    upvar $tfArrayName  toFromPartner
    upvar $gpArrayName  geminalPartner
    upvar $spArrayName  stereoPartner
    upvar $ztArrayName  saHasZeroTargets
    upvar $salVarName   saList

    #
    # if a PA is going to be a successful target for its SA and its PA partner's SA,
    # it needs to be compatible with any to-from, geminal, and stereo partners they have.  
    #
    # Note that geminal or stereo partner SAs with zero raw targets aren't counted against a given SA, 
    # since they could arise from mistakes in those assignments.  But if the to-from partner SA
    # has zero raw targets, we do count that against us, since that arises from the same entry 
    # in the shift table
    #

    
    set curPeak        [lindex $curPossTarget 0]

    foreach \
	saName   [list $curSAName                [lindex $curPossTarget 4]] \
	curProt  [list [lindex $curPossTarget 2] [lindex $curPossTarget 5]] \
	curHeavy [list [lindex $curPossTarget 3] [lindex $curPossTarget 6]] {

	    if {$toFromPartner($saName) == "NONE"} {
		set curToFromOK 1
	    } elseif {$saHasZeroTargets($toFromPartner($saName))} {
		set curToFromOK 0
	    } else {
		
		set curToFromOK 0

		foreach partner $possibleTargetsForSA($toFromPartner($saName)) {
	    
		    set partnerPeak  [lindex $partner 0]
		    set partnerProt  [lindex $partner 2]
		    set partnerHeavy [lindex $partner 3]
	    
		    if {($curPeak != $partnerPeak) && 
			[insideTolerance $curProt $partnerProt $protonTol] &&
			[insideTolerance $curHeavy $partnerHeavy $heavyTol]} {
		
			set curToFromOK 1
		    } 
		}
	    } 

	    if {! $curToFromOK} {
		return 0
	    }
    
    
	    if {$geminalPartner($saName) == "NONE"} {
		set curGeminalOK 1
	    } elseif {$saHasZeroTargets($geminalPartner($saName))} {
		set curGeminalOK 1
	    } else {
		
		set curGeminalOK 0

		foreach partner $possibleTargetsForSA($geminalPartner($saName)) {
	    
		    set partnerPeak  [lindex $partner 0]
		    set partnerHeavy [lindex $partner 3]
		    
		    if {($curPeak != $partnerPeak) && 
			[insideTolerance $curHeavy $partnerHeavy $heavyTol]} {
		
			set curGeminalOK 1
		    }
		}
	    } 

	    if {! $curGeminalOK} {
		return 0
	    }


	    if {$stereoPartner($saName) == "NONE"} {
		set curStereoOK 1
	    } elseif {$saHasZeroTargets($stereoPartner($saName))} {
		set curStereoOK 1
	    } else {
		
		set curStereoOK 0

		foreach partner $possibleTargetsForSA($stereoPartner($saName)) {
	    
		    set partnerPeak  [lindex $partner 0]
		    set partnerProt  [lindex $partner 2]

		    if {($curPeak != $partnerPeak) && 
			[outsideTolerance $curProt $partnerProt $protonTol]} {
			
			set curStereoOK 1
		    }
		}
	    }

	    if {! $curStereoOK} {
		return 0
	    }
	}

    return 1
}



proc filterTargetsForInitialReasonability {ptArrayName tfArrayName gpArrayName spArrayName ztArrayName salVarName protonTol heavyTol} {

    upvar $ptArrayName  possibleTargetsForSA 
    upvar $tfArrayName  toFromPartner
    upvar $gpArrayName  geminalPartner
    upvar $spArrayName  stereoPartner
    upvar $ztArrayName  saHasZeroTargets
    upvar $salVarName   saList

    for {set x 0} {$x < 10} {incr x} {

	foreach curSA [shuffledList $saList] {
	    
	    set reasonableTargs [list]
	    
	    foreach possTarg [shuffledList $possibleTargetsForSA($curSA)] {
		
		if {[targetIsReasonable $possTarg $curSA possibleTargetsForSA toFromPartner geminalPartner stereoPartner saHasZeroTargets saList $protonTol $heavyTol]} {
		    lappend reasonableTargs $possTarg
		}
	    }

	    set possibleTargetsForSA($curSA) $reasonableTargs
	}
    }
    
    return ""
}





proc calcExpectedPeaks {ipArrayName epArrayName ztArrayName salVarName} {

    upvar $ipArrayName intraPartnerSAs 
    upvar $epArrayName expectedPeaksForSA 
    upvar $ztArrayName saHasZeroTargets
    upvar $salVarName  saList

    #
    # expectedPeaksForSA is just a list of each SA's intraresidue partners, 
    # ignoring SAs that have zero targets available
    #
    # As such, it is the correct factor to compare the actual compatible 
    # partners lists against, in calculating each SA's individual score.
    #
 
    array unset expectedPeaksForSA

    foreach curSA $saList {
	
	set expectedPeaksForSA($curSA) [list]

	if {! $saHasZeroTargets($curSA)} {
	    foreach partnerSA $intraPartnerSAs($curSA) {
		if {! $saHasZeroTargets($partnerSA)} {
		    lappend expectedPeaksForSA($curSA) $partnerSA
		}
	    }
	}
    }

    return ""
}
	    
		
proc chooseTargets {ptArrayName ctArrayName tfArrayName gpArrayName spArrayName sitArrayName plVarName salVarName protonTol heavyTol} {

    upvar $ptArrayName  possibleTargetsForSA 
    upvar $ctArrayName  curTargetForSA
    upvar $tfArrayName  toFromPartner
    upvar $gpArrayName  geminalPartner
    upvar $spArrayName  stereoPartner
    upvar $sitArrayName saIsTo
    upvar $plVarName    peakList
    upvar $salVarName   saList

    #
    # Generate a set of targets for SAs that don't have any targets.
    #
    # Go through list of SAs in random order, choosing a target 
    # for each SA without a target at random from its list of possible targets.  
    #


    #
    # record all the peaks that have been used in existing targets
    #

    array unset peakUsedByFrom
    array unset peakUsedByTo

    foreach peakName $peakList {
	set peakUsedByFrom($peakName) 0
	set peakUsedByTo($peakName) 0
    }

    foreach curSA $saList {
	if {[info exists curTargetForSA($curSA)]} {
	    set curPeak [lindex $curTargetForSA($curSA) 0]
	    if {$curPeak != "NONE"} {
		if {$saIsTo($curSA)} {
		    set peakUsedByTo($curPeak) 1
		} else {
		    set peakUsedByFrom($curPeak) 1
		}
	    }
	}
    }

    #
    # create a list of SAs, with from SAs in randomized order first, and then to SAs in randomized order
    # (to avoid making bad choices based on the higher uncertainty of to SAs)
    #

    set fromSAlist [list]
    set toSAlist   [list]
    foreach curSA $saList {
	if {$saIsTo($curSA)} {
	    lappend toSAlist $curSA
	} else {
	    lappend fromSAlist $curSA
	}
    }

    set shuffledSAs [shuffledList $saList]

    foreach curSA $shuffledSAs {

	if {! [info exists curTargetForSA($curSA)]} {
	
	    if {$saIsTo($curSA)} {
		chooseTargetForSA possibleTargetsForSA curTargetForSA toFromPartner geminalPartner stereoPartner peakUsedByTo   $curSA $protonTol $heavyTol
	    } else {
		chooseTargetForSA possibleTargetsForSA curTargetForSA toFromPartner geminalPartner stereoPartner peakUsedByFrom $curSA $protonTol $heavyTol
	    }
	}
    }

    return ""
}




proc chooseTargetForSA {ptArrayName ctArrayName tfArrayName gpArrayName spArrayName puArrayName saName protonTol heavyTol} {

    #
    # Choose a target for a particular SA at random from its possible targets
    #
    # If the chosen target is incompatible with an already-selected target
    # value for its peakAssignment partner, to-from partner, or geminal
    # partner, choose another.  
    #
    # If the chosen SA's chosen target arises from a peak that has been 
    # used to provide a target to another SA of the same to-from sense,
    # choose another. 
    #
    # If the chosen target's proton shift is too close to the stereo partner's
    # already-selected target proton shift, choose another.
    #
    # If no target of a particular SA passes those filters, record it as a
    # list of NONEs.
    #

    upvar $ptArrayName possibleTargetsForSA
    upvar $ctArrayName curTargetForSA
    upvar $tfArrayName toFromPartner
    upvar $gpArrayName geminalPartner
    upvar $spArrayName stereoPartner
    upvar $puArrayName peakUsed

    if {($toFromPartner($saName) != "NONE") && ([info exists curTargetForSA($toFromPartner($saName))])} {
	
	set tfPartnerProt  [lindex $curTargetForSA($toFromPartner($saName)) 2]
	set tfPartnerHeavy [lindex $curTargetForSA($toFromPartner($saName)) 3]

    } else {
	set tfPartnerProt "NONE"
	set tfPartnerHeavy "NONE"
    }
	
    if {($geminalPartner($saName) != "NONE") && ([info exists curTargetForSA($geminalPartner($saName))])} {

	set gPartnerHeavy [lindex $curTargetForSA($geminalPartner($saName)) 3]

    } else {
	set gPartnerHeavy "NONE"
    }

    if {($stereoPartner($saName) != "NONE") && ([info exists curTargetForSA($stereoPartner($saName))])} {

	set stereoPartnerProt  [lindex $curTargetForSA($stereoPartner($saName)) 2]
	
    } else {
	set stereoPartnerProt  "NONE"
    }

    set curTargetForSA($saName) [list "NONE" "NONE" "NONE" "NONE" "NONE" "NONE" "NONE"]


    foreach possTarget [shuffledList $possibleTargetsForSA($saName)] {
	set curProt  [lindex $possTarget 2]
	set curHeavy [lindex $possTarget 3]

	set paPartner [lindex $possTarget 4]
	if {[info exists curTargetForSA($paPartner)]} {
	    set paPartnerProt  [lindex $curTargetForSA($paPartner) 2]
	    set paPartnerHeavy [lindex $curTargetForSA($paPartner) 3]
	} else {
	    set paPartnerProt "NONE"
	    set paPartnerHeavy "NONE"
	}
		
	if {[insideTolerance $curProt $tfPartnerProt $protonTol] &&
	    [insideTolerance $curHeavy $tfPartnerHeavy $heavyTol] &&
	    [insideTolerance $curHeavy $gPartnerHeavy $heavyTol] &&
	    [insideTolerance [lindex $possTarget 5] $paPartnerProt $protonTol] &&
	    [insideTolerance [lindex $possTarget 6] $paPartnerHeavy $heavyTol] &&
	    [outsideTolerance $curProt $stereoPartnerProt $protonTol]} {
	    
	    if {! $peakUsed([lindex $possTarget 0])} {
		set curTargetForSA($saName) $possTarget
		set peakUsed([lindex $possTarget 0]) 1
		break
	    }
	}
    }
    
    return ""
}


proc insideTolerance {targVal myVal tol} {

    if {($targVal == "NONE") || ($myVal == "NONE")} {
	return 1
    } else {
	set delta [expr abs($targVal - $myVal)]
	if {$delta <= $tol} {
	    return 1
	} else {
	    return 0
	}
    }
}

proc outsideTolerance {targVal myVal tol} {

    if {($targVal == "NONE") || ($myVal == "NONE")} {
	return 1
    } else {
	set delta [expr abs($targVal - $myVal)]
	if {$delta >= $tol} {
	    return 1
	} else {
	    return 0
	}
    }
}

proc allTargetsCompatibleWithATarget {curSA curTarg ptArrayName ctArrayName ztArrayName protonTol heavyTol} {

    upvar $ctArrayName curTargetForSA 
    upvar $ptArrayName possibleTargetsForSA
    upvar $ztArrayName saHasZeroTargets


    #
    # Given a specific target for a specific SA, 
    # return a list of all the targets for that SA which 
    # are compatible with it:  ie., have shifts that 
    # agree with the given target's (within tolerance), 
    # and have peakAssign partner targets that are compatible
    #

    set retVal [list]

    if {[lindex $curTarg 0] != "NONE"} {

	set curProtTarg  [lindex $curTarg 2]
	set curHeavyTarg [lindex $curTarg 3]

	foreach possTarg $possibleTargetsForSA($curSA) {
		
	    #
	    # make sure this possible target agrees with the current target value
	    #
		
	    if {[insideTolerance $curProtTarg  [lindex $possTarg 2] $protonTol] &&
		[insideTolerance $curHeavyTarg [lindex $possTarg 3] $heavyTol]} {
		
		#
		# make sure this possible target agrees with its peakAssign partner's target value
		#
		
		set papSA   [lindex $possTarg 4]
		set papTarg $curTargetForSA($papSA)

		#
		# count this as compatible if the peakAssign partner SA has zero possible targets
		#

		if {$saHasZeroTargets($papSA)} {

		    lappend retVal $possTarg

		} else {
		
		    #
		    # if the peakAssign partner SA has no target, don't count this possible target as compatible
		    #

		    if {[lindex $papTarg 0] != "NONE"} {
			
			set papProtTarg  [lindex $papTarg 2]
			set papHeavyTarg [lindex $papTarg 3]
			
			if {[insideTolerance [lindex $possTarg 5] $papProtTarg  $protonTol] &&
			    [insideTolerance [lindex $possTarg 6] $papHeavyTarg $heavyTol]} {
			    
			    lappend retVal $possTarg
			}
		    }
		}
	    }
	}
    }
		
    return $retVal
}




proc recordAllTargetsCompatibleWithCurrentChoices {ptArrayName ctArrayName protonTol heavyTol compatTArrayName salVarName ztArrayName} {

    upvar $ptArrayName      possibleTargetsForSA 
    upvar $ctArrayName      curTargetForSA 
    upvar $compatTArrayName compatibleTargetsForSA
    upvar $salVarName       saList
    upvar $ztArrayName      saHasZeroTargets 

    array unset compatibleTargetsForSA

    foreach curSA $saList {

	set compatibleTargetsForSA($curSA) [allTargetsCompatibleWithATarget $curSA $curTargetForSA($curSA) possibleTargetsForSA curTargetForSA saHasZeroTargets $protonTol $heavyTol]

    }
    
    return ""
}



proc updateShiftAssignmentValuesFromTargets {pot salVarName ctArrayName tsArrayName meanCorrections} {

    upvar $ctArrayName curTargetForSA 
    upvar $tsArrayName targetScores
    upvar $salVarName  saList


    foreach curSAname $saList {

	if {! [$pot hasShiftAssignmentNamed $curSAname]} {
	    continue
	}

	set curSA [$pot shiftAssignmentNamed $curSAname]
	ShiftAssignment -this $curSA 

	if {$targetScores($curSAname) == 0} {
	    if {[$curSA hasProtonShift] && [$curSA isFrom] && ([lindex $meanCorrections 0] != "NONE")} {
		set origShift [$curSA protonShift]
		set newShift [expr $origShift + [lindex $meanCorrections 0]]
		$curSA setProtonShift $newShift
		$curSA appendToNote [format "Updated proton shift %s --> %s" $origShift $newShift]
		$curSA appendToNote [format "   Estimated correct proton shift from mean correction of all from proton shifts:  %.3f +/- %.3f ppm" [lindex $meanCorrections 0] [lindex $meanCorrections 1]]
	    }

	    if {[$curSA hasHeavyatomShift] && [$curSA isFrom] && ([lindex $meanCorrections 2] != "NONE")} {
		set origShift [$curSA heavyatomShift]
		set newShift [expr $origShift + [lindex $meanCorrections 2]]
		$curSA setHeavyatomShift $newShift
		$curSA appendToNote [format "Updated heavyatom shift %s --> %s" $origShift $newShift]
		$curSA appendToNote [format "   Estimated correct heavyatom shift from mean correction of all from heavyatom shifts:  %.3f +/- %.3f ppm" [lindex $meanCorrections 2] [lindex $meanCorrections 3]]
	    }

	    if {[$curSA hasProtonShift] && [$curSA isTo] && ([lindex $meanCorrections 4] != "NONE")} {
		set origShift [$curSA protonShift]
		set newShift [expr $origShift + [lindex $meanCorrections 4]]
		$curSA setProtonShift $newShift
		$curSA appendToNote [format "Updated proton shift %s --> %s" $origShift $newShift]
		$curSA appendToNote [format "   Estimated correct proton shift from mean correction of all to proton shifts:  %.3f +/- %.3f ppm" [lindex $meanCorrections 4] [lindex $meanCorrections 5]]
	    }

	    if {[$curSA hasHeavyatomShift] && [$curSA isTo] && ([lindex $meanCorrections 6] != "NONE")} {
		set origShift [$curSA heavyatomShift]
		set newShift [expr $origShift + [lindex $meanCorrections 6]]
		$curSA setHeavyatomShift $newShift
		$curSA appendToNote [format "Updated heavyatom shift %s --> %s" $origShift $newShift]
		$curSA appendToNote [format "   Estimated correct heavyatom shift from mean correction of all to heavyatom shifts:  %.3f +/- %.3f ppm" [lindex $meanCorrections 6] [lindex $meanCorrections 7]]
	    }

	    rename $curSA ""
	    continue
	}


	if {[$curSA hasProtonShift]} {
	    set origProtShift [$curSA protonShift]
	} else {
	    set origProtShift "NONE"
	}

	if {[$curSA hasHeavyatomShift]} {
	    set origHeavyShift [$curSA heavyatomShift]
	} else {
	    set origHeavyShift "NONE"
	}

	set protTarg  [lindex $curTargetForSA($curSAname) 2]
	set heavyTarg [lindex $curTargetForSA($curSAname) 3]

	set curUpdated 0

	if {$protTarg != "NONE"} {
	    $curSA setProtonShift $protTarg
	    set curUpdated 1
	}

	if {$heavyTarg != "NONE"} {
	    $curSA setHeavyatomShift $heavyTarg
	    set curUpdated 1
	} 

	if {$curUpdated} {

	    set rpt ""
	    appendLineToString rpt [format "Updated chemical shift value(s) via stripe filter:"]
	    appendLineToString rpt [format "   proton shift %s  --> %s" $origProtShift $protTarg]
	    appendLineToString rpt [format "   heavyatom shift %s  --> %s" $origHeavyShift $heavyTarg]
	    appendLineToString rpt [format "   Values came from peak %s assigned with %s" \
					[lindex $curTargetForSA($curSAname) 0] \
					[lindex $curTargetForSA($curSAname) end]]
	    appendLineToString rpt [format "   Score for this target was %.2f%%" \
					[expr $targetScores($curSAname) * 100]]

	    $curSA appendToNote $rpt
	}

	rename $curSA ""
    }
    return ""
}


proc reportTargetScores {salVarName ztArrayName tsArrayName overallScore ctArrayName pot mcVarName protonTol heavyatomTol} {

    upvar $ztArrayName saHasZeroTargets
    upvar $tsArrayName targetScores
    upvar $mcVarName   meanCorrections
    upvar $salVarName  saList
    upvar $ctArrayName curTargetForSA

    set retVal ""

    appendLineToString retVal [format "Stripe filter, applied with proton tolerance of %f ppm and heavyatom tolerance of %f ppm" $protonTol $heavyatomTol]

    #
    # report overall score & histogram of individual SA targets' scores
    #

    set numZeroTargets 0
    set numZeroScores 0
    set totNumSAs 0
    set scoreList [list]

    appendLineToString retVal [format "Stripe filter overall score is %f" $overallScore]

    foreach saName $saList {

	incr totNumSAs

	if {$saHasZeroTargets($saName) == 1} {
	    incr numZeroTargets
	} else {
	    if {$targetScores($saName) == 0} {
		incr numZeroScores
	    }
	    lappend scoreList $targetScores($saName)
	}
    }

    appendLineToString retVal [format "%d of %d (%.2f%%) SAs had no possible targets" \
				   $numZeroTargets $totNumSAs [expr 100 * $numZeroTargets / double($totNumSAs)]]
    
    appendLineToString retVal [format "%d of %d (%.2f%%) SAs had possible targets but zero scores" \
				   $numZeroScores $totNumSAs [expr 100 * $numZeroScores / double($totNumSAs)]]
    
    appendLineToString retVal [histogram -data $scoreList -title "Scores of SAs with possible targets"]
    
    #
    # generate histogram of differences between target shifts and starting shifts, 
    # ignoring targets with zero scores
    #

    set fromProtCorrections [list]
    set fromHeavyCorrections [list]
    set toProtCorrections [list]
    set toHeavyCorrections [list]

    foreach saName $saList {
	if {$targetScores($saName) == 0} {
	    continue
	}
	if {! [$pot hasShiftAssignmentNamed $saName]} {
	    continue
	}
	set curSA [$pot shiftAssignmentNamed $saName]
	ShiftAssignment -this $curSA
	set curProtTarget  [lindex $curTargetForSA($saName) 2]
	set curHeavyTarget [lindex $curTargetForSA($saName) 3]
	if {[$curSA hasProtonShift] && ($curProtTarget != "NONE")} {
	    set curCorrection [expr $curProtTarget - [$curSA protonShift]]
	    if {[$curSA isFrom]} {
		lappend fromProtCorrections $curCorrection
	    } else {
		lappend toProtCorrections $curCorrection
	    }
	}
			 
	if {[$curSA hasHeavyatomShift] && ($curHeavyTarget != "NONE")} {
	    set curCorrection [expr $curHeavyTarget - [$curSA heavyatomShift]]
	    if {[$curSA isFrom]} {
		lappend fromHeavyCorrections $curCorrection
	    } else {
		lappend toHeavyCorrections $curCorrection
	    }
	}
	rename $curSA ""
    }

    set meanCorrections [list "NONE" "NONE" "NONE" "NONE" "NONE" "NONE" "NONE" "NONE"]

    if {[llength $fromProtCorrections] > 0} {
	appendLineToString retVal [histogram -data $fromProtCorrections  -title "Differences between starting and target from proton shifts"]
	set curMean [mean $fromProtCorrections]
	set curDev  [standardDeviation $fromProtCorrections]
	appendLineToString retVal [format "Mean correction along from proton dimension: %.3f +/- %.3f ppm\n" $curMean $curDev]
	set meanCorrections [lreplace $meanCorrections 0 1 $curMean $curDev]
    } 

    if {[llength $fromHeavyCorrections] > 0} {
	appendLineToString retVal [histogram -data $fromHeavyCorrections -title "Differences between starting and target from heavyatom shifts"]
	set curMean [mean $fromHeavyCorrections]
	set curDev  [standardDeviation $fromHeavyCorrections]
	appendLineToString retVal [format "Mean correction along from heavyatom dimension: %.3f +/- %.3f ppm\n" $curMean $curDev]
	set meanCorrections [lreplace $meanCorrections 2 3 $curMean $curDev]
    } 
    
    if {[llength $toProtCorrections] > 0} {
	appendLineToString retVal [histogram -data $toProtCorrections  -title "Differences between starting and target to proton shifts"]
	set curMean [mean $toProtCorrections]
	set curDev  [standardDeviation $toProtCorrections]
	appendLineToString retVal [format "Mean correction along to proton dimension: %.3f +/- %.3f ppm\n" $curMean $curDev]
	set meanCorrections [lreplace $meanCorrections 4 5 $curMean $curDev]

    } 

    if {[llength $toHeavyCorrections] > 0} {
	appendLineToString retVal [histogram -data $toHeavyCorrections -title "Differences between starting and target to heavyatom shifts"]
	set curMean [mean $toHeavyCorrections]
	set curDev  [standardDeviation $toHeavyCorrections]
	appendLineToString retVal [format "Mean correction along to heavyatom dimension: %.3f +/- %.3f ppm\n" $curMean $curDev]
	set meanCorrections [lreplace $meanCorrections 6 7 $curMean $curDev]
    } 
    
    return $retVal
}

    

proc calcNewTargetScores {ptArrayName ctArrayName tsArrayName osVarName salVarName protonTol heavyatomTol} {

    upvar $ptArrayName rawTargetsForSA 
    upvar $ctArrayName curTargetForSA 
    upvar $tsArrayName targetScores
    upvar $osVarName   overallScore 
    upvar $salVarName  saList 

    #
    # report fraction of possible targets each SA's current target is compatible with, along with 
    # the overall fraction.  Be sure to prevent double counting of targets with the same SA or the 
    # same peak
    #

    set totCompatTargets 0
    set totTargets 0

    foreach saName $saList {

	#
	# count how many of its possible targets are compatible with its chosen target
	#

	set curTarg $curTargetForSA($saName)

	set compatTargs [list]

	if {[lindex $curTarg 0] != "NONE"} {

	    set curProtTarget  [lindex $curTarg 2]
	    set curHeavyTarget [lindex $curTarg 3]       

	    foreach otherTarg $rawTargetsForSA($saName) {
		set otherProtTarget  [lindex $otherTarg 2]
		set otherHeavyTarget [lindex $otherTarg 3]
		
		if {[insideTolerance $curProtTarget $otherProtTarget $protonTol] &&
		    [insideTolerance $curHeavyTarget $otherHeavyTarget $heavyatomTol]} {
		    lappend compatTargs $otherTarg
		}
	    }
	}
	
	set saPartners [list]
	set peakPartners [list]
	foreach elem $compatTargs {
	    lappend saPartners [lindex $elem 4]
	    lappend peakPartners [lindex $elem 0]
	}
	
	set numSApartners [llength [lsort -unique $saPartners]]
	set numPeakPartners [llength [lsort -unique $peakPartners]]
	set numCompatTargets [min $numSApartners $numPeakPartners]
	
	set rawSAPartners [list]
	set rawPeakPartners [list]

	foreach elem $rawTargetsForSA($saName) {
	    lappend rawSAPartners [lindex $elem 4]
	    lappend rawPeakPartners [lindex $elem 0]
	}

	set numRawSAPartners [llength [lsort -unique $rawSAPartners]]
	set numRawPeakPartners [llength [lsort -unique $rawPeakPartners]]
	set numRawTargets [min $numRawSAPartners $numRawPeakPartners]

	
	if {$numRawTargets > 0} {
	    set targetScores($saName) [expr $numCompatTargets / double($numRawTargets)]
	} else {
	    # or NONE?
	    set targetScores($saName) 1   
	}
	
	set totCompatTargets [expr $totCompatTargets + $numCompatTargets]
	set totTargets [expr $totTargets + $numRawTargets]
    }


    if {$totTargets > 0} {
	set overallScore [expr $totCompatTargets / double($totTargets)]
    } else {
	set overallScore 1
    }

    return ""
}

	


proc calcTargetScores {ctArrayName epArrayName tsArrayName osVarName curTargArrayName plVarName salVarName sitArrayName} {

    upvar $ctArrayName      compatibleTargetsForSA
    upvar $epArrayName      expectedPeaksForSA
    upvar $tsArrayName      targetScores
    upvar $osVarName        overallScore
    upvar $curTargArrayName curTargetForSA
    upvar $plVarName        peakList
    upvar $salVarName       saList
    upvar $sitArrayName     saIsTo

    #
    # create arrays to record the number of times a particular peak 
    # is used as part of the current set of compatible targets for from SAs
    # and the current set of compatible targets for to SAs
    #

    array unset nTimesPeakUsedForFrom
    array unset nTimesPeakUsedForTo

    foreach peakName $peakList {
	set nTimesPeakUsedForFrom($peakName) 0
	set nTimesPeakUsedForTo($peakName) 0
    }


    #
    # create an array to record the number of times a particular pair
    # of SAs is used as part of the current set of compatible targets,
    # treating fromA-toB and toB-fromA as separate entities for now
    #

    array unset nTimesSAPairUsed

    foreach saA $saList {
	foreach saB $saList {
	    if {$saIsTo($saA) != $saIsTo($saB)} {
		set curSApair [join [list $saA $saB]]
		set nTimesSAPairUsed($curSApair) 0
	    }
	}
    }


    #
    # for each compatible target, record the number of times each 
    # peak is used, and the number of times each SA pair is used
    #

    foreach curSA $saList {
	 
	foreach compatTarg $compatibleTargetsForSA($curSA) {
	    
	    set curPeak [lindex $compatTarg 0]
	    if {$saIsTo($curSA)} {
		incr nTimesPeakUsedForTo($curPeak)
	    } else {
		incr nTimesPeakUsedForFrom($curPeak)
	    }
	    
	    set curSApair [lindex $compatTarg end]
	    incr nTimesSAPairUsed($curSApair) 
	    
	}
    }
    
    #
    # score for each SA is the minimum of its scaled number of peaks
    # and its scaled number of SA pairs, divided by the number of intrares partner SAs
    # (which is equal to both the number of expected peaks and the number of expected SA pairs)
    #
    # Overall score is the min of the total scaled number of peaks and 
    # the total scaled number of SA pairs in the current set of compatible 
    # targets
    #

    array unset targetScores

    set overallNumPeaks 0
    set overallNumSAPairs 0
    set overallNumExpected 0

    foreach curSA $saList {

	set scaledNumPeaks   0
	set scaledNumSAPairs 0

	foreach compatTarg $compatibleTargetsForSA($curSA) {
	    set curPeak [lindex $compatTarg 0]
	    if {$saIsTo($curSA)} {
		set scaledNumPeaks [expr $scaledNumPeaks + (1.0 / double($nTimesPeakUsedForTo($curPeak)))]
	    } else {
		set scaledNumPeaks [expr $scaledNumPeaks + (1.0 / double($nTimesPeakUsedForFrom($curPeak)))]
	    }

	    set curSApair [lindex $compatTarg end]
	    set scaledNumSAPairs [expr $scaledNumSAPairs + (1.0 / double($nTimesSAPairUsed($curSApair)))]
	}

	set curNumExpected [llength $expectedPeaksForSA($curSA)]
	if {$curNumExpected == 0} {
	    set curPeakScore 0
	    set curSAPairScore 0
	} else {

	    set curPeakScore   [expr $scaledNumPeaks / double($curNumExpected)]
	    set curSAPairScore [expr $scaledNumSAPairs / double($curNumExpected)]
	}


	set targetScores($curSA) [min $curPeakScore $curSAPairScore]

	set overallNumPeaks    [expr $overallNumPeaks + $scaledNumPeaks]
	set overallNumSAPairs  [expr $overallNumSAPairs + $scaledNumSAPairs]
	set overallNumExpected [expr $overallNumExpected + $curNumExpected]
    }

    #
    # Overall num expected is twice as large as the actual expected number of 
    # peaks because a peak should be visible for both its from and to SA targets.  
    #
    # But overall num expected is correct for the number of SA pairs, since 
    # I'm counting fromA-toB and toB-fromA as separate entities.
    #

    set overallPeakScore [expr $overallNumPeaks / (0.5 * $overallNumExpected)]
    set overallSAPairScore [expr $overallNumSAPairs / $overallNumExpected]
    set overallScore [min $overallPeakScore $overallSAPairScore]
	
    return ""
}


proc copyTargetsAndScores {intoTargArrName fromTargArrName intoScoreArrName fromScoreArrName salVarName} {

    upvar $salVarName   saList 
    upvar $intoTargArrName  intoTargetArray
    upvar $fromTargArrName  fromTargetArray
    upvar $intoScoreArrName intoScoreArray
    upvar $fromScoreArrName fromScoreArray
    array unset intoArray

    foreach curSA $saList {
	
	set intoTargetArray($curSA) $fromTargetArray($curSA)
	set  intoScoreArray($curSA)  $fromScoreArray($curSA)
    }

    return ""
}

proc resetSomeTargets {salVarName ctArrayName tfArrayName gpArrayName spArrayName tsArrayName charScore} {

    upvar $ctArrayName  curTargetForSA
    upvar $tfArrayName  toFromPartner
    upvar $gpArrayName  geminalPartner
    upvar $spArrayName  stereoPartner
    upvar $salVarName   saList
    upvar $tsArrayName  targetScores

    #
    # Decide to reset a shiftAssignment's target based on its score and the characteristic score.
    # 
    # If a shiftAssignment's target is reset, also reset its geminal partner's, 
    # its stereo partner's, and all of their to-from partners' 
    #

    foreach curSA $saList {

	if {[info exists targetScores($curSA)]} {
	    set curScore $targetScores($curSA)
	} else {
	    set curScore 0.0
	}

	#
	# characteristic score is the score at which the probability of reset has dropped to 1/e
	# 
	# If charScore is zero, we only reset SAs with zero scores
	#

	if {$charScore == 0} {
	    if {$curScore == 0} {
		set resetProb 1
	    } else {
		set resetProb 0
	    }
	} else {
	    set resetProb [expr exp(-(pow($curScore,2)) / double(pow($charScore,2)))]
	}

	#
	# if we've decided to reset this SA, 
	#

	if {[uniformRandom] <= $resetProb} {
	    
	    #
	    # erase its own target entry, 
	    # its geminal partner's target entry,
	    # its stereo partner's target entry,
	    # and its peak assign partner's target entry
	    #
	    
	    if {[info exists curTargetForSA($curSA)]} {
		set paPartner [lindex $curTargetForSA($curSA) 4]
	    } else {
		set paPartner "NONE"
	    }

	    catch {unset curTargetForSA($curSA)}		    
	    catch {unset curTargetForSA($geminalPartner($curSA))}
	    catch {unset curTargetForSA($stereoPartner($curSA))}
	    catch {unset curTargetForSA($paPartner)}

	    #
	    # erase all of their to-from partners' target entries as well
	    #

	    catch {unset curTargetForSA($toFromPartner($curSA))}		    
	    catch {unset curTargetForSA($toFromPartner($geminalPartner($curSA)))}
	    catch {unset curTargetForSA($toFromPartner($stereoPartner($curSA)))}
	    catch {unset curTargetForSA($toFromPartner($paPartner))}
	}
    }

    return ""
}


proc resetAllTargetsAndScores {ctArrayName tsArrayName salVarName} {

    upvar $ctArrayName  curTargetForSA
    upvar $tsArrayName  targetScores
    upvar $salVarName   saList

    array unset curTargetForSA

    foreach curSA $saList {
	set curTargetForSA($curSA) [list "NONE" "NONE" "NONE" "NONE" "NONE" "NONE" "NONE"]
	set targetScores($curSA) 0
    }

    return ""
}


proc acceptMCmove {curScore prevScore charDeltaScore} {

    if {$curScore > $prevScore} {
	return 1
    } elseif {$charDeltaScore == 0} {
	return 0
    } else {
	set delta [expr $prevScore - $curScore]
	set acceptProb [expr exp (-(pow($delta,2)) / double(pow($charDeltaScore,2)))]
	if {[uniformRandom] < $acceptProb} {
	    return 1
	} else {
	    return 0
	}
    }
}


proc removeSAsWithNoRawTargets {pot ztArrayName ipArrayName salVarName} {

    upvar $ztArrayName  saHasZeroTargets
    upvar $salVarName   saList
    upvar $ipArrayName  intraPartnerSAs 

    set retVal [list]

    foreach saName $saList {
	if {$saHasZeroTargets($saName) && ([llength $intraPartnerSAs($saName)] > 0)} {
	    $pot removeShiftAssignmentNamed $saName
	    lappend retVal $saName
	}
    }

    return $retVal
}


proc stripeCorrection args {

    set pot                            [requiredFlagVal $args -pot]
    set protonTol                      [flagVal $args [list -protonTolerance -stripeGenerationProtonTolerance] 0.02]
    set heavyatomTol                   [flagVal $args [list -heavyatomTolerance stripeGenerationHeavyatomTolerance] 0.2]
    set remVar                         [flagVal $args -remarksVariableName ""]
    set ignoreBackboneSeq              [flagExists $args -ignoreBackboneSequential]
    set numMCsteps                     [flagVal $args -numMonteCarloSteps 10]
    set numRuns                        [flagVal $args -numRuns 20]
    set charScore                      [flagVal $args -characteristicScore 0.3]
    set charDeltaScore                 [flagVal $args -characteristicDeltaScore 0.005]
    set maxMoveTries                   [flagVal $args -maxMoveTries 10]
    set remVar                         [flagVal $args -remarksVariableName ""]
    set applyPrefilter                 [flagExists $args -applyPrefilter]
    set oldScore                       [flagExists $args -oldScore]
    set dropSAsWithoutRawTargets       [flagExists $args -dropSAsWithoutRawTargets]


    #
    # Correct the chemical shift values in the shiftAssignments to agree with the 
    # unfolded NOESY peak positions
    #
    # 1.  For each SA, record a list of the other SAs w/ opposite to-from sense, 
    #     not including the to-from partner, and in the same residue (or, optionally, 
    #     bbn-seq to the SA of interest).
    #
    # 2.  Record all PAs (generated with a broad-tolerance match) that correspond to 
    #     these known-close SA pairs, giving a set of all possible target PAs.
    #
    # 3.  For each SA, choose one of its possible target PAs to be the current target
    #
    # 4.  Optimize those choices with Monte Carlo
    #


    # TODC--change loop to be a real MC. 
    #
    # If an MC iteration gives a lower overall score, revert to the last set of choices and try again.  
    # Accept iterations with higher scores part of the time.  
    # 

    #
    # record the names of each shift assignment
    #

    updateUser "Recording names of each shiftAssignment \r"
    
    recordShiftAssignmentNames $pot saList

    #
    # for each shiftAssignment, record its opposite sense partner name
    #

    updateUser "Recording each shiftAssignment's from-to partner \r"

    recordToFromPartners $pot toFromPartner

    #
    # for each shiftAssignment, record whether it's to or from
    #

    updateUser "Recording each shiftAssignment's to-from sense \r"

    recordToFromSense $pot saIsTo


    #
    # for each shiftAssignment, record all shiftAssignments that are in the same residue 
    # but have opposite to-from tags
    #

    updateUser "Gathering lists of intraresidue partner ShiftAssignments\r"

    recordIntraresPartners $pot intraPartnerSAs



    #
    # if we're including backbone sequential PAs when generating targets,
    # add them to the intraPartnerSAs arrays
    #
    
    if {! $ignoreBackboneSeq} {

	updateUser "Gathering lists of sequential backbone partner ShiftAssignments\r"

	recordBbnSeqPartners $pot intraPartnerSAs

    }

    #
    # filter out to-from partner SAs from the lists of intraresidue partner SAs
    #

    removeToFromPartnersFromIntraPartners intraPartnerSAs toFromPartner saList


    #
    # for every SA with a geminal partner (ie., its heavyatom selection is identical to its partner),
    # record the partner's name if they both have heavyatom shifts
    #
    # Note that this allows only one geminal partner per SA
    #

    updateUser "Recording geminal partners among the shiftAssignments\r"

    recordGeminalPartners $pot geminalPartner

    #
    # record the name of each SA's stereopartner.  Use this to force chosen targets to be 
    # sufficiently different for stereopartners' proton shifts.  Helps to ensure that the 
    # same stripe isn't used for both stereopartners, as can often happen if they're geminal 
    # partners as well.  
    #
    # Reason for using stereopartners specifically is that this works in cases like the 
    # to dimension of a 3dN, where there are no geminal partners because of the lack of 
    # heavy atoms, but we still want to avoid re-using the same stripe
    #

    updateUser "Recording stereo partners among the shiftAssignments \r"

    recordStereoPartners $pot stereoPartner

    #
    # for each PeakAssignment,
    # record Peak, PeakAssign, from and to shiftAssign names, and unfolded peak positions in a list,
    # with an array index based on the names of its from and to shiftAssigns
    #

    updateUser "Recording each peak assignment's unfolded position\r"

    recordAllPossibleTargetPAs $pot possibleTargetsForSA intraPartnerSAs peakList saList

    #
    # record peaks with no possible targets at all
    #

    recordSAsWithNoTargets possibleTargetsForSA  saHasZeroRawTargets saList 
    
    #
    # remove possible targets that don't have expected compatible partners
    #

    if {$applyPrefilter} {
	filterTargetsForInitialReasonability possibleTargetsForSA toFromPartner geminalPartner stereoPartner saHasZeroRawTargets saList $protonTol $heavyatomTol
    }
    
    #
    # record peaks with no reasonable targets at all
    #
    
    recordSAsWithNoTargets possibleTargetsForSA  saHasZeroFilteredTargets saList 
    
    #
    # calc max number of intrares peaks I could see, for use in calculating score, below
    #
    
    calcExpectedPeaks intraPartnerSAs expectedPeaksForSA saHasZeroFilteredTargets saList

    set bestOverallScore -1

    #
    # make several independent attempts at target choice
    #

    for {set runCount 0} {$runCount < $numRuns} {incr runCount} {

	
	resetAllTargetsAndScores curTargetForSA targetScores saList
       
	copyTargetsAndScores prevTargets curTargetForSA prevScores targetScores saList
	set prevScore -1
	

	#
	# begin Monte Carlo optimization of target choices
	#

	for {set mcCount 0} {$mcCount < $numMCsteps} {incr mcCount} {

	    set moveAccepted 0
	    set numMoveTries 0

	    while {! $moveAccepted} {

		incr numMoveTries

		updateUser [format "New stripe run %d MC step %d try %d curScore %f  bestScore %f \r" $runCount $mcCount $numMoveTries $prevScore $bestOverallScore]

		#
		# begin MC step at previously-accepted targets  Need to copy previous targets' scores as well, so resetSomeTargets will work
		#

		copyTargetsAndScores curTargetForSA prevTargets targetScores prevScores saList

		#
		# select some targets and remove them
		#
		
		resetSomeTargets saList curTargetForSA toFromPartner geminalPartner stereoPartner targetScores $charScore

		#
		# generate targets for all SAs that are missing them
		#
		
		chooseTargets possibleTargetsForSA curTargetForSA toFromPartner geminalPartner stereoPartner saIsTo peakList saList $protonTol $heavyatomTol
		
		#
		# record all possible targets that are within tolerance of the current set of chosen targets,
		# including both the SA's own target value and its peak assign partner's target value
		#		
		
		recordAllTargetsCompatibleWithCurrentChoices possibleTargetsForSA curTargetForSA $protonTol $heavyatomTol compatibleTargetsForSA saList saHasZeroFilteredTargets
		
		#
		# count the number of unique peaks covered by the possible targets compatible with the current target choices
		#
		
		if {$oldScore} {
		    calcTargetScores compatibleTargetsForSA expectedPeaksForSA  targetScores overallScore curTargetForSA peakList saList saIsTo
		} else {
		    calcNewTargetScores possibleTargetsForSA curTargetForSA targetScores overallScore saList $protonTol $heavyatomTol
		}
		
		#
		# see if I accept this move
		#
		
		if {[acceptMCmove $overallScore $prevScore $charDeltaScore]} {
		    
		    set prevScore $overallScore
		    copyTargetsAndScores prevTargets curTargetForSA prevScores targetScores saList
		    set moveAccepted 1
		    
		} 
		
		if {$numMoveTries > $maxMoveTries} {
		    set moveAccepted 1
		}
		
		#
		# save these targets if they're the best so far
		#
		
		if {$overallScore > $bestOverallScore} {
		    set bestOverallScore $overallScore
		    copyTargetsAndScores bestTargets curTargetForSA bestScores targetScores saList
		}
		
	    }
	}
    }

    if {$dropSAsWithoutRawTargets} {
	set removedSAs [removeSAsWithNoRawTargets $pot saHasZeroRawTargets intraPartnerSAs saList]
    }


    #
    # report what I did
    #
    
    if {$remVar != ""} {
	
	upvar $remVar tempRem
    }

    set r [reportTargetScores saList saHasZeroRawTargets bestScores $bestOverallScore bestTargets $pot meanCorrections $protonTol $heavyatomTol]
    
    if {$dropSAsWithoutRawTargets} {
	appendLineToString r [format "Removed %d shift assignments that had expected intraresidue partners but no target shifts.  They are:" [llength $removedSAs]]
	
	foreach elem $removedSAs {
	    appendLineToString r $elem
	}
    }
    
    lappend tempRem $r

    #
    # update the shift assignments' values for the shifts, based on the best
    # targets
    #

    updateShiftAssignmentValuesFromTargets $pot saList bestTargets \
	bestScores $meanCorrections
	
    

    return ""
}

}

namespace import StripeFilter::*
