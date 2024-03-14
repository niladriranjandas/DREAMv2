package provide marvin 1.0

#
# tools for dealing with atom selections
#


#
# whether two selections refer to atoms in the same residue
#

proc selectionsInSameResidue {sel1 sel2} {

    set resList1 [residuesInSelection $sel1]
    set resList2 [residuesInSelection $sel2]

    foreach elem1 $resList1 {
	if {[lsearch $resList2 $elem1] != -1} {
	    return 1
	}
    }

    return 0
}

# 
# Returns a list of {segName resNum} pairs, with no repeats, 
# for atoms in the given selection
#

proc residuesInSelection {tempSel} {

    set resList [list]
    
    for {set x 0} {$x < [$tempSel size]} {incr x} {

	set curAtom [$tempSel atomAtPos $x]
	Atom -this $curAtom

	set curSeg     [string trim [$curAtom segmentName]]
	set curResNum  [string trim [$curAtom residueNum]]
	set curResName [string trim [$curAtom residueName]]

	if {$curSeg != ""} {
	    lappend resList [join [list $curSeg $curResNum $curResName]]
	} else {
	    lappend resList [join [list $curResNum $curResName]]
	}

	rename $curAtom ""
    }
    
    return [lsort -dictionary -unique $resList]
}

#
# Given a residue code like residuesInSelection produces, 
# or a list in the form of {segname, resnum}, or a resnum,
# return a residue code for the previous residue
#

proc previousResidue args {

    set curSeg ""
    set curResName ""
    set curResNum ""
    
    # called with either just an integer resNum, or
    # a joined-list residue code
    
    if {[llength $args] == 1} {
	set temp [split [lindex $args 0]]

	switch [llength $temp] {
	    1 { set curResNum [lindex $temp 0] }
	    2 { set curResNum [lindex $temp 0]
		set curResName [lindex $temp 1] }
	    3 {
		set curSeg [lindex $temp 0]
		set curResNum [lindex $temp 1]
		set curResName [lindex $temp 2]
	    }
	    default { error "can't parse $args" }
	}
    } else {
	set curSeg [lindex $args 0]
	set curResNum [lindex $args 1]
    }

    if {$curSeg == ""} {
	set tempSel [AtomSel -args [format "(resid %d)" [expr $curResNum - 1]]]
	set retVal [residuesInSelection $tempSel]
	rename $tempSel ""
    } else {
	set tempSel [AtomSel -args [format "(segid %s and resid %d)" $curSeg [expr $curResNum - 1]]]
	set retVal [residuesInSelection $tempSel]
	rename $tempSel ""
    } 

    switch [llength $retVal] {
	0 { return [list]} 
	1 { return [lindex $retVal 0] }
	default { error [format "previousResidue of %s finds more than one residue: %s" $args $retVal]}
    }
}


proc nextResidue args {

    set curSeg ""
    set curResName ""
    set curResNum ""
    
    # called with either just an integer resNum, or
    # a joined-list residue code
    
    if {[llength $args] == 1} {
	set temp [split [lindex $args 0]]

	switch [llength $temp] {
	    1 { set curResNum [lindex $temp 0] }
	    2 { set curResNum [lindex $temp 0]
		set curResName [lindex $temp 1] }
	    3 {
		set curSeg [lindex $temp 0]
		set curResNum [lindex $temp 1]
		set curResName [lindex $temp 2]
	    }
	    default { error "can't parse $args" }
	}
    } else {
	set curSeg [lindex $args 0]
	set curResNum [lindex $args 1]
    }

    if {$curSeg == ""} {
	set tempSel [AtomSel -args [format "(resid %d)" [expr $curResNum + 1]]]
	set retVal [residuesInSelection $tempSel]
	rename $tempSel ""
    } else {
	set tempSel [AtomSel -args [format "(segid %s and resid %d)" $curSeg [expr $curResNum + 1]]]
	set retVal [residuesInSelection $tempSel]
	rename $tempSel ""
    } 

    switch [llength $retVal] {
	0 { return [list]} 
	1 { return [lindex $retVal 0] }
	default { error [format "nextResidue of %s finds more than one residue: %s" $args $retVal]}
    }
}

#
# how many different residues are in a particular selection?
# Used by setupRgyr and in determining NOE precision
#


proc numResiduesInSelection {tempSel} {

    return [llength [residuesInSelection $tempSel]]
}


#
# break up selections of the form ((subsel1) OR (subsel2) OR ...)
# into lists of the subselections
#

# 
# grab the next piece of a selection string that's surrounded by parens
#

proc nextSubSel args {

    set aSelString  [requiredFlagVal $args -selString]
    set startPos    [flagVal $args -startPos 0]

    #
    # eat leading text until we hit an open-paren
    #

    while {1} {

	if {[string index $aSelString $startPos] == "("} {
	    break
	}

	if {$startPos == [string length $aSelString]} {
	    error "No open paren found"
	} 

	incr startPos
    }

    set numOpenParens 1
    set endPos [expr $startPos + 1]

    #
    # keep eating characters until the count of
    # open parens reaches zero
    #

    while {($numOpenParens > 0)} {

	    
	if {[string index $aSelString $endPos] == "("} {
		
	    incr numOpenParens
	}

	if {[string index $aSelString $endPos] == ")"} {
		
	    incr numOpenParens -1
	}

	if {$endPos == [string length $aSelString]} {
	    error "String ends with un-matched open parens"
	}
	
	incr endPos
    }

    # 
    # undo the last increment to endPos
    #

    incr endPos -1
    
    if {[flagExists $args -positions]} {
	return [list $startPos $endPos]
    } else {
	return [string range $aSelString $startPos $endPos]
    }

}

proc splitSelection args {

    if {[flagExists $args -sel]} {

	set tempSel [flagVal $args -sel]
	AtomSel -this $tempSel
	set origString [$tempSel string]
	rename $tempSel ""

    } elseif {[flagExists $args -string]} {
	
	set origString [flagVal $args -string]

    } else {

	error "required flags -sel or -string not found"
    }

    set subSels [list]

    #
    # if this string has parens around the whole thing, 
    # get rid of them.
    #

    if {[nextSubSel -selString $origString] == $origString} {
	set origString [string range $origString 1 end-1]
    }

    #
    # eliminate any leading or trailing whitespace
    #

    set origString [string trim $origString]
    
    #
    # Make sure this string is of the form (subsel 1) OR (subsel2) OR ...
    # with NOTHING ELSE
    #

    set curPos [nextSubSel -selString $origString -startPos 0 -positions]
    if {[lindex $curPos 0] != 0} {
	error "selection string doesn't begin with a subselection"
    }
    
    while {[lindex $curPos 1] != [expr [string length $origString] - 1]} {
	
	set nextPos [nextSubSel -selString $origString -startPos [lindex $curPos 1] -positions]       

	set intervening [string range $origString [expr [lindex $curPos 1] + 1] [expr [lindex $nextPos 0] - 1]]

	if {[string toupper [string trim $intervening]] != "OR"} {
	    error "subselections not separated by OR"
	}

	lappend subSels [string range $origString [lindex $curPos 0] [lindex $curPos 1]]
	set curPos $nextPos
    }

    lappend subSels [string range $origString [lindex $curPos 0] [lindex $curPos 1]]

    return $subSels
}


#
# Given two selections, see whether any atoms in the first
# are bonded to any atoms in the second.
#

proc cacheBondArray {bondArrayName} {

    upvar $bondArrayName bondArray

    for {set count 0} {$count < [xplorSim numBonds]} {incr count} {
	set curBond [xplorSim bondPairByID $count]
	set curFrom [lindex $curBond 0]
	set curTo   [lindex $curBond 1]

	lappend bondArray($curFrom) $curTo
	lappend bondArray($curTo) $curFrom
    }
}

proc selectionsAreBonded {sel1 sel2 bondArrayName} {

    upvar $bondArrayName bondArray

    set l1 [$sel1 indices]
    set l2 [$sel2 indices]

    foreach i1 $l1 {
	set curBondPartners $bondArray($i1)
	foreach i2 $l2 {
	    if {[lsearch -exact $curBondPartners $i2] != -1} {
		return 1
	    }
	}
    }

    return 0
}




#
# stereospecific selection conversion
#

proc stereopartnerList {} {

    #
    # Just a convenient place to keep a list of all the stereopartners
    # Each entry gives residue type, stereospecific name, 
    # stereospecific partner name, and nonstereospecific name
    #

    # FIX:  I need to add stereopartners for DNA/RNA

    set stereopartners [list]

    lappend stereopartners [list gly ha1  ha2  ha*]
    lappend stereopartners [list gly ha2  ha1  ha*]

    lappend stereopartners [list val hg1* hg2* hg*]
    lappend stereopartners [list val hg2* hg1* hg*]
    lappend stereopartners [list val cg1  cg2  cg*]
    lappend stereopartners [list val cg2  cg1  cg*]

    lappend stereopartners [list leu hb1  hb2  hb*]
    lappend stereopartners [list leu hb2  hb1  hb*]
    lappend stereopartners [list leu hd1* hd2* hd*]
    lappend stereopartners [list leu hd2* hd1* hd*]
    lappend stereopartners [list leu cd1  cd2  cd*]
    lappend stereopartners [list leu cd2  cd1  cd*]

    lappend stereopartners [list ile hg1  hg2  hg*]
    lappend stereopartners [list ile hg2  hg1  hg*]

    lappend stereopartners [list phe hb1  hb2  hb*]
    lappend stereopartners [list phe hb2  hb1  hb*]
    lappend stereopartners [list phe hd1  hd2  hd*]
    lappend stereopartners [list phe hd2  hd1  hd*]
    lappend stereopartners [list phe he1  he2  he*]
    lappend stereopartners [list phe he2  he1  he*]
    lappend stereopartners [list phe cd1  cd2  cd*]
    lappend stereopartners [list phe cd2  cd1  cd*]
    lappend stereopartners [list phe ce1  ce2  ce*]
    lappend stereopartners [list phe ce2  ce1  ce*]

    lappend stereopartners [list tyr hb1  hb2  hb*]
    lappend stereopartners [list tyr hb2  hb1  hb*]
    lappend stereopartners [list tyr hd1  hd2  hd*]
    lappend stereopartners [list tyr hd2  hd1  hd*]
    lappend stereopartners [list tyr he1  he2  he*]
    lappend stereopartners [list tyr he2  he1  he*]
    lappend stereopartners [list tyr cd1  cd2  cd*]
    lappend stereopartners [list tyr cd2  cd1  cd*]
    lappend stereopartners [list tyr ce1  ce2  ce*]
    lappend stereopartners [list tyr ce2  ce1  ce*]

    lappend stereopartners [list trp hb1  hb2  hb*]
    lappend stereopartners [list trp hb2  hb1  hb*]

    lappend stereopartners [list his hb1  hb2  hb*]
    lappend stereopartners [list his hb2  hb1  hb*]

    lappend stereopartners [list ser hb1  hb2  hb*]
    lappend stereopartners [list ser hb2  hb1  hb*]

    lappend stereopartners [list cys hb1  hb2  hb*]
    lappend stereopartners [list cys hb2  hb1  hb*]

    lappend stereopartners [list met hb1  hb2  hb*]
    lappend stereopartners [list met hb2  hb1  hb*]
    lappend stereopartners [list met hg1  hg2  hg*]
    lappend stereopartners [list met hg2  hg1  hg*]

    lappend stereopartners [list asn hb1  hb2  hb*]
    lappend stereopartners [list asn hb2  hb1  hb*]
    lappend stereopartners [list asn hd21 hd22 hd2*]
    lappend stereopartners [list asn hd22 hd21 hd2*]

    lappend stereopartners [list gln hb1  hb2  hb*]
    lappend stereopartners [list gln hb2  hb1  hb*]
    lappend stereopartners [list gln hg1  hg2  hg*]
    lappend stereopartners [list gln hg2  hg1  hg*]
    lappend stereopartners [list gln he21 he22  he2*]
    lappend stereopartners [list gln he22 he21  he2*]

    lappend stereopartners [list asp hb1  hb2  hb*]
    lappend stereopartners [list asp hb2  hb1  hb*]

    lappend stereopartners [list glu hb1  hb2  hb*]
    lappend stereopartners [list glu hb2  hb1  hb*]
    lappend stereopartners [list glu hg1  hg2  hg*]
    lappend stereopartners [list glu hg2  hg1  hg*]

    lappend stereopartners [list lys hb1  hb2  hb*]
    lappend stereopartners [list lys hb2  hb1  hb*]
    lappend stereopartners [list lys hg1  hg2  hg*]
    lappend stereopartners [list lys hg2  hg1  hg*]
    lappend stereopartners [list lys hd1  hd2  hd*]
    lappend stereopartners [list lys hd2  hd1  hd*]
    lappend stereopartners [list lys he1  he2  he*]
    lappend stereopartners [list lys he2  he1  he*]

    lappend stereopartners [list arg hb1  hb2  hb*]
    lappend stereopartners [list arg hb2  hb1  hb*]
    lappend stereopartners [list arg hg1  hg2  hg*]
    lappend stereopartners [list arg hg2  hg1  hg*]
    lappend stereopartners [list arg hd1  hd2  hd*]
    lappend stereopartners [list arg hd2  hd1  hd*]

    lappend stereopartners [list pro hb1  hb2  hb*]
    lappend stereopartners [list pro hb2  hb1  hb*]
    lappend stereopartners [list pro hg1  hg2  hg*]
    lappend stereopartners [list pro hg2  hg1  hg*]
    lappend stereopartners [list pro hd1  hd2  hd*]
    lappend stereopartners [list pro hd2  hd1  hd*]

    return $stereopartners
}


#
# Converts, eg. (leu hd1*) to (leu hd*)
#

proc makeSelectionNonStereoSpecific {aSelString} {

    set aSel [AtomSel -args $aSelString]

    #
    # complain if this isn't a single-residue selection
    #

    switch [numResiduesInSelection $aSel] {
	0 { rename $aSel ""
	    error [format "makeSelectionNonStereoSpecific called on an empty selection: %s" $aSelString] } 
	1 {} 
	default { rename $aSel ""
	    error [format "makeSelectionNonStereoSpecific called on a selection that extends over >1 residue: %s. \n%s" \
		       $aSelString "Are the restraints' selections segid-specific?"] }
    }

    #
    # check if the original selection string includes a segid tag.  If not, 
    # don't create one in the output selection string
    #

    if {[string match -nocase "*segid*" $aSelString]} {
	set origHasSegid 1
    } else {
	set origHasSegid 0
    }

    #
    # grab the segid, resnum, and resname from one of the selected atoms
    #
    
    set curAtom [$aSel atomAtPos 0]
    Atom -this $curAtom

    set curSeg     [$curAtom segmentName]
    set curNum     [$curAtom residueNum]
    set curResType [$curAtom residueName]
    rename $curAtom ""

    set retVal $aSelString

    foreach stereopartner [stereopartnerList] {

	#
	# don't even bother checking further if the residue types aren't the same
	#
	
	if {[string tolower $curResType] != [lindex $stereopartner 0]} {
	    
	    continue
	}

	#
	# build a selection around the first stereopartner atom.  If it's equal to 
	# the given selection, we have a match.  Return a selection string based on 
	# the nonspecific atom name in the stereopartner entry
	#
	# Note that selections which cover the entire non-specific pair already won't 
	# be caught by this filter, but they don't need anything done to them anyway
	#

	set temp [AtomSel -args [format "(segid \"%s\" and resid %d and resn %s and name %s)" \
				     $curSeg $curNum $curResType [lindex $stereopartner 1]]]

	if {[$aSel isEqualTo [$temp cget -this]]} {

	    if {$origHasSegid} {
		set retVal [format "(segid \"%s\" and resid %d and resn %s and name %s)" \
				$curSeg $curNum $curResType [lindex $stereopartner 3]]
	    } else {
		set retVal [format "(resid %d and resn %s and name %s)" \
				$curNum $curResType [lindex $stereopartner 3]]
	    }

	    rename $temp ""
	    break
	}

	rename $temp ""
    }

    rename $aSel ""
    
    return $retVal
}



#
# converts, eg., (leu hd1*) to (leu hd2*)
#


proc stereoPartnerForSelection {aSelString} {

    set aSel [AtomSel -args $aSelString]

    #
    # complain if this isn't a single-residue selection
    #

    switch [numResiduesInSelection $aSel] {
	0 { rename $aSel ""
	    error [format "stereoPartnerForSelection called on an empty selection: %s" $aSelString] } 
	1 {}
	default { rename $aSel ""
	    error [format "stereoPartnerForSelection called on a selection that extends over >1 residue: %s. \n%s" \
		       $aSelString "Are the restraints' selections segid-specific?"] }
    }

    #
    # check if the original selection string includes a segid tag.  If not, 
    # don't create one in the output selection string
    #

    if {[string match -nocase "*segid*" $aSelString]} {
	set origHasSegid 1
    } else {
	set origHasSegid 0
    }

    #
    # grab the segid, resnum, and resname from one of the selected atoms
    #
    
    set curAtom [$aSel atomAtPos 0]
    Atom -this $curAtom

    set curSeg     [$curAtom segmentName]
    set curNum     [$curAtom residueNum]
    set curResType [$curAtom residueName]
    rename $curAtom ""

    set retVal $aSelString

    foreach stereopartner [stereopartnerList] {

	#
	# don't even bother checking further if the residue types
	# aren't the same
	#
	
	if {[string tolower $curResType] != [lindex $stereopartner 0]} {
	    
	    continue
	}

	#
	# build a selection around the first stereopartner atom.  If
	# it's equal to the given selection, we have a match.  Return
	# a selection string based on the stereopartner's atom name in
	# the stereopartner entry
	#
	# Note that selections which cover the entire non-specific
	# pair already won't be caught by this filter, but they don't
	# need anything done to them anyway
	#

	set temp [AtomSel \
		      args [format \
		"(segid \"%s\" and resid %d and resn %s and name %s)" \
		     $curSeg $curNum $curResType [lindex $stereopartner 1]]]

	if {[$aSel isEqualTo [$temp cget -this]]} {

	    if {$origHasSegid} {
		set retVal \
		    [format \
		"(segid \"%s\" and resid %d and resn %s and name %s)" \
			$curSeg $curNum $curResType [lindex $stereopartner 2]]
	    } else {
		set retVal [format "(resid %d and resn %s and name %s)" \
				$curNum $curResType [lindex $stereopartner 2]]
	    }
		
	    rename $temp ""
	    break
	}

	rename $temp ""
    }

    rename $aSel ""
    
    return $retVal
}


#
# Just a place to keep a selection for methyl groups.
# FIX:  need to add nucleic acids
#

proc methylSelectionString {} {
    
    return "((resn ala and name hb*) or
	     (resn val and name hg1*) or
	     (resn val and name hg2*) or
	     (resn ile and name hd*) or 
	     (resn ile and name hg2*) or
	     (resn leu and name hd1*) or
	     (resn leu and name hd2*) or 
	     (resn met and name he*) or 
	     (resn thr and name hg2*))"
}


#
# a quick way of converting (segid a and resid 4 and resn ala and name ca) to a_4_ala_ca
# for use in shiftAssignment names
#

proc compressedSelectionString {aSelString} {

    regsub -all "segidentifier " $aSelString "" temp
    regsub -all "segid " $aSelString "" temp
    regsub -all "residue " $temp "" temp
    regsub -all "resid " $temp "" temp
    regsub -all "resname " $temp "" temp
    regsub -all "resn " $temp "" temp
    regsub -all "name " $temp "" temp
    regsub -all "atom " $temp "" temp
    regsub -all "and " $temp "" temp

    regsub -all " " $temp "_" temp

    return $temp
} 
