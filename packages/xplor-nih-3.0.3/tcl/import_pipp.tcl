#
# TCL routines for processing PIPP .PCK and .DEG files
#

package provide pipp 1.0
package require marvin
package require aeneas

# need aeneas package for foldedOrAliasedPosition, used by process3dAssignedPippPeakTable, below


proc readPippShiftTable args {

    set fileName [requiredFlagVal $args -fileName]
    set remVar   [flagVal $args -remarksVariableName ""]

    global errorInfo

    if { [catch {set inUnit [open $fileName r]}] } {
    	error "Error opening input file $fileName"
    }

    set retVal [list]

    while {![eof $inUnit]} {

	if {[catch {eatNormalTextUpTo $inUnit "RES_ID"}]} {
	    set errorInfo ""
	    break
	}

	set curResNum [nextWord $inUnit]

	if { [catch {set curLine [nextLine $inUnit]}] } {
	    set errorInfo ""
	    break
	}

	#
	# we're inside of a residue definition.
	# process it
	#

	while {1} {

	    set curFirstWord [lindex $curLine 0]

	    switch -glob $curFirstWord {

		"RES_TYPE" {set curResType [lindex $curLine 1]}
		"PREV_TYPE" -
		"HETEROGENEITY" -
		"CUR_TYPE" -
		"SPIN_SYSTEM_ID" {}
		\#* {}
		*:-1 {}
		"END_RES_DEF" { break }
		default { 
		    set curAtomRawName [lindex $curLine 0]
		    set curShift [lindex $curLine 1]

		    # skip lines with atom names of (null)

		    if {$curAtomRawName == "(null)"} {
			break
		    }

		    # skip entries with res_id of * or res_type of UNK

		    if {($curResType == "UNK") || ($curResNum == "*")} {
			break
		    }
		    
		    # check for cases of A|B names
		    set curAtomNames [split $curAtomRawName "|"]

		    if {[llength $curAtomNames] == 1} {
			set curSelString [format "(resid %d and resn %s and name %s)" \
					      $curResNum $curResType [lindex $curAtomNames 0]]
			
			lappend retVal [list $curShift $curSelString]

		    } elseif {[llength $curAtomNames] == 2} {

			set curEntry [list $curShift]
			
			foreach curAtomName $curAtomNames {
			    set curSelString [format "(resid %d and resn %s and name %s)" \
						  $curResNum $curResType $curAtomName]
			    
			    lappend curEntry $curSelString 
			}

			lappend retVal $curEntry
		    }
		}
	    }
	    
	    if { [catch {set curLine [nextLine $inUnit]}] } {
		set errorInfo ""
		break
	    }
	}
    }

    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "%d shifts were read from PIPP formatted file %s" [llength $retVal] $fileName]
    }
    puts stdout [format "%d shifts were read from PIPP formatted file %s" \
	   [llength $retVal] $fileName]

    close $inUnit
    return $retVal
}



#
# reads the parameters in a PIPP .PCK table and 
# returns the center of the selected axis in PPM 
# and its width in PPM
#

proc readSpectralRangeFromPippPeakTable args {

    set fname [requiredFlagVal $args -fileName]
    set axis  [requiredFlagVal $args -axis]
    
    global errorInfo

    # open file
    
    if {[catch {set inUnit [open $fname r]}]} {
	error "Error opening input file $fname"
    }
    
    # eat text until we hit a line that starts with DATA <column name>
    set found 0
    while {! $found} {
	
	if {[catch {set l [nextLine $inUnit]}]} {
	    set errorInfo ""
	    break
	}
	
	if {([lindex $l 0] == "DATA") && ([lindex $l 1] == $axis)} {
	    set found 1
	}
    }
    
    close $inUnit

    if {! $found} {
	error "DATA line for axis named $axis not found in file $fname"
    }

    #
    # PIPP header format is 
    # DATA <axisName> <nPoints> <specWid Hz> <specFreq MHz> <freqFirstPt Hz> <ppmFlag> <widthEst ppm>
    #
    # where ppmFlag is "PPM" if the values in the peak lines are given in PPM.  
    # For now, just complain if they're not
    #
    # Dan Garrett says that PIPP can sometimes mess up the value of the linewidth estimate, 
    # so I won't use it
    #
    
    set nPoints          [lindex $l 2]
    set specWidHz        [lindex $l 3]
    set spectrometerFreq [lindex $l 4]
    set firstPtFreqHz    [lindex $l 5]
    set ppmFlag          [lindex $l 6]
    set lineWidEstPPM    [lindex $l 7]
    
    if {$ppmFlag != "PPM"} {
	error "PIPP header in file $fname indicates peak positions for axis $axis aren't in ppm."
    }
    
    #
    # the spectrum actually begins 1/2 of a point before the first point's frequency.
    # 
    
    set digitalRes        [expr $specWidHz / double($nPoints)]
    set spectrumStartFreq [expr $firstPtFreqHz + (0.5 * $digitalRes)]
    set spectrumEndFreq   [expr $spectrumStartFreq - $specWidHz]
    
    #
    # return spectral range in PPM
    #
	
    set specStartPPM [expr $spectrumStartFreq / double($spectrometerFreq)]
    set specEndPPM   [expr $spectrumEndFreq   / double($spectrometerFreq)]
    
    return [list $specStartPPM $specEndPPM]
}


proc readPippPeakTable args {

    set fname [requiredFlagVal $args -fileName]
    set colNames  [requiredFlagVal $args -cols]

    global errorInfo

    set retVal [list]

    # open file

    if {[catch {set inUnit [open $fname r]}]} {
	error "Error opening input file $fname"
    }

    #
    # read/parse the VARS line
    #

    # eat text until we hit a line that starts with VARS
    while {1} {

	set l [nextLine $inUnit]

	# skip all comments

	if {[string match \#* $l]} {
	    continue
	}
	
	if {[lindex $l 0] == "VARS"} {
	    break
	}
    }

    set columnPosns [list]

    foreach colName $colNames {
	lappend columnPosns [expr [lsearch -exact $l $colName] - 1]
    }

    # make sure none of these are missing

    foreach pos $columnPosns colName $colNames {
	if {$pos == -2} {
	    error "Can't find column $colName in VARS line: $l" 
	}
    }

    # read each line's data

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
	# skip FORMAT, NULLVALUE, NULLSTRING lines
	#

	if {([lindex $l 0] == "FORMAT") || ([lindex $l 0] == "NULLVALUE") || ([lindex $l 0] == "NULLSTRING")} {
	    continue
	}

	#
	# skip commented lines
	#

	if {[lindex $l 0] == "\#"} {
	    continue
	}

	set curVal [list]

	#
	# extract the columns' data
	#
	
	foreach pos $columnPosns {
	    lappend curVal [lindex $l $pos]
	}
	    
	#
	# make sure each column is defined
	#

	foreach elem $curVal {
	    if {$elem == ""} {
		error "Can't parse line $l from file $fname"
	    }
	}


	lappend retVal $curVal
    }

    close $inUnit
    return $retVal
}


# 
# for processing assigned PCK/DEG tables only.  
#

proc pippSel2xplor {selString} {

    # 
    # PIPP selections are in the form K25.CB,HB1;
    # M124.CG,HG1; L31.HD1#
    #

    #
    # if the string is ****, it means no selection
    # was given.  
    #

    if {$selString == "****"} {
	error "No selection given."
    }

    # First, lose the trailing semicolon if it's there

    if {[string match "*;" $selString]} {
	
	set l [string length $selString]
	set selString [string range $selString 0 [expr $l - 2]]

    }

    # Second, grab the residue number and type
    
    set pos [string first "." $selString]

    set resInfo [string range $selString 0 [expr $pos - 1] ]

    set resType [oneToThreeLetterResName [string range $resInfo 0 0]]
    set resNum  [string range $resInfo 1 end]

    # Third, grab the atom name.  If there's a comma in it,
    # then pull out the last section. (which is the proton's name)

    set name [string range $selString [expr $pos + 1] end]

    if {[string match "*,*" $name]} {

	set pos [string last "," $name]
	set name [string range $name [expr $pos + 1] end]
    }	

    # Fourth, break names of the form A|B into their 
    # constituent elements
    
    set nameList [split $name "|"]
    
    for {set nameCount 0} {$nameCount < [llength $nameList]} {incr nameCount} {

	set curName [lindex $nameList $nameCount]
	
	if {$nameCount == 0} {
	    set marvSelString [format "(resid %d and resname %s and (name %s" \
				   $resNum $resType $curName]
	} else {
	    set marvSelString [format "%s or name %s" $marvSelString $curName]
	}
    }
    
    set marvSelString [format "%s))" $marvSelString]

    return $marvSelString
}



#
# procs for reading assigned PIPP tables.  These rely on ShiftAssignments that are created 
# separately, from the shift table, rather than ones that are invented along the way, as is
# necessary for reading xplor restraints
#
# Note that non-stereo processing of ShiftAssigns should be done after this step, so that 
# the selections can be matched properly
#

proc process3dAssignedPippPeakTable args {

    global errorInfo

    set fnames                [requiredFlagVal $args -fileName]
    set pot                   [requiredFlagVal $args -pot]
    set fromProtShiftCol      [requiredFlagVal $args -fromProtonColumnName]
    set fromHeavyShiftCol     [requiredFlagVal $args -fromHeavyatomColumnName]
    set toProtShiftCol        [requiredFlagVal $args -toProtonColumnName]
    set fromSelCol            [requiredFlagVal $args -fromSelectionColumnName]
    set toSelCol              [requiredFlagVal $args -toSelectionColumnName]
    set namePrefix            [requiredFlagVal $args -namePrefix]
    set idCol                 [flagVal $args -peakIDcolumnName "PkID"]
    set intensityCol          [flagVal $args -intensityColumnName "Intensity"]
    set remVar                [flagVal $args -remarksVariableName ""]
    set fromProtTol           [flagVal $args -fromProtonTolerance 0.02]
    set fromHeavyTol          [flagVal $args -fromHeavyatomTolerance 0.1]
    set toProtTol             [flagVal $args -toProtonTolerance 0.02]
    set fromProtonRange       [requiredFlagVal $args -fromProtonSpectralRangePPM]
    set fromHeavyRange        [requiredFlagVal $args -fromHeavyatomSpectralRangePPM]
    set toProtonRange         [requiredFlagVal $args -toProtonSpectralRangePPM]
    set foldedAlongFromProton [flagExists $args -shiftsFoldedAlongFromProtonDimension]
    set foldedAlongFromHeavy  [flagExists $args -shiftsFoldedAlongFromHeavyatomDimension]
    set foldedAlongToProton   [flagExists $args -shiftsFoldedAlongToProtonDimension]
    set fromProtonSignChanges [flagExists $args [list -signChangesUponFoldingFromProtonDimension \
						     -signChangesUponAliasingFromProtonDimension]]

    set fromHeavySignChanges  [flagExists $args [list -signChangesUponFoldingFromHeavyatomDimension \
						     -signChangesUponAliasingFromHeavyatomDimension]]

    set toProtonSignChanges   [flagExists $args [list -signChangesUponFoldingToProtonDimension \
						     -signChangesUponAliasingToProtonDimension]]

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
    # extract spectral top/bottom from range flags in order to record unfolded peak positions
    #

    set fromProtonTop [max $fromProtonRange]
    set fromProtonBot [min $fromProtonRange]
    set fromHeavyTop  [max $fromHeavyRange]
    set fromHeavyBot  [min $fromHeavyRange]
    set toProtonTop   [max $toProtonRange]
    set toProtonBot   [min $toProtonRange]

    #
    # Need to record unfolded positions of each shiftAssign for use in 
    # explicit inverse exceptions.  Do it here because I have the spectral 
    # ranges and so forth
    #
    
    recordUnfoldedPositions \
	-pot $pot \
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
    # read the PCK and/or DEG table
    #

    set peaks [list]

    foreach fname $fnames {
	set peaks [concat $peaks [readPippPeakTable -fileName $fname -cols [list $idCol $intensityCol $fromProtShiftCol $fromHeavyShiftCol $toProtShiftCol $fromSelCol $toSelCol]]]
    }

    set nPeaksAdded 0
    set nPeakAssignsAdded 0
    set peaksWithNoPAs [list]
    set peaksWithNonMatchingPAs [list]
    set peaksWithMissingShiftAssignments [list]
    set peaksWithMultipleShiftAssignments [list]


    marvinPyth command {import pasd}
    set nMono [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.nMono} \
			   {} {out}] \
		      0] 1]
    set aveExp [lindex \
		    [lindex \
			 [marvinPyth command {out=pasd.aveExp} \
			      {} {out}] \
			 0] 1]

    set peakCount 0

    foreach peak $peaks {

	set curID       [lindex $peak 0]
	set curPeakName [format "%s%s" $namePrefix $curID]

	updateUser [format "Processing assigned PIPP peak %s (%d of %d)\r" $curPeakName [incr peakCount] [llength $peaks]]

	#
	# See if this is a new peak.  If so, create it.
	#

	if {! [$pot hasPeakNamed $curPeakName]} {
	    
	    set newPeak [Peak -args $curPeakName]
	    $newPeak setIntensity          [lindex $peak 1]
	    $newPeak setFromProtonShift    [lindex $peak 2]
	    $newPeak setFromHeavyatomShift [lindex $peak 3]
	    $newPeak setToProtonShift      [lindex $peak 4]
	    $newPeak appendToNote          [format "from assigned PIPP file %s, peak %s" $fname $curID]
	    $pot addPeak [$newPeak cget -this]
	    $newPeak -disown
	    rename $newPeak ""
	    incr nPeaksAdded
	}

	#
	# look up the peak with the current name.  Note that
	# this assumes that lines referring to peaks with existing 
	# names are additional peakAssigns for those peaks (as happens
	# in DEG tables).  
	#
	
	set curPeak [$pot peakNamed $curPeakName]
	Peak -this $curPeak
	
	#
	# complain if this peak's location is too far off from the location in the PIPP file line I just read,
	# but still create the PeakAssignment
	#
	
	set deltaFromProt  [expr abs([$curPeak fromProtonShift]    - [lindex $peak 2])]
	set deltaFromHeavy [expr abs([$curPeak fromHeavyatomShift] - [lindex $peak 3])]
	set deltaToProt    [expr abs([$curPeak toProtonShift]      - [lindex $peak 4])]
	
	set isBad 0

	if {$deltaFromProt > $fromProtTol} {
	    $curPeak appendToNote [format "WARNING:  large difference in from proton shift in multiple assignments for this peak: %f %f" \
				       [$curPeak fromProtonShift] [lindex $peak 2]]
	    set isBad 1
	}

	if {$deltaFromHeavy > $fromHeavyTol} {
	    $curPeak appendToNote [format "WARNING:  large difference in from heavyatom shift in multiple assignments for this peak: %f %f" \
				       [$curPeak fromHeavyatomShift] [lindex $peak 3]]
	    set isBad 1
	}

	if {$deltaToProt > $toProtTol} {
	    $curPeak appendToNote [format "WARNING:  large difference in to proton shift in multiple assignments for this peak: %f %f" \
				       [$curPeak toProtonShift] [lindex $peak 4]]
	    set isBad 1
	}

	if {$isBad} {
	    lappend peaksWithNonMatchingPAs [$curPeak name]
	}

	#
	# try converting the selections given in the PIPP table to xplor format.
	# Note that an unassigned peak will throw an error, in which case we want 
	# to give up on any peakAssignments
	#

	if {[catch {set fromProtonSelString [pippSel2xplor [lindex $peak 5]]}]} {
	    $curPeak appendToNote [format "WARNING:  No peak assignments, since from assignment given in PCK table was %s" [lindex $peak 5]]
	    lappend peaksWithNoPAs [$curPeak name]
	    rename $curPeak ""
	    set errorInfo ""
	    continue
	}
	    
	if {[catch {set toProtonSelString [pippSel2xplor [lindex $peak 6]]}]} {
	    $curPeak appendToNote [format "WARNING:  No peak assignments, since to assignment given in PCK table was %s" [lindex $peak 6]]
	    lappend peaksWithNoPAs [$curPeak name]
	    rename $curPeak ""
	    set errorInfo ""
	    continue
	}

	set fromProtonSel [AtomSel -args $fromProtonSelString]
	set toProtonSel   [AtomSel -args $toProtonSelString]

	#
	# search for matching shiftAssignments
	#

	set matchingFromSANames [shiftAssignmentsInvolvingSelectedAtoms -shiftAssignments [$pot fromShiftAssignments] -sel $fromProtonSel]
	set matchingToSANames   [shiftAssignmentsInvolvingSelectedAtoms -shiftAssignments [$pot toShiftAssignments]   -sel $toProtonSel]

	#
	# complain if I have no matching SAs
	#

	if {[llength $matchingFromSANames] == 0} {
	    $curPeak appendToNote [format "ERROR:  No ShiftAssignment matches given PIPP from proton selection %s (converted to %s)" [lindex $peak 5] [$fromProtonSel string]]
	    lappend peaksWithMissingShiftAssignments [$curPeak name]
	}

	if {[llength $matchingToSANames] == 0} {
	    $curPeak appendToNote [format "ERROR:  No ShiftAssignment matches given PIPP to proton selection %s (converted to %s)" [lindex $peak 6] [$toProtonSel string]]
	    lappend peaksWithMissingShiftAssignments [$curPeak name]
	}


	#
	# I can have multiple matching SAs if my shift table has non-stereo entries (or other degeneracy) in it.  
	# If so, filter out any matching shiftAssignments whose shifts don't fit the peak within tolerance.
	#

	if {[llength $matchingFromSANames] > 1} {
	
	    set shiftMatches [list]

	    foreach matchSAname $matchingFromSANames {
		set matchSA [$pot shiftAssignmentNamed $matchSAname]
		ShiftAssignment -this $matchSA
		
		set unfoldedFromProtonPPs [unfoldPeakPosnToMatch [$curPeak fromProtonShift]    [$matchSA protonShift]    $foldedAlongFromProton $fromProtonTop $fromProtonBot $fromProtonSignChanges] 
		set unfoldedFromHeavyPPs  [unfoldPeakPosnToMatch [$curPeak fromHeavyatomShift] [$matchSA heavyatomShift] $foldedAlongFromHeavy  $fromHeavyTop  $fromHeavyBot  $fromHeavySignChanges] 

		foreach protonPP $unfoldedFromProtonPPs {
		    if {[unfoldedPeakPosnMatches $protonPP [$matchSA protonShift] $fromProtTol]} {
			foreach heavyPP $unfoldedFromHeavyPPs {
			    if {[unfoldedPeakPosnMatches $heavyPP [$matchSA heavyatomShift] $fromHeavyTol]} {
				lappend shiftMatches $matchSAname
			    }
			}
		    }
		}
		
		rename $matchSA ""
	    }

	    set shiftMatches [lsort -unique $shiftMatches]
	    
	    #
	    # if any shiftAssignments survived the filter, 
	    # only work with those
	    #

	    if {[llength $shiftMatches] > 0} {
		set matchingFromSANames $shiftMatches
	    }
	}

	if {[llength $matchingToSANames] > 1} {

	    set shiftMatches [list]
	    foreach matchSAname $matchingToSANames {
		set matchSA [$pot shiftAssignmentNamed $matchSAname]
		ShiftAssignment -this $matchSA
		
		set unfoldedToProtonPPs [unfoldPeakPosnToMatch [$curPeak toProtonShift] [$matchSA protonShift] $foldedAlongToProton $toProtonTop $toProtonBot $toProtonSignChanges] 

		foreach protonPP $unfoldedToProtonPPs {
		    if {[unfoldedPeakPosnMatches $protonPP [$matchSA protonShift] $toProtTol]} {
			lappend shiftMatches $matchSAname
		    }
		}

		rename $matchSA ""
	    }

	    set shiftMatches [lsort -unique $shiftMatches]
	    
	    #
	    # if any shiftAssignments survived the filter, 
	    # only work with those
	    #
	    
	    if {[llength $shiftMatches] > 0} {
		set matchingToSANames $shiftMatches
	    }
	}

	#
	# record any complaints about multiple shiftAssignments
	#

	if {[llength $matchingFromSANames] > 1} {
	    $curPeak appendToNote [format "NOTE:  %d ShiftAssignments match given PIPP from proton selection %s" [llength $matchingFromSANames] [lindex $peak 5]]
	    lappend peaksWithMultipleShiftAssignments [$curPeak name]
	}

	if {[llength $matchingToSANames] > 1} {
	    $curPeak appendToNote [format "NOTE:  %d ShiftAssignments match given PIPP to proton selection %s" [llength $matchingToSANames] [lindex $peak 6]]
	    lappend peaksWithMultipleShiftAssignments [$curPeak name]
	}

	#
	# create the PeakAssignment(s)
	#

	foreach fromSAname $matchingFromSANames {
	    foreach toSAname $matchingToSANames {

		set newPAname [format "%s_%d" [$curPeak name] [$curPeak numPeakAssignments]]
		set newPA [PeakAssignment -args $newPAname]
		$newPA setNMono $nMono
		$newPA setAveExp $aveExp

		$curPeak addPeakAssignment [$newPA cget -this]

		set fromSA [$pot shiftAssignmentNamed $fromSAname]
		set toSA   [$pot shiftAssignmentNamed $toSAname]
	
		ShiftAssignment -this $fromSA
		ShiftAssignment -this $toSA

		$newPA setFromAssignment [$fromSA cget -this]
		$newPA setToAssignment   [$toSA cget -this]

		$newPA setPreviousLikelihood 1.0

		#
		# Unfold the NOE peak's position and record it, so I can do stripe likelihood
		#
		# To do this correctly, I need to match the peak's position roughly and then match the 
		# peak's sign.  I'll use specWid/2 as the tolerance
		#

		set unfoldedFromProtonPPs [unfoldPeakPosnToMatch [$curPeak fromProtonShift]    [$fromSA protonShift]    $foldedAlongFromProton $fromProtonTop $fromProtonBot $fromProtonSignChanges] 
		set unfoldedFromHeavyPPs  [unfoldPeakPosnToMatch [$curPeak fromHeavyatomShift] [$fromSA heavyatomShift] $foldedAlongFromHeavy  $fromHeavyTop  $fromHeavyBot  $fromHeavySignChanges] 
		set unfoldedToProtonPPs   [unfoldPeakPosnToMatch [$curPeak toProtonShift]      [$toSA   protonShift]    $foldedAlongToProton   $toProtonTop   $toProtonBot   $toProtonSignChanges] 

		set fromProtSpecWid2  [expr ($fromProtonTop - $fromProtonBot) / double(2)]
		set fromHeavySpecWid2 [expr ($fromHeavyTop  - $fromHeavyBot)  / double(2)]
		set toProtSpecWid2    [expr ($toProtonTop   - $toProtonBot)   / double(2)]

		foreach fromProtonPP $unfoldedFromProtonPPs {
		    if {[unfoldedPeakPosnMatches $fromProtonPP [$fromSA protonShift] $fromProtSpecWid2]} {
			foreach fromHeavyPP $unfoldedFromHeavyPPs {
			    if {[unfoldedPeakPosnMatches $fromHeavyPP [$fromSA heavyatomShift] $fromHeavySpecWid2]} {
				foreach toProtonPP $unfoldedToProtonPPs {
				    if {[unfoldedPeakPosnMatches $toProtonPP [$toSA protonShift] $toProtSpecWid2]} {
					
					set expectedSign [expr [lindex $fromProtonPP 1] * [lindex $fromHeavyPP 1] * [lindex $toProtonPP 1] * $unfoldedSign]
					if {(($expectedSign == [sign [$curPeak intensity]]) || $ignoreSign)} {

					    $newPA setUnfoldedFromProtonPeakPosition    [lindex $fromProtonPP 0]
					    $newPA setUnfoldedFromHeavyatomPeakPosition [lindex $fromHeavyPP 0]
					    $newPA setUnfoldedToProtonPeakPosition      [lindex $toProtonPP 0]
					}
				    }
				}
			    }
			}
		    }
		}
	    
		$newPA -disown
		rename $newPA ""
		rename $fromSA ""
		rename $toSA ""

		incr nPeakAssignsAdded
	    }
	}
	
	rename $fromProtonSel ""
	rename $toProtonSel ""
	rename $curPeak ""
    }

    set peaksWithNoPAs                    [lsort -dictionary -unique $peaksWithNoPAs]
    set peaksWithNonMatchingPAs           [lsort -dictionary -unique $peaksWithNonMatchingPAs]
    set peaksWithMissingShiftAssignments  [lsort -dictionary -unique $peaksWithMissingShiftAssignments]
    set peaksWithMultipleShiftAssignments [lsort -dictionary -unique $peaksWithMultipleShiftAssignments]


    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "Created %d peaks and %d peakAssignments from assigned PIPP file %s" $nPeaksAdded $nPeakAssignsAdded $fname]

	if {[llength $peaksWithNoPAs] > 0} {
	    set rpt [format "%d peaks had no peakAssignments (usually from uninterpretable PIPP selections).  They are:" [llength $peaksWithNoPAs]]
	    foreach elem $peaksWithNoPAs {
		appendLineToString rpt [format "   %s" $elem]
	    }
	    lappend tempRem $rpt
	}


	if {[llength $peaksWithNonMatchingPAs] > 0} {
	    set rpt [format "%d peaks had multiple peakAssignments whose peak locations disagreed.  They are:" [llength $peaksWithNonMatchingPAs]]
	    foreach elem $peaksWithNonMatchingPAs {
		appendLineToString rpt [format "   %s" $elem]
	    }
	    lappend tempRem $rpt
	}

	if {[llength $peaksWithMissingShiftAssignments] > 0} {
	    set rpt [format "%d peaks had assignments whose selections did not match any ShiftAssignment.  They are:" [llength $peaksWithMissingShiftAssignments]]
	    foreach elem $peaksWithMissingShiftAssignments {
		appendLineToString rpt [format "   %s" $elem]
	    }
	    lappend tempRem $rpt
	}

	if {[llength $peaksWithMultipleShiftAssignments] > 0} {
	    set rpt [format "%d peaks had assignments whose selections matched > 1 ShiftAssignment.  They are:" [llength $peaksWithMultipleShiftAssignments]]
	    foreach elem $peaksWithMultipleShiftAssignments {
		appendLineToString rpt [format "   %s" $elem]
	    }
	    lappend tempRem $rpt
	}
    }
    
    return ""
}
    







proc process4dAssignedPippPeakTable args {

    global errorInfo

    set fnames                [requiredFlagVal $args -fileName]
    set pot                   [requiredFlagVal $args -pot]
    set fromProtShiftCol      [requiredFlagVal $args -fromProtonColumnName]
    set fromHeavyShiftCol     [requiredFlagVal $args -fromHeavyatomColumnName]
    set toProtShiftCol        [requiredFlagVal $args -toProtonColumnName]
    set toHeavyShiftCol       [requiredFlagVal $args -toHeavyatomColumnName]
    set fromSelCol            [requiredFlagVal $args -fromSelectionColumnName]
    set toSelCol              [requiredFlagVal $args -toSelectionColumnName]
    set namePrefix            [requiredFlagVal $args -namePrefix]
    set idCol                 [flagVal $args -peakIDcolumnName "PkID"]
    set intensityCol          [flagVal $args -intensityColumnName "Intensity"]
    set remVar                [flagVal $args -remarksVariableName ""]
    set fromProtTol           [flagVal $args -fromProtonTolerance 0.02]
    set fromHeavyTol          [flagVal $args -fromHeavyatomTolerance 0.1]
    set toProtTol             [flagVal $args -toProtonTolerance 0.02]
    set toHeavyTol            [flagVal $args -toHeavyatomTolerance 0.1]
    set fromProtonRange       [requiredFlagVal $args -fromProtonSpectralRangePPM]
    set fromHeavyRange        [requiredFlagVal $args -fromHeavyatomSpectralRangePPM]
    set toProtonRange         [requiredFlagVal $args -toProtonSpectralRangePPM]
    set toHeavyRange          [requiredFlagVal $args -toHeavyatomSpectralRangePPM]
    set foldedAlongFromProton [flagExists $args -shiftsFoldedAlongFromProtonDimension]
    set foldedAlongFromHeavy  [flagExists $args -shiftsFoldedAlongFromHeavyatomDimension]
    set foldedAlongToProton   [flagExists $args -shiftsFoldedAlongToProtonDimension]
    set foldedAlongToHeavy    [flagExists $args -shiftsFoldedAlongToHeavyatomDimension]
    set fromProtonSignChanges [flagExists $args [list -signChangesUponFoldingFromProtonDimension \
						     -signChangesUponAliasingFromProtonDimension]]

    set fromHeavySignChanges  [flagExists $args [list -signChangesUponFoldingFromHeavyatomDimension \
						     -signChangesUponAliasingFromHeavyatomDimension]]

    set toProtonSignChanges   [flagExists $args [list -signChangesUponFoldingToProtonDimension \
						     -signChangesUponAliasingToProtonDimension]]

    set toHeavySignChanges  [flagExists $args [list -signChangesUponFoldingToHeavyatomDimension \
						   -signChangesUponAliasingToHeavyatomDimension]]

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
    # extract spectral top/bottom from range flags in order to record unfolded peak positions
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
    # Need to record unfolded positions of each shiftAssign for use in 
    # explicit inverse exceptions.  Do it here because I have the spectral 
    # ranges and so forth
    #
    
    recordUnfoldedPositions \
	-pot $pot \
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
	-shiftsFoldedAlongToHeavyatomDimension   $foldedAlongToHeavy 


     
    #
    # read the PCK and/or DEG table
    #

    set peaks [list]

    foreach fname $fnames {
	set peaks [concat $peaks [readPippPeakTable -fileName $fname -cols [list $idCol $intensityCol $fromProtShiftCol $fromHeavyShiftCol $toProtShiftCol $toHeavyShiftCol $fromSelCol $toSelCol]]]
    }

    set nPeaksAdded 0
    set nPeakAssignsAdded 0
    set peaksWithNoPAs [list]
    set peaksWithNonMatchingPAs [list]
    set peaksWithMissingShiftAssignments [list]
    set peaksWithMultipleShiftAssignments [list]

    marvinPyth command {import pasd}
    set nMono [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.nMono} \
			   {} {out}] \
		      0] 1]
    set aveExp [lindex \
		    [lindex \
			 [marvinPyth command {out=pasd.aveExp} \
			      {} {out}] \
			 0] 1]

    set peakCount 0

    foreach peak $peaks {

	set curID       [lindex $peak 0]
	set curPeakName [format "%s%s" $namePrefix $curID]

	updateUser [format "Processing assigned PIPP peak %s (%d of %d)\r" $curPeakName [incr peakCount] [llength $peaks]]

	#
	# See if this is a new peak.  If so, create it.
	#

	if {! [$pot hasPeakNamed $curPeakName]} {
	    
	    set newPeak [Peak -args $curPeakName]
	    $newPeak setIntensity          [lindex $peak 1]
	    $newPeak setFromProtonShift    [lindex $peak 2]
	    $newPeak setFromHeavyatomShift [lindex $peak 3]
	    $newPeak setToProtonShift      [lindex $peak 4]
	    $newPeak setToHeavyatomShift   [lindex $peak 5]
	    $newPeak appendToNote          [format "from assigned PIPP file %s, peak %s" $fname $curID]
	    $pot addPeak [$newPeak cget -this]
	    $newPeak -disown
	    rename $newPeak ""
	    incr nPeaksAdded
	}

	#
	# look up the peak with the current name.  Note that
	# this assumes that lines referring to peaks with existing 
	# names are additional peakAssigns for those peaks (as happens
	# in DEG tables).  
	#
	
	set curPeak [$pot peakNamed $curPeakName]
	Peak -this $curPeak
	
	#
	# complain if this peak's location is too far off from the location in the PIPP file line I just read,
	# but still create the PeakAssignment
	#
	
	set deltaFromProt  [expr abs([$curPeak fromProtonShift]    - [lindex $peak 2])]
	set deltaFromHeavy [expr abs([$curPeak fromHeavyatomShift] - [lindex $peak 3])]
	set deltaToProt    [expr abs([$curPeak toProtonShift]      - [lindex $peak 4])]
	set deltaToHeavy   [expr abs([$curPeak toHeavyatomShift]   - [lindex $peak 5])]
	
	set isBad 0

	if {$deltaFromProt > $fromProtTol} {
	    $curPeak appendToNote [format "WARNING:  large difference in from proton shift in multiple assignments for this peak: %f %f" \
				       [$curPeak fromProtonShift] [lindex $peak 2]]
	    set isBad 1
	}

	if {$deltaFromHeavy > $fromHeavyTol} {
	    $curPeak appendToNote [format "WARNING:  large difference in from heavyatom shift in multiple assignments for this peak: %f %f" \
				       [$curPeak fromHeavyatomShift] [lindex $peak 3]]
	    set isBad 1
	}

	if {$deltaToProt > $toProtTol} {
	    $curPeak appendToNote [format "WARNING:  large difference in to proton shift in multiple assignments for this peak: %f %f" \
				       [$curPeak toProtonShift] [lindex $peak 4]]
	    set isBad 1
	}

	if {$deltaToHeavy > $toHeavyTol} {
	    $curPeak appendToNote [format "WARNING:  large difference in to heavyatom shift in multiple assignments for this peak: %f %f" \
				       [$curPeak toHeavyatomShift] [lindex $peak 5]]
	    set isBad 1
	}

	if {$isBad} {
	    lappend peaksWithNonMatchingPAs [$curPeak name]
	}

	#
	# try converting the selections given in the PIPP table to xplor format.
	# Note that an unassigned peak will throw an error, in which case we want 
	# to give up on any peakAssignments
	#

	if {[catch {set fromProtonSelString [pippSel2xplor [lindex $peak 6]]}]} {
	    $curPeak appendToNote [format "WARNING:  No peak assignments, since from assignment given in PCK table was %s" [lindex $peak 6]]
	    lappend peaksWithNoPAs [$curPeak name]
	    rename $curPeak ""
	    set errorInfo ""
	    continue
	}
	    
	if {[catch {set toProtonSelString [pippSel2xplor [lindex $peak 7]]}]} {
	    $curPeak appendToNote [format "WARNING:  No peak assignments, since to assignment given in PCK table was %s" [lindex $peak 7]]
	    lappend peaksWithNoPAs [$curPeak name]
	    rename $curPeak ""
	    set errorInfo ""
	    continue
	}

	set fromProtonSel [AtomSel -args $fromProtonSelString]
	set toProtonSel   [AtomSel -args $toProtonSelString]

	#
	# search for matching shiftAssignments
	#

	set matchingFromSANames [shiftAssignmentsInvolvingSelectedAtoms -shiftAssignments [$pot fromShiftAssignments] -sel $fromProtonSel]
	set matchingToSANames   [shiftAssignmentsInvolvingSelectedAtoms -shiftAssignments [$pot toShiftAssignments]   -sel $toProtonSel]

	#
	# complain if I have no matching SAs
	#

	if {[llength $matchingFromSANames] == 0} {
	    $curPeak appendToNote [format "ERROR:  No ShiftAssignment matches given PIPP from proton selection %s (converted to %s)" [lindex $peak 6] [$fromProtonSel string]]
	    lappend peaksWithMissingShiftAssignments [$curPeak name]
	}

	if {[llength $matchingToSANames] == 0} {
	    $curPeak appendToNote [format "ERROR:  No ShiftAssignment matches given PIPP to proton selection %s (converted to %s)" [lindex $peak 7] [$toProtonSel string]]
	    lappend peaksWithMissingShiftAssignments [$curPeak name]
	}


	#
	# I can have multiple matching SAs if my shift table has non-stereo entries (or other degeneracy) in it.  
	# If so, filter out any matching shiftAssignments whose shifts don't fit the peak within tolerance.
	#

	if {[llength $matchingFromSANames] > 1} {
	
	    set shiftMatches [list]

	    foreach matchSAname $matchingFromSANames {
		set matchSA [$pot shiftAssignmentNamed $matchSAname]
		ShiftAssignment -this $matchSA
		
		set unfoldedFromProtonPPs [unfoldPeakPosnToMatch [$curPeak fromProtonShift]    [$matchSA protonShift]    $foldedAlongFromProton $fromProtonTop $fromProtonBot $fromProtonSignChanges] 
		set unfoldedFromHeavyPPs  [unfoldPeakPosnToMatch [$curPeak fromHeavyatomShift] [$matchSA heavyatomShift] $foldedAlongFromHeavy  $fromHeavyTop  $fromHeavyBot  $fromHeavySignChanges] 

		foreach protonPP $unfoldedFromProtonPPs {
		    if {[unfoldedPeakPosnMatches $protonPP [$matchSA protonShift] $fromProtTol]} {
			foreach heavyPP $unfoldedFromHeavyPPs {
			    if {[unfoldedPeakPosnMatches $heavyPP [$matchSA heavyatomShift] $fromHeavyTol]} {
				lappend shiftMatches $matchSAname
			    }
			}
		    }
		}
		
		rename $matchSA ""
	    }

	    set shiftMatches [lsort -unique $shiftMatches]
	    
	    #
	    # if any shiftAssignments survived the filter, 
	    # only work with those
	    #

	    if {[llength $shiftMatches] > 0} {
		set matchingFromSANames $shiftMatches
	    }
	}

	if {[llength $matchingToSANames] > 1} {

	    set shiftMatches [list]
	    foreach matchSAname $matchingToSANames {
		set matchSA [$pot shiftAssignmentNamed $matchSAname]
		ShiftAssignment -this $matchSA
		
		set unfoldedToProtonPPs [unfoldPeakPosnToMatch [$curPeak toProtonShift]    [$matchSA protonShift]    $foldedAlongToProton $toProtonTop $toProtonBot $toProtonSignChanges] 
		set unfoldedToHeavyPPs  [unfoldPeakPosnToMatch [$curPeak toHeavyatomShift] [$matchSA heavyatomShift] $foldedAlongToHeavy  $toHeavyTop  $toHeavyBot  $toHeavySignChanges] 

		foreach protonPP $unfoldedToProtonPPs {
		    if {[unfoldedPeakPosnMatches $protonPP [$matchSA protonShift] $toProtTol]} {
			foreach heavyPP $unfoldedFromHeavyPPs {
			    if {[unfoldedPeakPosnMatches $heavyPP [$matchSA heavyatomShift] $toHeavyTol]} {
				lappend shiftMatches $matchSAname
			    }
			}
		    }
		}

		rename $matchSA ""
	    }

	    set shiftMatches [lsort -unique $shiftMatches]
	    
	    #
	    # if any shiftAssignments survived the filter, 
	    # only work with those
	    #
	    
	    if {[llength $shiftMatches] > 0} {
		set matchingToSANames $shiftMatches
	    }
	}

	#
	# record any complaints about multiple shiftAssignments
	#

	if {[llength $matchingFromSANames] > 1} {
	    $curPeak appendToNote [format "NOTE:  %d ShiftAssignments match given PIPP from proton selection %s" [llength $matchingFromSANames] [lindex $peak 6]]
	    lappend peaksWithMultipleShiftAssignments [$curPeak name]
	}

	if {[llength $matchingToSANames] > 1} {
	    $curPeak appendToNote [format "NOTE:  %d ShiftAssignments match given PIPP to proton selection %s" [llength $matchingToSANames] [lindex $peak 7]]
	    lappend peaksWithMultipleShiftAssignments [$curPeak name]
	}

	#
	# create the PeakAssignment(s)
	#

	foreach fromSAname $matchingFromSANames {
	    foreach toSAname $matchingToSANames {

		set newPAname [format "%s_%d" [$curPeak name] [$curPeak numPeakAssignments]]
		set newPA [PeakAssignment -args $newPAname]
		$newPA setNMono $nMono
		$newPA setAveExp $aveExp

		$curPeak addPeakAssignment [$newPA cget -this]

		set fromSA [$pot shiftAssignmentNamed $fromSAname]
		set toSA   [$pot shiftAssignmentNamed $toSAname]
	
		ShiftAssignment -this $fromSA
		ShiftAssignment -this $toSA

		$newPA setFromAssignment [$fromSA cget -this]
		$newPA setToAssignment   [$toSA cget -this]

		$newPA setPreviousLikelihood 1.0

		#
		# Unfold the NOE peak's position and record it, so I can do stripe likelihood
		#
		# To do this correctly, I need to match the peak's position roughly and then match the 
		# peak's sign.  I'll use specWid/2 as the tolerance
		#

		set unfoldedFromProtonPPs [unfoldPeakPosnToMatch [$curPeak fromProtonShift]    [$fromSA protonShift]    $foldedAlongFromProton $fromProtonTop $fromProtonBot $fromProtonSignChanges] 
		set unfoldedFromHeavyPPs  [unfoldPeakPosnToMatch [$curPeak fromHeavyatomShift] [$fromSA heavyatomShift] $foldedAlongFromHeavy  $fromHeavyTop  $fromHeavyBot  $fromHeavySignChanges] 
		set unfoldedToProtonPPs   [unfoldPeakPosnToMatch [$curPeak toProtonShift]      [$toSA   protonShift]    $foldedAlongToProton   $toProtonTop   $toProtonBot   $toProtonSignChanges] 
		set unfoldedToHeavyPPs    [unfoldPeakPosnToMatch [$curPeak toHeavyatomShift]   [$toSA heavyatomShift]   $foldedAlongToHeavy    $toHeavyTop    $toHeavyBot    $toHeavySignChanges] 

		set fromProtSpecWid2  [expr ($fromProtonTop - $fromProtonBot) / double(2)]
		set fromHeavySpecWid2 [expr ($fromHeavyTop  - $fromHeavyBot)  / double(2)]
		set toProtSpecWid2    [expr ($toProtonTop   - $toProtonBot)   / double(2)]
		set toHeavySpecWid2   [expr ($toHeavyTop    - $toHeavyBot)    / double(2)]

		foreach fromProtonPP $unfoldedFromProtonPPs {
		    if {[unfoldedPeakPosnMatches $fromProtonPP [$fromSA protonShift] $fromProtSpecWid2]} {
			foreach fromHeavyPP $unfoldedFromHeavyPPs {
			    if {[unfoldedPeakPosnMatches $fromHeavyPP [$fromSA heavyatomShift] $fromHeavySpecWid2]} {
				foreach toProtonPP $unfoldedToProtonPPs {
				    if {[unfoldedPeakPosnMatches $toProtonPP [$toSA protonShift] $toProtSpecWid2]} {
					foreach toHeavyPP $unfoldedToHeavyPPs {
					    if {[unfoldedPeakPosnMatches $toHeavyPP [$toSA heavyatomShift] $toHeavySpecWid2]} {

						set expectedSign [expr [lindex $fromProtonPP 1] * [lindex $fromHeavyPP 1] * [lindex $toProtonPP 1] * [lindex $toHeavyPP 1] * $unfoldedSign]
						if {(($expectedSign == [sign [$curPeak intensity]]) || $ignoreSign)} {

						    $newPA setUnfoldedFromProtonPeakPosition    [lindex $fromProtonPP 0]
						    $newPA setUnfoldedFromHeavyatomPeakPosition [lindex $fromHeavyPP 0]
						    $newPA setUnfoldedToProtonPeakPosition      [lindex $toProtonPP 0]
						    $newPA setUnfoldedToHeavyatomPeakPosition   [lindex $toHeavyPP 0]
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    
		$newPA -disown
		rename $newPA ""
		rename $fromSA ""
		rename $toSA ""

		incr nPeakAssignsAdded
	    }
	}
	
	rename $fromProtonSel ""
	rename $toProtonSel ""
	rename $curPeak ""
    }

    set peaksWithNoPAs                    [lsort -dictionary -unique $peaksWithNoPAs]
    set peaksWithNonMatchingPAs           [lsort -dictionary -unique $peaksWithNonMatchingPAs]
    set peaksWithMissingShiftAssignments  [lsort -dictionary -unique $peaksWithMissingShiftAssignments]
    set peaksWithMultipleShiftAssignments [lsort -dictionary -unique $peaksWithMultipleShiftAssignments]


    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "Created %d peaks and %d peakAssignments from assigned PIPP file %s" $nPeaksAdded $nPeakAssignsAdded $fname]

	if {[llength $peaksWithNoPAs] > 0} {
	    set rpt [format "%d peaks had no peakAssignments (usually from uninterpretable PIPP selections).  They are:" [llength $peaksWithNoPAs]]
	    foreach elem $peaksWithNoPAs {
		appendLineToString rpt [format "   %s" $elem]
	    }
	    lappend tempRem $rpt
	}


	if {[llength $peaksWithNonMatchingPAs] > 0} {
	    set rpt [format "%d peaks had multiple peakAssignments whose peak locations disagreed.  They are:" [llength $peaksWithNonMatchingPAs]]
	    foreach elem $peaksWithNonMatchingPAs {
		appendLineToString rpt [format "   %s" $elem]
	    }
	    lappend tempRem $rpt
	}

	if {[llength $peaksWithMissingShiftAssignments] > 0} {
	    set rpt [format "%d peaks had assignments whose selections did not match any ShiftAssignment.  They are:" [llength $peaksWithMissingShiftAssignments]]
	    foreach elem $peaksWithMissingShiftAssignments {
		appendLineToString rpt [format "   %s" $elem]
	    }
	    lappend tempRem $rpt
	}

	if {[llength $peaksWithMultipleShiftAssignments] > 0} {
	    set rpt [format "%d peaks had assignments whose selections matched > 1 ShiftAssignment.  They are:" [llength $peaksWithMultipleShiftAssignments]]
	    foreach elem $peaksWithMultipleShiftAssignments {
		appendLineToString rpt [format "   %s" $elem]
	    }
	    lappend tempRem $rpt
	}
    }
    
    return ""
}











#
# top-level import routines for un-assigned PIPP peak tables
#

proc process2dPippPeakTable args {

    set fname      [requiredFlagVal $args -fileName]
    set pot        [requiredFlagVal $args -pot]
    set H1name     [requiredFlagVal $args -fromProtonColumnName]
    set H2name     [requiredFlagVal $args -toProtonColumnName]
    set idName     [flagVal $args -peakIDcolumnName "PkID"]
    set intName    [flagVal $args -intensityColumnName "Intensity"]
    set namePrefix [flagVal $args -namePrefix "3d"]
    set remVar     [flagVal $args -remarksVariableName ""]
    set intenThresh [flagVal $args -intensityThreshold 0]
     
    set peaks [readPippPeakTable -fileName $fname -cols \
	       [list $idName $intName $H1name $H2name]]
    set nPeaksAdded 0

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
	$tempPeak appendToNote [format "from file %s, peak %s" $fname $curID]

	# Testing (14-Jun-2019, GAB)	
	#puts [format "%s %s %s %s %s" $curName [lindex $peak 0] [lindex $peak 2] [lindex $peak 3] [lindex $peak 1]]
	
	if {[expr abs([$tempPeak intensity])] >= $intenThresh} {
	     $pot addPeak [$tempPeak cget -this]
	     incr nPeaksAdded
	}
	$tempPeak -disown
	rename $tempPeak ""
    }

    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "Created %d peaks from file %s" $nPeaksAdded $fname]
    }
    return ""
}
 
proc process3dPippPeakTable args {

    set fname      [requiredFlagVal $args -fileName]
    set pot        [requiredFlagVal $args -pot]
    set H1name     [requiredFlagVal $args -fromProtonColumnName]
    set C1name     [requiredFlagVal $args [list -fromHeavyatomColumnName -fromCarbonColumnName -fromNitrogenColumnName]]
    set H2name     [requiredFlagVal $args -toProtonColumnName]
    set idName     [flagVal $args -peakIDcolumnName "PkID"]
    set intName    [flagVal $args -intensityColumnName "Intensity"]
    set namePrefix [flagVal $args -namePrefix "3d"]
    set remVar     [flagVal $args -remarksVariableName ""]
     
    set peaks [readPippPeakTable -fileName $fname \
		   -cols [list $idName $intName $H1name $H2name $C1name]]
    set nPeaksAdded 0

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
	$tempPeak appendToNote [format "from file %s, peak %s" $fname $curID]
	
	$pot addPeak [$tempPeak cget -this]
	$tempPeak -disown
	rename $tempPeak ""
	incr nPeaksAdded
    }

    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "Created %d peaks from file %s" $nPeaksAdded $fname]
    }
    return ""
}
 
proc process4dPippPeakTable args {

    set fname      [requiredFlagVal $args -fileName]
    set pot        [requiredFlagVal $args -pot]
    set H1name     [requiredFlagVal $args -fromProtonColumnName]
    set C1name     [requiredFlagVal $args [list -fromHeavyatomColumnName -fromCarbonColumnName -fromNitrogenColumnName]]
    set H2name     [requiredFlagVal $args -toProtonColumnName]
    set C2name     [requiredFlagVal $args [list -toHeavyatomColumnName -toCarbonColumnName -toNitrogenColumnName]]
    set idName     [flagVal $args -peakIDcolumnName "PkID"]
    set intName    [flagVal $args -intensityColumnName "Intensity"]
    set namePrefix [flagVal $args -namePrefix "4d"]
    set remVar     [flagVal $args -remarksVariableName ""]
     
    set peaks [readPippPeakTable -fileName $fname -cols [list $idName $intName $H1name $H2name $C1name $C2name]]

    foreach peak $peaks {

	set curID   [lindex $peak 0]
	set curName [format "%s%s" $namePrefix $curID]
	set tempPeak [Peak -args $curName]
	$tempPeak setIntensity [lindex $peak 1]
	$tempPeak setFromProtonShift [lindex $peak 2]
	$tempPeak setToProtonShift [lindex $peak 3]
	$tempPeak setFromHeavyatomShift [lindex $peak 4]
	$tempPeak setToHeavyatomShift [lindex $peak 5]
	$tempPeak appendToNote [format "from file %s, peak %s" $fname $curID]

	$pot addPeak [$tempPeak cget -this]
	$tempPeak -disown
	rename $tempPeak ""
    }

    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "Created %d peaks from file %s" [llength $peaks] $fname]
    }
    return ""
}

# 
# these are here just for back-compatability with scripts that call the 
# experiment-specific proc
#


proc process3dNPippPeakTable args {

    set remVar [flagVal $args -remarksVariableName ""]
    if {$remVar != ""} {
	upvar $remVar $remVar
    }

    if {! [flagExists $args -namePrefix]} {
	lappend args "-namePrefix" "3dn"
    }

    eval process3dPippPeakTable $args
    return ""
}

proc process3dCPippPeakTable args {

    set remVar [flagVal $args -remarksVariableName ""]
    if {$remVar != ""} {
	upvar $remVar $remVar
    }

    if {! [flagExists $args -namePrefix]} {
	lappend args "-namePrefix" "3dc"
    }

    eval process3dPippPeakTable $args
    return ""
}
	
proc process4dNCPippPeakTable args {

    set remVar [flagVal $args -remarksVariableName ""]
    if {$remVar != ""} {
	upvar $remVar $remVar
    }

    if {! [flagExists $args -namePrefix]} {
	lappend args "-namePrefix" "4dnc"
    }

    eval process4dPippPeakTable $args
    return ""
}

proc process4dCCPippPeakTable args {

    set remVar [flagVal $args -remarksVariableName ""]
    if {$remVar != ""} {
	upvar $remVar $remVar
    }

    if {! [flagExists $args -namePrefix]} {
	lappend args "-namePrefix" "4dcc"
    }

    eval process4dPippPeakTable $args
    return ""
}


proc writePippPeakTable args {

    set fileName       [requiredFlagVal $args -fileName]
    set namePrefix     [flagVal $args -namePrefix]
    set pot            [flagVal $args -pot]
    set peakList       [flagVal $args -peakList]
    set idName         [flagVal $args -peakIDcolumnName "PkID"]
    set intensityName  [flagVal $args -intensityColumnName "Intensity"]
    set fromProtonName [flagVal $args -fromProtonColumnName]
    set fromHeavyName  [flagVal $args -fromHeavyatomColumnName]
    set toProtonName   [flagVal $args -toProtonColumnName]
    set toHeavyName    [flagVal $args -toHeavyatomColumnName]
    set remarksList    [flagVal $args -remarks]

    if {($pot == "") && ($peakList == "")} {
	error "Neither -pot or -peakList defined in call to writeMarvinPeaks"
    }
    
    if {$pot != ""} {
	set peakList [$pot peaks]
    }
    

    if {[catch {set outUnit [open $fileName w]}]} {
	error "Error opening output file $fileName"
    }

    #
    # write remarks
    #

    puts $outUnit [format "\#"]
    puts $outUnit [format "\# Writing %d peaks to file %s" [llength $peakList] $fileName]
    puts $outUnit [format "\#"]
    
    if {[llength $remarksList] > 0} {
	puts $outUnit [format "\#"]
	foreach rem $remarksList {
	    foreach l [split $rem "\n"] {
		puts $outUnit [format "\# %s" $l]
	    }
	puts $outUnit [format "\#"]
	}
    }


    #
    # write the FORMAT and VARS lines
    #

    set formatLine "FORMAT"
    set varsLine   "VARS  "

    foreach colName [list $idName $fromProtonName $fromHeavyName $toProtonName $toHeavyName $intensityName] \
	colFormat [list "%4d" "%10.4f" "%10.4f" "%10.4f" "%10.4f" "%+8.2e"]  {

	    if {$colName != ""} {

		set formatLine [format "%s   %s" $formatLine $colFormat]
		set varsLine   [format "%s      %s" $varsLine   $colName]
	    }
	}

    puts $outUnit $formatLine
    puts $outUnit $varsLine

    #
    # write the individual peak data lines
    #

    foreach curPeak $peakList {
	Peak -this $curPeak 

	set curLine "      "

	if {$idName != ""} {
	    scan [$curPeak name] "$namePrefix%d" curPeakID
	    set curLine [format "%s %4d" $curLine $curPeakID]
	}

	if {$fromProtonName != ""} {
	    set curLine [format "%s %10.4f" $curLine [$curPeak fromProtonShift]]
	}

	if {$fromHeavyName != ""} {
	    set curLine [format "%s %10.4f" $curLine [$curPeak fromHeavyatomShift]]
	}

	if {$toProtonName != ""} {
	    set curLine [format "%s %10.4f" $curLine [$curPeak toProtonShift]]
	}

	if {$toHeavyName != ""} {
	    set curLine [format "%s %10.4f" $curLine [$curPeak toHeavyatomShift]]
	}

	if {$intensityName != ""} {
	    set curLine [format "%s %+8.2e" $curLine [$curPeak intensity]]
	}

	puts $outUnit $curLine

	rename $curPeak ""
    }

    close $outUnit

    return ""
}
