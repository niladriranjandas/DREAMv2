#
# routines to read TALOS output into xplor-nih
#

package provide marvin 1.0


proc readTalosN args {

    set predFile     [flagVal $args -predFile pred.tab]
    set predAllFile  [flagVal $args -predAllFile predAll.tab]
    set outName      [flagVal $args -outFileName]
    set scale        [flagVal $args -scale 200.0]
    set widthPadding [flagVal $args -widthPadding 5.0]
    set minWidth     [flagVal $args -minWidth 20.0]
    set segNames     [flagVal $args -segName [list ""]]	
    set useWeights   [flagVal $args -useWeights True]
    set omitResids   [flagVal $args -omitResids [list]]	

    marvinPyth command "from talosTools import convertTalosN"

    set string ""
    foreach segid $segNames {
	set cmd "str=convertTalosN('$predFile',
				   '$predAllFile',
				   $widthPadding,
				   $minWidth,
				   segName='$segid',
                                   useWeights=$useWeights,
                                   omitResids=omitResids.split())"
	set outVars [marvinPyth command $cmd \
			 [list [list omitResids $omitResids]] {{str}}]
	foreach item $outVars {
	    
	    set name [lindex $item 0]
	    set val  [lindex $item 1]
	    if {$name == "str"} {
		set string "$string $val \n"
	    }
	}
    }

    if {$outName == ""} {
	marvinPyth command "import protocol"
	set cmd [format "protocol.initDihedrals(string=str,
                                                scale=%f,
                                                useDefaults=False,
                                                )" $scale]
	marvinPyth command $cmd [list [list str $string]]
    } else {
	# 
	# write out the corresponding dihedral restraints
	#

	if {[catch {set outUnit [open $outName w]}]} {
	    error "Error opening output file $fname"
	}

	puts outUnit $str
	close $outUnit
    }
}
  

proc readIndividualResidueTalos args {

    set fnames       [requiredFlagVal $args -fileNames]
    set outName      [flagVal $args -outFileName]
    set scale        [flagVal $args -scale 200.0]
    set widthPadding [flagVal $args -widthPadding 5.0]
    set minWidth     [flagVal $args -minWidth 20.0]
    set segNames     [flagVal $args -segName [list ""]]	
    
    #
    # Reads the residue-by-residue files (typically pred/res001.tab, etc),
    # Defines angle ranges using the maximum and minimum values found among the database hits, 
    # Pads that range by the given amount on either end
    # Makes sure that the range is at least minWdith wide in each dimension
    # Creates restraints for each segment specified
    #

    global errorInfo

    #
    # for each file,
    #

    foreach fname $fnames {
	
	if {[catch {set inUnit [open $fname r]}]} {
	    error "Error opening input file $fname"
	}

	set header         [readIndividualResidueTalosHeader $inUnit]
	set curResID       [lindex $header 0]
	set columnPosnList [lindex $header 1]

	#
	# now read each database hit line
	#

	set hitList [list]

	while {![eof $inUnit]} {
	    
	    if {![catch {set curHit [readTalosLine $inUnit $columnPosnList]} msg]} {
		
		lappend hitList $curHit
		
	    } else {
		
		set savedInfo $errorInfo
		if {$msg == "blank line"} {
		    set errorInfo ""
		} else {
		    error $msg $savedInfo
		}
	    }
	}
	
	close $inUnit

	#
	# Create phi and psi bounds from the database hits.  
	# Use the overall range of hits along each dimension, plus the 
	# width padding.
	# 

	set phiList [list]
	set psiList [list]

	foreach hit $hitList {
	    set phiHit [lindex $hit 0]
	    set psiHit [lindex $hit 1]
	    set weightHit [lindex $hit 2]

	    if {$weightHit != 0} {
		lappend phiList $phiHit
		lappend psiList $psiHit
	    }
	}

	# 
	# handle wraparounds by picking the first talos hit
	# to be the target, and processing every other talos 
	# hit so that its angular difference is < 180 degrees
	#

	if {[llength $phiList] > 0} {

	    set targ [lindex $phiList 0]
	    set wrappedAngleList [list] 
	    foreach unwrappedAngle $phiList {
		set delta1 [expr abs($targ - $unwrappedAngle)]
		set delta2 [expr abs($targ - ($unwrappedAngle - 360))]
		set delta3 [expr abs($targ - ($unwrappedAngle + 360))]
		set minDelta [min $delta1 $delta2 $delta3]
		if {$delta1 == $minDelta} {
		    lappend wrappedAngleList $unwrappedAngle
		} elseif {$delta2 == $minDelta} {
		    lappend wrappedAngleList [expr $unwrappedAngle - 360]
		} else {
		    lappend wrappedAngleList [expr $unwrappedAngle + 360]
		}
	    }

	    set minPhi [expr [min $wrappedAngleList] - $widthPadding]
	    set maxPhi [expr [max $wrappedAngleList] + $widthPadding]

	    set targ [lindex $psiList 0]
	    set wrappedAngleList [list] 
	    foreach unwrappedAngle $psiList {
		set delta1 [expr abs($targ - $unwrappedAngle)]
		set delta2 [expr abs($targ - ($unwrappedAngle - 360))]
		set delta3 [expr abs($targ - ($unwrappedAngle + 360))]
		set minDelta [min $delta1 $delta2 $delta3]
		if {$delta1 == $minDelta} {
		    lappend wrappedAngleList $unwrappedAngle
		} elseif {$delta2 == $minDelta} {
		    lappend wrappedAngleList [expr $unwrappedAngle - 360]
		} else {
		    lappend wrappedAngleList [expr $unwrappedAngle + 360]
		}
	    }

	    set minPsi [expr [min $wrappedAngleList] - $widthPadding]
	    set maxPsi [expr [max $wrappedAngleList] + $widthPadding]

	    set phiWidth  [expr ($maxPhi - $minPhi) * 0.5]
	    set psiWidth  [expr ($maxPsi - $minPsi) * 0.5]
	    set phiTarget [expr $minPhi + $phiWidth]
	    set psiTarget [expr $minPsi + $psiWidth]

	    set phiWidth [max $phiWidth $minWidth]
	    set psiWidth [max $psiWidth $minWidth]

	    puts [list $curResID $phiTarget $psiTarget $phiWidth $psiWidth]
	    
	    lappend restraints [list $curResID $phiTarget $psiTarget $phiWidth $psiWidth "Good"]
	}
    }

    if {$outName == ""} {

	#
	# create the corresponding dihedral restraints
	#
	
	XplorCommand "set message off echo off end"
	
	#
	# allocate room for the restraints
	#
	
	set nSegs [max [llength $segNames] 1]
	set nRestraints [expr [llength $restraints] * 2 * $nSegs]
	XplorCommand "restraints dihedral nass $nRestraints end"
	
	#
	# set the overall force constant scale factor
	#
	
	XplorCommand "restraints dihedral scale $scale end"
	
	#
	# create each restraint
	#
	
	set nTalosRestraints 0
	
	foreach segName $segNames {
	    foreach curRes $restraints {
		incr nTalosRestraints [createTalosRestraint 1.0 $curRes $segName]
	    }
	}
	
	updateUser [format "%d residues got phi/psi restraints from TALOS file $fname\n" $nTalosRestraints $fname]

    } else {
	
	# 
	# write out the corresponding dihedral restraints
	#

	if {[catch {set outUnit [open $outName w]}]} {
	    error "Error opening output file $fname"
	}

	set nTalosRestraints 0

	foreach segName $segNames {
	    foreach curRes $restraints {
		incr nTalosRestraints [writeTalosRestraint 1.0 $curRes $segName $outUnit]
	    }
	}

	close $outUnit

	updateUser [format "%d residues got phi/psi restraints from TALOS file %s written to %s \n" $nTalosRestraints $fname $outName]

    }
}

	    

proc readIndividualResidueTalosHeader {inUnit} {

    global errorInfo

    # eat text until we get to a line that begins DATA RESIDS

    while {1} {

	if {[catch {set l [nextLine $inUnit]}]} {
	    set errorInfo ""
	    error "Can't find DATA RESIDS line in Talos file"
	}
	
	if {([lindex $l 0] == "DATA") && ([lindex $l 1] == "RESIDS")} {
	    break
	}
    }

    # record the residue number for the central residue 

    set curResid [lindex $l 3]

    # rewind to beginning of file

    seek $inUnit 0

    # eat text until we get to a line that begins with VARS

    while {1} {

	if {[catch {set l [nextLine $inUnit]}]} {
	    set errorInfo ""
	    error "Can't find VARS line in Talos file"
	}
	
	if {[lindex $l 0] == "VARS"} {
	    break
	}
    }
   	
    # grab the column numbers for phi, psi, and weight

    set phiPos     [expr [lsearch -exact $l "PHI"] - 1]
    set psiPos     [expr [lsearch -exact $l "PSI"] - 1]
    set weightPos  [expr [lsearch -exact $l "W"] - 1]

    if {$phiPos == -2} {
	error "PHI column not found in Talos VARS line"
    }

    if {$psiPos == -2} {
	error "PSI column not found in Talos VARS line"
    }

    if {$weightPos == -2} {
	error "W column not found in Talos VARS line"
    }

    return [list $curResid [list $phiPos $psiPos $weightPos]]
}


proc readTalos args {

    set fname       [requiredFlagVal $args -fileName]
    set scale       [flagVal $args -scale 200.0]
    set widthScale  [flagVal $args -widthScale 1]
    set segNameList [flagVal $args -segName [list ""]]	

    global errorInfo

    if {[catch {set inUnit [open $fname r]}]} {
	error "Error opening input file $fname"
    }

    set columnPosnList [readTalosHeader $inUnit]

    set resList [list]

    #
    # now read each data line
    #

    while {![eof $inUnit]} {

	if {![catch {set curRes [readTalosLine $inUnit $columnPosnList]} msg]} {

	    lappend resList $curRes

	} else {

	    set savedInfo $errorInfo
	    if {$msg == "blank line"} {
		set errorInfo ""
	    } else {
		error $msg $savedInfo
	    }
	}
    }

    close $inUnit

    XplorCommand "set message off echo off end"

    #
    # allocate room for the restraints
    #

    set nSegs [max [llength $segNameList] 1]
    set nRestraints [expr [llength $resList] * 2 * $nSegs]
    XplorCommand "restraints dihedral nass $nRestraints end"

    #
    # set the overall force constant scale factor
    #

    XplorCommand "restraints dihedral scale $scale end"

    #
    # create each restraint
    #

    set nTalosRestraints 0

    foreach segName $segNameList {
	foreach curRes $resList {
	    incr nTalosRestraints [createTalosRestraint $widthScale $curRes $segName]
	}
    }

    updateUser [format "%d residues got phi/psi restraints from TALOS file $fname\n" $nTalosRestraints $fname]
}


proc readTalosHeader {inUnit} {

    global errorInfo

    # eat text until we get to a line that begins with VARS

    while {1} {

	if {[catch {set l [nextLine $inUnit]}]} {
	    set errorInfo ""
	    error "Can't find VARS line in Talos file"
	}
	
	if {[lindex $l 0] == "VARS"} {
	    break
	}
    }
   	
    # grab the column numbers for resid, phi, psi, phierror, psierror, and class

    set residPos   [expr [lsearch -exact $l "RESID"] - 1]
    set phiPos     [expr [lsearch -exact $l "PHI"] - 1]
    set psiPos     [expr [lsearch -exact $l "PSI"] - 1]
    set dPhiPos    [expr [lsearch -exact $l "DPHI"] - 1]
    set dPsiPos    [expr [lsearch -exact $l "DPSI"] - 1]
    set classPos   [expr [lsearch -exact $l "CLASS"] - 1]

    if {$residPos == -2} {
	error "RESID column not found in Talos VARS line"
    } 

    if {$phiPos == -2} {
	error "PHI column not found in Talos VARS line"
    }

    if {$psiPos == -2} {
	error "PSI column not found in Talos VARS line"
    }

    if {$dPhiPos == -2} {
	error "DPHI column not found in Talos VARS line"
    }

    if {$dPsiPos == -2} {
	error "DPSI column not found in Talos VARS line"
    }

    if {$classPos == -2} {
	error "CLASS column not found in Talos VARS line"
    }

    return [list $residPos $phiPos $psiPos $dPhiPos $dPsiPos $classPos]
}


proc readTalosLine {inUnit columnPosnList} {

    global errorInfo

    #
    # grab the line
    #

    set l [nextLine $inUnit]

    if {($l == "") || ([string index $l 0] == "\#") || ([lindex $l 0] == "FORMAT")} {
	error "blank line"
    }

    #
    # try extracting out the columns' data
    #

    if {[catch {

	set resid   [lindex $l [lindex $columnPosnList 0]]
	set phi     [lindex $l [lindex $columnPosnList 1]]
	set psi     [lindex $l [lindex $columnPosnList 2]]
	set dPhi    [lindex $l [lindex $columnPosnList 3]]
	set dPsi    [lindex $l [lindex $columnPosnList 4]]
	set class   [lindex $l [lindex $columnPosnList 5]]

    }]} {
	
	set savedInfo $errorInfo
	error "can't parse line: $l" $savedInfo
    } 

    #
    # return the columns' data
    #

    return [list $resid $phi $psi $dPhi $dPsi $class]
}

proc createTalosRestraint {widthScale curRes segName} {
    
    set resid   [lindex $curRes 0]
    set phi     [lindex $curRes 1]
    set psi     [lindex $curRes 2]
    set dPhi    [lindex $curRes 3]
    set dPsi    [lindex $curRes 4]
    set class   [lindex $curRes 5]

    if {$class == "Good"} {

	set prevRes [expr $resid - 1]
	set curRes $resid 
	set nextRes [expr $resid + 1]

	# create the phi restraint

	if {$segName == ""} {
	    set s1 [format "(resid %d and name c)"  $prevRes]
	    set s2 [format "(resid %d and name n)"  $curRes]
	    set s3 [format "(resid %d and name ca)" $curRes]
	    set s4 [format "(resid %d and name c)"  $curRes]
	} else {
	    set s1 [format "(segid %s and resid %d and name c)"  $segName $prevRes]
	    set s2 [format "(segid %s and resid %d and name n)"  $segName $curRes]
	    set s3 [format "(segid %s and resid %d and name ca)" $segName $curRes]
	    set s4 [format "(segid %s and resid %d and name c)"  $segName $curRes]
	}

	set cmd [format "restraints dihedral assign %s\n %s\n %s\n %s\n 1.0 %f %f 2 end" \
		     $s1 $s2 $s3 $s4 $phi [expr $dPhi * $widthScale]]

	XplorCommand $cmd

	# create the psi restraint

	if {$segName == ""} {
	    set s1 [format "(resid %d and name n)"  $curRes]
	    set s2 [format "(resid %d and name ca)" $curRes]
	    set s3 [format "(resid %d and name c)"  $curRes]
	    set s4 [format "(resid %d and name n)"  $nextRes]
	} else {
	    set s1 [format "(segid %s and resid %d and name n)"  $segName $curRes]
	    set s2 [format "(segid %s and resid %d and name ca)" $segName $curRes]
	    set s3 [format "(segid %s and resid %d and name c)"  $segName $curRes]
	    set s4 [format "(segid %s and resid %d and name n)"  $segName $nextRes]
	}

	set cmd [format "restraints dihedral assign %s\n %s\n %s\n %s\n 1.0 %f %f 2 end" \
		     $s1 $s2 $s3 $s4 $psi [expr $dPsi * $widthScale]]

	XplorCommand $cmd

	return 1
    } else {
	return 0
    }
}


proc writeTalosRestraint {widthScale curRes segName outUnit} {
    
    set resid   [lindex $curRes 0]
    set phi     [lindex $curRes 1]
    set psi     [lindex $curRes 2]
    set dPhi    [lindex $curRes 3]
    set dPsi    [lindex $curRes 4]
    set class   [lindex $curRes 5]

    if {$class == "Good"} {

	set prevRes [expr $resid - 1]
	set curRes $resid 
	set nextRes [expr $resid + 1]

	# create the phi restraint

	if {$segName == ""} {
	    set s1 [format "(resid %d and name c)"  $prevRes]
	    set s2 [format "(resid %d and name n)"  $curRes]
	    set s3 [format "(resid %d and name ca)" $curRes]
	    set s4 [format "(resid %d and name c)"  $curRes]
	} else {
	    set s1 [format "(segid %s and resid %d and name c)"  $segName $prevRes]
	    set s2 [format "(segid %s and resid %d and name n)"  $segName $curRes]
	    set s3 [format "(segid %s and resid %d and name ca)" $segName $curRes]
	    set s4 [format "(segid %s and resid %d and name c)"  $segName $curRes]
	}


	puts $outUnit [format "assign %s\n %s\n %s\n %s\n 1.0 %f %f 2" $s1 $s2 $s3 $s4 $phi [expr $dPhi * $widthScale]]


	# create the psi restraint

	if {$segName == ""} {
	    set s1 [format "(resid %d and name n)"  $curRes]
	    set s2 [format "(resid %d and name ca)" $curRes]
	    set s3 [format "(resid %d and name c)"  $curRes]
	    set s4 [format "(resid %d and name n)"  $nextRes]
	} else {
	    set s1 [format "(segid %s and resid %d and name n)"  $segName $curRes]
	    set s2 [format "(segid %s and resid %d and name ca)" $segName $curRes]
	    set s3 [format "(segid %s and resid %d and name c)"  $segName $curRes]
	    set s4 [format "(segid %s and resid %d and name n)"  $segName $nextRes]
	}

	puts $outUnit [format "assign %s\n %s\n %s\n %s\n 1.0 %f %f 2" $s1 $s2 $s3 $s4 $psi [expr $dPsi * $widthScale]]

	return 1
    } else {
	return 0
    }
}




#
# for convenience, read in a TALOS file & write it out as a dihedral restraints file
#

proc convertTalosToDihedral args {

    set inFname    [requiredFlagVal $args -inFileName]
    set outFname   [requiredFlagVal $args -outFileName]
    set scale      [flagVal $args -scale 200.0]
    set widthScale [flagVal $args -widthScale 1]
    set segName    [flagVal $args -segName]

    global errorInfo

    if {[catch {set inUnit [open $inFname r]}]} {
	error "Error opening input file $inFname"
    }

    set columnPosnList [readTalosHeader $inUnit]

    set resList [list]

    #
    # now read each data line
    #

    while {![eof $inUnit]} {

	if {![catch {set curRes [readTalosLine $inUnit $columnPosnList]} msg]} {

	    lappend resList $curRes

	} else {

	    set savedInfo $errorInfo
	    if {$msg == "blank line"} {
		set errorInfo ""
	    } else {
		error $msg $savedInfo
	    }
	}
    }

    close $inUnit


    if {[catch {set outUnit [open $outFname w]}]} {
	error "Error opening output file $outFname"
    }

    puts $outUnit [format "! Dihedral restraints file created from TALOS file %s" $inFname]

    #
    # create each restraint
    #

    set nTalosRestraints 0

    foreach curRes $resList {

	set resid   [lindex $curRes 0]
	set phi     [lindex $curRes 1]
	set psi     [lindex $curRes 2]
	set dPhi    [lindex $curRes 3]
	set dPsi    [lindex $curRes 4]
	set class   [lindex $curRes 5]

	if {$class == "Good"} {

	    set prevRes [expr $resid - 1]
	    set curRes $resid 
	    set nextRes [expr $resid + 1]

	    # create the phi restraint
	    
	    if {$segName == ""} {
		set s1 [format "(resid %d and name c)"  $prevRes]
		set s2 [format "(resid %d and name n)"  $curRes]
		set s3 [format "(resid %d and name ca)" $curRes]
		set s4 [format "(resid %d and name c)"  $curRes]
	    } else {
		set s1 [format "(segid %s and resid %d and name c)"  $segName $prevRes]
		set s2 [format "(segid %s and resid %d and name n)"  $segName $curRes]
		set s3 [format "(segid %s and resid %d and name ca)" $segName $curRes]
		set s4 [format "(segid %s and resid %d and name c)"  $segName $curRes]
	    }

	    set cmd [format "assign %s\n %s\n %s\n %s\n 1.0 %f %f 2 \n" \
			 $s1 $s2 $s3 $s4 $phi [expr $dPhi * $widthScale]]

	    puts $outUnit $cmd

	    # create the psi restraint
	    
	    if {$segName == ""} {
		set s1 [format "(resid %d and name n)"  $curRes]
		set s2 [format "(resid %d and name ca)" $curRes]
		set s3 [format "(resid %d and name c)"  $curRes]
		set s4 [format "(resid %d and name n)"  $nextRes]
	    } else {
		set s1 [format "(segid %s and resid %d and name n)"  $segName $curRes]
		set s2 [format "(segid %s and resid %d and name ca)" $segName $curRes]
		set s3 [format "(segid %s and resid %d and name c)"  $segName $curRes]
		set s4 [format "(segid %s and resid %d and name n)"  $segName $nextRes]
	    }

	    set cmd [format "assign %s\n %s\n %s\n %s\n 1.0 %f %f 2 \n" \
			 $s1 $s2 $s3 $s4 $psi [expr $dPsi * $widthScale]]

	    puts $outUnit $cmd
	}
    }

    close $outUnit
}
