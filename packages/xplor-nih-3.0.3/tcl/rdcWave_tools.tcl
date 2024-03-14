#
# procedures to read and write RDCWave restraints
#
# JJK 2/3/04
#

package provide rdcwave 1.0
package require rdcwavepot
package require marvin

proc readRDCWaveRestraints args {

    global errorInfo
    
    set fileName    [requiredFlagVal $args -fileName]
    set potential   [requiredFlagVal $args -pot]
    set defPhiWid   [flagVal $args -defaultPhiWidth 0.0000000000]
    set defThetaWid [flagVal $args -defaultThetaWidth 0.000000000]
    
    if { [catch {set inUnit [open $fileName r]}] } {
    	error "Error opening input file $fileName"
    }
    
    XplorCommand "set message off echo off end"
    
    while {![eof $inUnit]} {
	
	if { [catch {set curRes [readOneRDCWaveRestraint $inUnit $defPhiWid $defThetaWid]} msg]} {

	    if {$msg == "No more restraints"} {
		set errorInfo ""
		break
	    } else {
		error [format "Restraint reading error in file %s:\n %s" $fileName $msg]
	    }
	}

	processOneRDCWaveRestraint $curRes $potential
	updateUser [format "read %s\r" [lindex $curRes 0]]
    }
    
    close $inUnit
}

proc readOneRDCWaveRestraint {fileID defPhiWid defThetaWid} {

    global errorInfo

    set Pi 3.1415926535897932384626433832795028841971693993

    # 
    # Try to find the beginning of another restraint.  
    #
    
    if {[catch {eatNormalTextUpTo $fileID "restraint"}]} {
	
	set errorInfo ""
	error "No more restraints"
    }

    # defaults for optional fields 
    
    set phiTarget 0.0  
    set thetaTarget 0.0
    set phiWidth $defPhiWid
    set thetaWidth $defThetaWid
    set curMode "none"

    # grab the required fields

    set curName [nextWord $fileID]
    if {$curName == ""} {
	error "Missing restraint name"
    }
    
    # search for optional fields

    set done 0
    while {! $done} {
	
	set curWord [nextWord $fileID]
	
	switch -exact -- $curWord {
	    
	    "" {
		error [format "unexpected end of file at restraint %s" $curName]
	    }
	    
	    "end" {
		set done 1
	    }
	    
	    "restraint" {
		error [format "missing end statement after restraint %s" $curName]
	    }
	    
	    "-phiTarget" {

		set phiTarget [nextWord $fileID]
		
		if {! [isNumber $phiTarget]} {
		    error [format "Missing or invalid phiTarget for restraint %s" $curName]
		}

		if {($phiTarget > $Pi) || ($phiTarget < [expr -1 * $Pi])} {
		    if {$phiTarget <= [expr 2 * $Pi]} {
			puts stderr [format "Warning: phiTarget appears to be in range 0..2Pi for restraint %s" $curName]
			set phiTarget [expr $phiTarget - (2 * $Pi)]
		    } else {
			error [format "Warning: phiTarget not in range -Pi..Pi for restraint %s" $curName]
		    }
		}

	    }

	    "-thetaTarget" {

		set thetaTarget [nextWord $fileID]
		
		if {! [isNumber $thetaTarget]} {
		    error [format "Missing or invalid thetaTarget for restraint %s" $curName]
		}

		if {($thetaTarget > $Pi) || ($thetaTarget < 0)} {
		    error [format "Warning: thetaTarget not in range 0..Pi for restraint %s" $curName]
		}

	    }

	    "-phiWidth" {

		set phiWidth [nextWord $fileID]
		
		if {! [isNumber $phiWidth]} {
		    error [format "Missing or invalid phiWidth for restraint %s" $curName]
		}

		if {($phiWidth < 0) || ($phiWidth > [expr 2 * $Pi])} {
		    error [format "Warning: phiWidth not in range 0..2*Pi for restraint %s" $curName]
		}

	    }

	    "-thetaWidth" {

		set thetaWidth [nextWord $fileID]
		
		if {! [isNumber $thetaWidth]} {
		    error [format "Missing or invalid thetaWidth for restraint %s" $curName]
		}

		if {($thetaWidth < 0) || ($thetaWidth > $Pi)} {
		    error [format "Warning: thetaWidth not in range 0..Pi for restraint %s" $curName]
		}

	    }



	    "-simple" {
		
		set curMode "simple"

		set curSelA [nextSel $fileID]
		if {$curSelA == ""} {
		    error [format "Missing or invalid 1st selection for simple restraint %s" $curName]
		}

		set curSelB [nextSel $fileID]
		if {$curSelB == ""} {
		    error [format "Missing or invalid 2nd selection for simple restraint %s" $curName]
		}

		set curSelC ""
		set curSelD ""
	    }

	    "-norm" {
		
		set curMode "norm"

		set curSelA [nextSel $fileID]
		if {$curSelA == ""} {
		    error [format "Missing or invalid 1st selection for norm restraint %s" $curName]
		}

		set curSelB [nextSel $fileID]
		if {$curSelB == ""} {
		    error [format "Missing or invalid 2nd selection for norm restraint %s" $curName]
		}

		set curSelC [nextSel $fileID]
		if {$curSelC == ""} {
		    error [format "Missing or invalid 3rd selection for norm restraint %s" $curName]
		}

		set curSelD [nextSel $fileID]
		if {$curSelD == ""} {
		    error [format "Missing or invalid 4th selection for norm restraint %s" $curName]
		}
	    }
	}
    }
    
    #
    # make sure we set a mode
    #

    if {$curMode == "none"} {
	error [format "Missing -simple or -norm flag for restraint %s" $curName]
    }

    return [list $curName $phiTarget $thetaTarget $phiWidth $thetaWidth $curMode $curSelA $curSelB $curSelC $curSelD]
}


proc processOneRDCWaveRestraint {restraintData aPotential} {
    
    #
    # takes the data read by readOneRDCWaveRestraint
    # and actually creates an RDCWave restraint
    # and attaches it to the given RDCWavePot
    #
    
    set curName        [lindex $restraintData 0]
    set curPhiTarget   [lindex $restraintData 1]
    set curThetaTarget [lindex $restraintData 2]
    set curPhiWidth    [lindex $restraintData 3]
    set curThetaWidth  [lindex $restraintData 4]
    set curMode        [lindex $restraintData 5]
    set curSelA        [lindex $restraintData 6]
    set curSelB        [lindex $restraintData 7]
    set curSelC        [lindex $restraintData 8]
    set curSelD        [lindex $restraintData 9]
    
    #
    # create the restraint
    # note that these selections are attached to the currentSimulation
    #

    RDCWaveRestraint newRestraint $curName

    newRestraint setPhiTarget   $curPhiTarget
    newRestraint setThetaTarget $curThetaTarget
    newRestraint setPhiWidth    $curPhiWidth
    newRestraint setThetaWidth  $curThetaWidth

    if {$curMode == "simple"} {

	newRestraint makeSimpleRestraint
	
	AtomSel temp $curSelA
	newRestraint setAtomA [temp cget -this]
	temp -disown
	rename temp ""

	AtomSel temp $curSelB
	newRestraint setAtomB [temp cget -this]
	temp -disown
	rename temp ""

    } else {

	newRestraint makeNormRestraint

	AtomSel temp $curSelA
	newRestraint setAtomA [temp cget -this]
	temp -disown
	rename temp ""

	AtomSel temp $curSelB
	newRestraint setAtomB [temp cget -this]
	temp -disown
	rename temp ""

	AtomSel temp $curSelC
	newRestraint setAtomC [temp cget -this]
	temp -disown
	rename temp ""

	AtomSel temp $curSelD
	newRestraint setAtomD [temp cget -this]
	temp -disown
	rename temp ""
    }
    
    #
    # add the restraint to the specified RDCWave potential
    #
    
    $aPotential addToRestraints [newRestraint cget -this]
    
    newRestraint -disown
    rename newRestraint ""
}

proc stringAppend {s1name s2} {

    upvar $s1name s1

    set s1 [format "%s\n%s" $s1 $s2]
}



proc RDCWaveReport args {
    
    set retVal ""

    set pot [requiredFlagVal $args -pot]

    set resPtrList [$pot restraints]

    set o [$pot atomO]
    set x [$pot atomX]
    set y [$pot atomY]
    set z [$pot atomZ]

    Atom -this $o
    stringAppend retVal  [format "%s %d %s %s %s" [$o segmentName] [$o residueNum] \
			      [$o residueName] [$o atomName] [$o pos]]
    Atom -this $x
    stringAppend retVal  [format "%s %d %s %s %s" [$x segmentName] [$x residueNum] \
			      [$x residueName] [$x atomName] [$x pos]]
    Atom -this $y
    stringAppend retVal  [format "%s %d %s %s %s" [$y segmentName] [$y residueNum] \
			      [$y residueName] [$y atomName] [$y pos]]
    Atom -this $z
    stringAppend retVal  [format "%s %d %s %s %s" [$z segmentName] [$z residueNum] \
			      [$z residueName] [$z atomName] [$z pos]]
    rename $o ""
    rename $x ""
    rename $y ""
    rename $z ""

    stringAppend retVal [format "energy: %f" [$pot calcEnergy]]

    stringAppend retVal \
 "Restraint name   targTheta calcTheta deltaTheta --- targPhi calcPhi deltaPhi"

    foreach rPtr $resPtrList {
	
	RDCWaveRestraint temp -this $rPtr

	set targPhi  [temp phiTarget]
	set curPhi   [temp calcPhi  $o $x $y $z]
	set deltaPhi [temp deltaPhi $o $x $y $z]

	set targTheta  [temp thetaTarget]
	set curTheta   [temp calcTheta $o $x $y $z]
	set deltaTheta [temp deltaTheta $o $x $y $z]

	stringAppend retVal [format "%s   %5.3f %5.3f %5.3f  --- %5.3f %5.3f %5.3f " \
				 [temp name] $targTheta $curTheta $deltaTheta \
				 $targPhi $curPhi $deltaPhi]
	
	rename temp ""
    }

    return $retVal
}



#
# for a given RDCWavePot, 
# adds axis atoms to the simulation,
# modifies the IVM to allow their proper rotation,
# and sets the axis atoms in the RDCWavePot
#

proc setupRDCWaveAxes args {

    set pot [requiredFlagVal $args -pot]

    namespace eval MarvinDefaults {
	variable IVMObjectName
    }

    set ivm $MarvinDefaults::IVMObjectName

    marvinPyth command "import varTensorTools"
    marvinPyth command "from varTensorTools import addAxisAtoms"

    set outVars [marvinPyth command "sr = addAxisAtoms()" {} {sr}]

    #
    # check the Python variables I extracted for one
    # named 'sr'
    #
    
    foreach item $outVars {
	
	set name [lindex $item 0]
	set val  [lindex $item 1]
	
	if {$name == "sr"} {
	    scan $val  "(%s %d)" segName resNum
	    set segName [string range $segName 1 end-2]
	}
    }

    #
    # delete the OO2, PA1, and PA2 axis atoms
    #

    XplorCommand "dele sele (segid $segName and resid $resNum and 
                             name OO2 or name PA1 or name PA2) end"

    #
    # grab the axes atoms
    #

    AtomSel oSel "segid $segName and resid $resNum and name OO"
    AtomSel xSel "segid $segName and resid $resNum and name X"
    AtomSel ySel "segid $segName and resid $resNum and name Y"
    AtomSel zSel "segid $segName and resid $resNum and name Z"

    set o [oSel atomAtPos 0]
    set x [xSel atomAtPos 0]
    set y [ySel atomAtPos 0]
    set z [zSel atomAtPos 0]

    Atom -this $o
    Atom -this $x
    Atom -this $y
    Atom -this $z

    set oID   [$o index]
    set xID   [$x index]
    set yID   [$y index]
    set zID   [$z index]

    rename $o ""
    rename $x ""
    rename $y ""
    rename $z ""

    #
    # set their masses to that their sum is equal to that of the 
    # rest of the system
    #

    XplorCommand "vector show sum (mass) (not (segid $segName and resid $resNum))"
    set totMass [XplorVariableNamed "result"]
    set newMass [expr $totMass / 4.0]
    XplorCommand "vector do (mass = $newMass) (segid $segName and resid $resNum)"

    #
    # break the ox and oz bonds
    #

    marvinPyth command "$ivm.setBondExclude( $ivm.bondExclude() + \[($xID, $oID)\] )"
    marvinPyth command "$ivm.setBondExclude( $ivm.bondExclude() + \[($zID, $oID)\] )"
    
    #
    # set up the rigid grouping
    #
    
    marvinPyth command "$ivm.group( ($xID, $yID, $zID, $oID) )"  
    
    #
    # add a new base atom 
    #
    
    # marvinPyth command "$ivm.setBaseAtoms( $ivm.baseAtoms() + \[$zID,\] )"

    #
    # add a new hinge
    #

    marvinPyth command "$ivm.setHingeList( $ivm.hingeList() + \[('rotate', ($xID, $yID, $zID, $oID))\] )"

    marvinPyth command "$ivm.autoTorsion()"

    #
    # tell the RDCWavePot what the axes atoms are
    #
    
    $pot setAlignmentTensorAtoms \
	[oSel cget -this] \
	[xSel cget -this] \
	[ySel cget -this] \
	[zSel cget -this] 

}
