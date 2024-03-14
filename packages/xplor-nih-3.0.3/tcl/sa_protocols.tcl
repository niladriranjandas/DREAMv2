#
# a standard annealing schedule for random collapse 
# from an extended strand, using Rg and Ca..Ca vdW
# energies only
#

package provide marvin 1.0

proc makeRandomCoords args {
    
    set aSel [requiredFlagVal $args -selection] 
    
    for {set count 0} {$count < [$aSel size]} {incr count} {
	
	set curAtom [$aSel atomAtPos $count]
	Atom -this $curAtom
	
	$curAtom setPos [list \
			     [expr [uniformRandom] * 10.0] \
			     [expr [uniformRandom] * 10.0] \
			     [expr [uniformRandom] * 10.0]]
	rename $curAtom ""
    }
}

#
# creates a randomized structure with good covalent geometry 
# and a large Rgyr
#

proc makeRandomExpandedStruct args {

    setupIVM

    setupRgyr \
	-restraints [list [AtomSel -args "all"]] \
	-correction 15

    set rgyr [create_XplorPot COLL]
    set bond [create_XplorPot BOND]
    set angl [create_XplorPot ANGL]
    set impr [create_XplorPot IMPR]
    
    randomCollapse \
	-bathTemp 4000 \
	-potsToUse [list $bond $angl $impr $rgyr] \
	-mini \
	-allVDW 

    [$bond __deref__] -delete
    $bond -delete
    
    [$angle __deref__] -delete
    $angle -delete
    
    [$impr __deref__] -delete
    $impr -delete
    
    [$rgyr __deref__] -delete
    $rgyr -delete
}


proc makeRandomCollapsedStruct args {

    if {[flagExists $args -setupIVM]} {
	setupIVM
    }

    setupRgyr \
	-restraints [list [AtomSel -args "all"]] \
	-correction -1.0

    set rgyr [create_XplorPot COLL]
    set bond [create_XplorPot BOND]
    set angl [create_XplorPot ANGL]
    set impr [create_XplorPot IMPR]
    
    randomCollapse \
	-bathTemp 4000 \
	-potsToUse [list $bond $angl $impr $rgyr] \
	-mini \
	-allVDW 

    [$bond __deref__] -delete
    $bond -delete
    
    [$angle __deref__] -delete
    $angle -delete
    
    [$impr __deref__] -delete
    $impr -delete
    
    [$rgyr __deref__] -delete
    $rgyr -delete
}

proc makeRandomStructure args {

    if {[flagExists $args -setupIVM]} {
	setupIVM
    }

    set bond [create_XplorPot BOND]
    set angl [create_XplorPot ANGL]
    set impr [create_XplorPot IMPR]
    
    randomCollapse \
	-bathTemp 4000 \
	-potsToUse [list $bond $angl $impr] \
	-mini \
	-allVDW 

    [$bond __deref__] -delete
    $bond -delete
    
    [$angle __deref__] -delete
    $angle -delete
    
    [$impr __deref__] -delete
    $impr -delete
}

proc makeExtendedStrand args {

    set targPhi  [flagVal $args -targetPhi -139] 
    set targPsi  [flagVal $args -targetPsi  135] 

    if {[flagExists $args -helix]} {
	set targPhi -60
	set targPsi -45
    }

    AtomSel tempSel "name ca"
    set tempSim [tempSel simulation]
    Simulation -this $tempSim

    set preResNum 0
    set prePreResNum 0

    XplorCommand "vector do (rmsd = decode(resid)) (all)"

    foreach curIndex [tempSel indices] {
	
	set curResNum [$tempSim residueNum $curIndex]
	
	if {$prePreResNum > 0} {
    
	    XplorCommand "pick dihedral 
                             (resid $prePreResNum and name c)
                             (resid $preResNum and name n)
                             (resid $preResNum and name ca)
                             (resid $preResNum and name c)
                             geom"

	    set curPhi [XplorVariableNamed result]
	    set needToMovePhi [expr $curPhi - $targPhi]

	    XplorCommand "pick dihedral 
                             (resid $preResNum and name n)
                             (resid $preResNum and name ca)
                             (resid $preResNum and name c)
                             (resid $curResNum and name n)
                             geom"
	    
    	    set curPsi [XplorVariableNamed result]
	    set needToMovePsi [expr $curPsi - $targPsi]

	    XplorCommand "coor rotate 
                             sele ((attribute rmsd < $preResNum) 
                                   or (name hn and (resid $preResNum)))
                             center (head (name n and resid $preResNum)) 
                             axis (tail = (name ca and resid $preResNum)
                                   head = (name n  and resid $preResNum))
                             $needToMovePhi
                          end"

	    XplorCommand "coor rotate 
                             sele ((attribute rmsd < $preResNum) 
                                   or (resid $preResNum and not (name c or name o)))
                             center (head (name ca and resid $preResNum)) 
                             axis (tail = (name c  and resid $preResNum)
                                   head = (name ca and resid $preResNum))
                             $needToMovePsi
                          end"
	}
	
	set prePreResNum $preResNum
	set preResNum $curResNum
    }

    rename tempSel ""
    rename $tempSim ""
}

#
# these are convenience functions for the main protocol procedures, below
#

proc randomCollapse args {

    
    #
    # grab the flags
    #
    
    set bathTemp               [flagVal $args -bathTemp 4000] 
    set timeLength             [flagVal $args -timeLength 5.0] 
    set convergedRgEner        [flagVal $args -convergedRgEner 2000.0]
    set convergedVdwPerResidue [flagVal $args -convergedVdwPerResidue 4.0]
    set stepsBetweenChecks     [flagVal $args -stepsBetweenChecks 100]
    set maxConvergeTries       [flagVal $args -maxConvergeTries 50]
    set potsToUse              [flagVal $args -potsToUse [list]]
    
    #
    # setup--note that the marvinPyth PyInterp object is created at the top of sa_tools.tcl
    #
    
    namespace eval MarvinDefaults {
	variable IVMObjectName
    }
    
    if {$MarvinDefaults::IVMObjectName == {}} {
	error "No IVM object name defined in randomCollapse.  Did you run setupIVM first?"
    }
        
    marvinPyth command "from monteCarlo import randomizeTorsions"

    #
    # determine the number of residues for eVdw/residue calculations
    #

    set tempSel [AtomSel -args "(all)"]
    set nResidues [llength [residuesInSelection $tempSel]]
    $tempSel -delete

    #
    # add covalent geometry potentials to potsToUse if they're not already there
    #
	
    set hasBond 0
    set hasAngle 0
    set hasImproper 0
    set hasVdw 0
	
    foreach elem $potsToUse {
	if {[$elem instanceName] == "BOND"} {
	    set hasBond 1
	}
	if {[$elem instanceName] == "ANGL"} {
	    set hasAngle 1
	}
	if {[$elem instanceName] == "IMPR"} {
	    set hasImproper 1
	}
	if {[$elem instanceName] == "VDW"} {
	    set hasVdw 1
	}
    }
	
    if {! $hasBond} {
	set bondPot [create_XplorPot "BOND"]
	lappend potsToUse $bondPot
    }
	
    if {! $hasAngle} {
	set anglePot [create_XplorPot "ANGL"]
	lappend potsToUse $anglePot
    }
    
    if {! $hasImproper} {
	set improperPot [create_XplorPot "IMPR"]
	lappend potsToUse $improperPot
    }
    
    if {! $hasVdw} {
	set vdwPot [create_XplorPot VDW]
	# don't append vdwPot to potsToUse--it's turned on and off separately
    }


    set successful 0

    while {! $successful} {

	#
	# randomize the torsions
	#

	marvinPyth command "randomizeTorsions($MarvinDefaults::IVMObjectName)"

	#
	# the torsion randomizer can mess up certain types of covalent geometry
	# (for example, disulfides), so re-clean them here
	#
	# If I can't get them to clean up, repeat with another set of random torsions
	#

	if {[catch cleanCovalentGeom]} {
	    continue
	}
	
	#
	# clean up everything but the vdw
	#
	
	potList removeAll
	foreach pot $potsToUse {
	    #	if {[catch {potList byName [$pot instanceName]}]} {
	    potList add [$pot cget -this]
	    #	}
	}
	
	if {[flagExists $args -mini]} {
	    
	    runMinimizer \
		-numSteps [expr 10 * $stepsBetweenChecks] \
		-stepsize 0.01 
	}
    
	#cleanCovlentGeom
   
	
	#
	# set up the vdws
	#
	
	marvinPyth command {import protocol}
	if {[flagExists $args -allVdw]} {
	    
	    marvinPyth command {protocol.initNBond(nbxmod=4,rcon=1)}	    
	    
	} else {
	    
	    marvinPyth command \
                 {protocol.initNBond(nbxmod=4,onlyCA=True,rcon=1,repel=1.2)}	    
	    
	}
	
	#    if {[catch {potList byName [$vdwPot instanceName]}]} {
	potList add [$vdwPot cget -this]
	#    }
	
		    
	#
	# find the Rgyr pot, if it exists in the potsToUse
	#
	
	set rgPot ""
	
	foreach elem $potsToUse {
	    if {[$elem instanceName] == "COLL"} {
		set rgPot $elem
	    }
	}
	
	
	for {set convergeCount 1} \
	    {$convergeCount < $maxConvergeTries} \
	    {incr convergeCount} {
	    
	    runMinimizer \
		-numSteps $stepsBetweenChecks \
		-stepsize 0.01 
	    
	    set vdwEner [$vdwPot calcEnergy]
	    set vdwPerRes [expr $vdwEner / double($nResidues)]	    
	    
	    puts -nonewline "convergence test:  "
	    
	    if {($vdwPerRes > $convergedVdwPerResidue)} {
		puts -nonewline "Ca..Ca VDW too large ($vdwPerRes). "
	    } else {
		puts -nonewline "Ca..Ca VDW ok. "
	    }
	    
	    if {$rgPot != ""} {
		set rgEner  [$rgPot calcEnergy]
		
		
		if {($rgEner > $convergedRgEner)} {
		    puts "Rg energy too large ($rgEner)."
		} else {
		    puts "Rg energy ok."
		}
	    } else {
		set rgEner 0
		puts "No Rg test."
	    }
	    
	    if {($vdwPerRes < $convergedVdwPerResidue) &&
		($rgEner < $convergedRgEner)} {
		
		set successful 1
		break
	    }
	}

	if {! $successful} {
	    puts "Minimizer couldn't converge.  Trying again with new random torsion angles."
	}
	
	

    }
    
    #
    # another mini without vdw to ensure good starting dihedral restraints, etc
    #
    
    potList removeAll
    foreach pot $potsToUse {
	#	    if {[catch {potList byName [$pot instanceName]}]} {
	potList add [$pot cget -this]
	#	    }
    }
    
    runMinimizer \
	-numSteps [expr 10 * $stepsBetweenChecks] \
	-stepsize 0.01 
    
    #
    # clean up  -- Note that there's something wrong with cleanup of smart pointers, 
    # since it leaves the command object for the underlying pointer lying around.  
    #

    potList removeAll


    if {[info exists bondPot]} {
	[$bondPot __deref__] -delete
	$bondPot -delete
    }
    if {[info exists anglePot]} {
	[$anglePot __deref__] -delete
	$anglePot -delete
    }
    if {[info exists improperPot]} {
	[$improperPot __deref__] -delete
	$improperPot -delete
    }
    if {[info exists vdwPot]} {
	[$vdwPot __deref__] -delete
	$vdwPot -delete
    }
}



#
# A standard high-temperature phase
#

proc highTemp args {

    set radius                        [flagVal $args -radius 1.2]
    set kvdw                          [flagVal $args -kvdw   1.0]
    set knoe                          [flagVal $args -knoe   30.0]
    set kdih                          [flagVal $args -kdih   2.0]
    set krama                         [flagVal $args -krama  0.02]
    set kinvnoe                       [flagVal $args -kInverseNoe 1.0]
    set characteristicDistViol        [requiredFlagVal $args -characteristicDistViol]
    set characteristicNoeCompleteness [requiredFlagVal $args -characteristicNoeCompleteness]
    set characteristicScatter         [requiredFlagVal $args -characteristicScatter]
    set DVweight                      [requiredFlagVal $args -distanceViolationWeight]
    set PLweight                      [requiredFlagVal $args -previousLikelihoodWeight]
    set NCweight                      [requiredFlagVal $args -noeCompletenessWeight]
    set PSweight                      [requiredFlagVal $args -scatterWeight]
    set characteristicDeltaDV         [flagVal $args -characteristicDeltaDV   0.1]
    set characteristicDeltaPL         [flagVal $args -characteristicDeltaPL   0.1]
    set characteristicDeltaNC         [flagVal $args -characteristicDeltaNC   0.1]
    set characteristicDeltaPS         [flagVal $args -characteristicDeltaPS   0.001]
    set bathTemp                      [flagVal $args -bathTemp 4000.01]
    set timeLength                    [flagVal $args -timeLength 50.0]
    set numNOEreshuffles              [flagVal $args -numNOEreshuffles 4]
    set potsToUse                     [flagVal $args -potsToUse {}]
    set noePotList                    [flagVal $args -noePot]
    set completeNoePotList            [flagVal $args -completeNoePot]
    set numMCsteps                    [flagVal $args -numMCsteps 1]
    set CDF                           [flagVal $args -characteristicDeltaFactor 1.0]
    set aft                           [flagVal $args -activationFilenameTemplate]
    set eToleranceFactor              [flagVal $args -eToleranceFactor 0.001]
    
    #
    # various vdW possibilities
    #
    
    if {[flagExists $args -CaVDW]} {

	XplorCommand \
	    "constraints 
                interaction (name ca or name p) (name ca or name p) weights * 1.0 end
                interaction (not (name ca or name p)) (all) weights * 1.0 vdw 0.0 end
             end"

	XplorCommand \
	    "param nbonds atom 
                nbxmod 4 
                wmin 0.01 
                cutnb 100.0
                tolerance 45.0 
                rexp 2 
                irex 2 
             end end"

	set vdwPot  [create_XplorPot VDW]
	lappend potsToUse $vdwPot

    } elseif {[flagExists $args -noVDW]} {

	XplorCommand \
	    "constraints 
                interaction (all) (all) weights * 1.0 vdw 0.0 end
             end"
	
    } else {
	
	XplorCommand \
	    "constraints 
                interaction (all) (all) weights * 1.0 end
             end"
	
	XplorCommand \
	    "param nbonds atom 
                nbxmod 4 
                wmin 0.01 
                cutnb 6.5
                tolerance 2.0 
                rexp 2 
                irex 2 
             end end"

	set vdwPot  [create_XplorPot VDW]
	lappend potsToUse $vdwPot
    }
    
    #
    # add covalent geometry potentials if they're not already there
    #
    
    set hasBond 0
    set hasAngle 0
    set hasImproper 0
	
    foreach elem $potsToUse {
	if {[$elem instanceName] == "BOND"} {
	    set hasBond 1
	}
	if {[$elem instanceName] == "ANGL"} {
	    set hasAngle 1
	}
	if {[$elem instanceName] == "IMPR"} {
	    set hasImproper 1
	}
    }

    if {! $hasBond} {
	set bondPot [create_XplorPot "BOND"]
	lappend potsToUse $bondPot
    }

    if {! $hasAngle} {
	set anglePot [create_XplorPot "ANGL"]
	lappend potsToUse $anglePot
    }

    if {! $hasImproper} {
	set improperPot [create_XplorPot "IMPR"]
	lappend potsToUse $improperPot
    }

    #
    # add active potentials to the potList
    #
    
    potList removeAll
    foreach pot $potsToUse {
	#	if {[catch {potList byName [$pot instanceName]}]} {
	potList add [$pot cget -this]
	#	}
    }

    #
    # set force constants, etc
    #

    # FIX: need to automatically calc cutnb here

    XplorCommand "param nbonds repel $radius rcon $kvdw end end"
    XplorCommand "rama scale $krama end"
    XplorCommand "restraints dihed scale $kdih end"
    
    foreach noePot $noePotList {
	
	$noePot setForceConst                          $knoe
	$noePot setInverseForceConst                   0
	$noePot setCharacteristicDistanceViol          $characteristicDistViol 
	$noePot setCharacteristicNoeCompleteness       0
	$noePot setCharacteristicScatter               9999.9
	$noePot setDistanceViolationWeight             $DVweight
	$noePot setPreviousLikelihoodWeight            $PLweight             
	$noePot setNoeCompletenessWeight               $NCweight             
	$noePot setScatterWeight                       $PSweight
	$noePot setCharacteristicDeltaDV               $characteristicDeltaDV  
	$noePot setCharacteristicDeltaPL               $characteristicDeltaPL  
	$noePot setCharacteristicDeltaNoeCompleteness  $characteristicDeltaNC
	$noePot setCharacteristicDeltaScatter          $characteristicDeltaPS
    }

    foreach completeNoePot $completeNoePotList {
	
	$completeNoePot setInverseForceConst                   $kinvnoe
	$completeNoePot setCharacteristicNoeCompleteness       $characteristicNoeCompleteness
	$completeNoePot setCharacteristicScatter               $characteristicScatter

    }

    if {$numNOEreshuffles < 1} {
	set numNOEreshuffles 1
    }

    set timePerCycle [expr $timeLength / double($numNOEreshuffles)]
    
    for {set curCycle 0} {$curCycle < $numNOEreshuffles} {incr curCycle} {
	
	puts "high temp cycle $curCycle of $numNOEreshuffles"

	foreach noePot $noePotList {
	    
	    updateActivationByMonteCarlo \
		-numMCsteps $numMCsteps \
		-characteristicDeltaFactor $CDF \
		-pot $noePot \
		-verbose
	    
	    if {$aft != ""} {
		writeActivePeakAssigns -pot $noePot -fileName [format $aft $curCycle]
	    }
	}

	runDynamics \
	    -timeLength $timePerCycle \
	    -bathTemp $bathTemp \
	    -energyToleranceFactor $eToleranceFactor
    }
    
    potList removeAll

    if {[info exists bondPot]} {
	[$bondPot __deref__] -delete
	$bondPot -delete
    }
    if {[info exists anglePot]} {
	[$anglePot __deref__] -delete
	$anglePot -delete
    }
    if {[info exists improperPot]} {
	[$improperPot __deref__] -delete
	$improperPot -delete
    }
    if {[info exists vdwPot]} {
	[$vdwPot __deref__] -delete
	$vdwPot -delete
    }

}


#
# A standard cooling phase
#


proc cooling args {

    #
    # grab the flags 
    # those with two vals specify coolstart and coolend points.
    #
    
    set radius                      [flagVal $args -radius {0.9 0.9}]
    set kvdw                        [flagVal $args -kvdw   {0.04 4.0}]
    set knoe                        [flagVal $args -knoe   {2.0 30.0}]
    set kdih                        [flagVal $args -kdih   {2.0 200.0}]
    set krama                       [flagVal $args -krama  {0.02 10.0}]
    set kinvnoe                     [flagVal $args -kInverseNoe {1.0 1.0}]
    set characteristicDistViol      [requiredFlagVal $args -characteristicDistViol]
    set characteristicNoeCompl      [requiredFlagVal $args -characteristicNoeCompleteness]
    set characteristicScatter       [requiredFlagVal $args -characteristicScatter]
    set DVweight                    [requiredFlagVal $args -distanceViolationWeight]
    set PLweight                    [requiredFlagVal $args -previousLikelihoodWeight]
    set NCweight                    [requiredFlagVal $args -noeCompletenessWeight]
    set PSweight                    [requiredFlagVal $args -scatterWeight]
    set characteristicDeltaDV       [flagVal $args -characteristicDeltaDV   {0.1 0.1}]
    set characteristicDeltaPL       [flagVal $args -characteristicDeltaPL   {0.1 0.1}]
    set characteristicDeltaNC       [flagVal $args -characteristicDeltaNC   {0.1 0.1}]
    set characteristicDeltaPS       [flagVal $args -characteristicDeltaPS   {0.001 0.001}]
    set bathTemp                    [flagVal $args -bathTemp {4000.01 100.0}]
    set timeLength                  [flagVal $args -timeLength 250.0]
    set numNOEreshuffles            [flagVal $args -numNOEreshuffles 32]
    set tempstep                    [flagVal $args -tempstep 25]
    set potsToUse                   [flagVal $args -potsToUse {}]
    set noePotList                  [flagVal $args -noePot]
    set completeNoePotList          [flagVal $args -completeNoePot]
    set numMCsteps                  [flagVal $args -numMCsteps 1]
    set CDF                         [flagVal $args -characteristicDeltaFactor 1.0]

    #
    # vdW setup
    #

    marvinPyth command {protocol.initNBond(nbxmod=4,repel=1)}	    
    
    #
    # add covalent geometry potentials if they're not already there
    #
    
    set hasBond 0
    set hasAngle 0
    set hasImproper 0
    set hasVdw 0
	
    foreach elem $potsToUse {
	if {[$elem instanceName] == "BOND"} {
	    set hasBond 1
	}
	if {[$elem instanceName] == "ANGL"} {
	    set hasAngle 1
	}
	if {[$elem instanceName] == "IMPR"} {
	    set hasImproper 1
	}
	if {[$elem instanceName] == "VDW"} {
	    set hasVdw 1
	}
    }

    if {! $hasBond} {
	set bondPot [create_XplorPot "BOND"]
	lappend potsToUse $bondPot
    }

    if {! $hasAngle} {
	set anglePot [create_XplorPot "ANGL"]
	lappend potsToUse $anglePot
    }

    if {! $hasImproper} {
	set improperPot [create_XplorPot "IMPR"]
	lappend potsToUse $improperPot
    }

    if {! $hasVdw} {
	set vdwPot [create_XplorPot "VDW"]
	lappend potsToUse $vdwPot
    }
    
    #
    # add active potentials to the potList
    #

    potList removeAll
    foreach pot $potsToUse {
	#	if {[catch {potList byName [$pot instanceName]}]} {
	potList add [$pot cget -this]
	#	}
    }
    
    set numCoolCycles [expr round(([firstItem $bathTemp] - [lastItem $bathTemp]) / double($tempstep))]
    set reshuffleProb [expr $numNOEreshuffles / double($numCoolCycles)]
    set timePerCycle  [expr $timeLength / double($numCoolCycles)]

    set bath [firstItem $bathTemp]
    for {set curCycle 0} {$curCycle < $numCoolCycles} {incr curCycle} {
	
	puts [format "cooling cycle %d of %d" $curCycle $numCoolCycles]
	
	set bath [expr $bath - $tempstep]

	set rad  [powerScale $radius  $numCoolCycles $curCycle]
	set kv   [powerScale $kvdw    $numCoolCycles $curCycle]
	set kr   [powerScale $krama   $numCoolCycles $curCycle]
	set kd   [powerScale $kdih    $numCoolCycles $curCycle]

	XplorCommand "param nbonds repel $rad rcon $kv end end"
	XplorCommand "rama scale $kr end"
	XplorCommand "restraints dihed scale $kd end"
	
	foreach noePot $noePotList {
	    $noePot setForceConst        [powerScale $knoe $numCoolCycles $curCycle]
	    $noePot setInverseForceConst 0
	}

	foreach completeNoePot $completeNoePotList {
	    $completeNoePot setInverseForceConst [powerScale $kinvnoe $numCoolCycles $curCycle]
	}

	if {[uniformRandom] <  $reshuffleProb} {

	    foreach completeNoePot $completeNoePotList {
		
		$completeNoePot  setCharacteristicNoeCompleteness       [linearScale $characteristicNoeCompl      $numCoolCycles $curCycle]
		$completeNoePot  setCharacteristicScatter               [linearScale $characteristicScatter       $numCoolCycles $curCycle]
	    }

	    foreach noePot $noePotList {
		
		$noePot setCharacteristicDistanceViol          [linearScale $characteristicDistViol      $numCoolCycles $curCycle]
		$noePot setDistanceViolationWeight             [linearScale $DVweight                    $numCoolCycles $curCycle]
		$noePot setPreviousLikelihoodWeight            [linearScale $PLweight                    $numCoolCycles $curCycle]
		$noePot setNoeCompletenessWeight               [linearScale $NCweight                    $numCoolCycles $curCycle]
		$noePot setScatterWeight                       [linearScale $PSweight                    $numCoolCycles $curCycle]
		$noePot setCharacteristicDeltaDV               [linearScale $characteristicDeltaDV       $numCoolCycles $curCycle]
		$noePot setCharacteristicDeltaPL               [linearScale $characteristicDeltaPL       $numCoolCycles $curCycle]
		$noePot setCharacteristicDeltaNoeCompleteness  [linearScale $characteristicDeltaNC       $numCoolCycles $curCycle]
		$noePot setCharacteristicDeltaScatter          [linearScale $characteristicDeltaPS       $numCoolCycles $curCycle]

		
		updateActivationByMonteCarlo \
		    -numMCsteps $numMCsteps \
		    -characteristicDeltaFactor $CDF \
		    -pot $noePot \
		    -verbose
	    }
	    
	}

	runDynamics \
	    -timeLength $timePerCycle \
	    -bathTemp $bath 
    }
    
    potList removeAll

    if {[info exists bondPot]} {
	[$bondPot __deref__] -delete
	$bondPot -delete
    }
    if {[info exists anglePot]} {
	[$anglePot __deref__] -delete
	$anglePot -delete
    }
    if {[info exists improperPot]} {
	[$improperPot __deref__] -delete
	$improperPot -delete
    }
    if {[info exists vdwPot]} {
	[$vdwPot __deref__] -delete
	$vdwPot -delete
    }
}


proc finalMinimization args {

    set potsToUse [flagVal $args -potsToUse {}]
    set numSteps  [flagVal $args -numSteps 200]
    set radius    [flagVal $args -radius 0.9]
    set kvdw      [flagVal $args -kvdw 4.0]

    #
    # add covalent geometry potentials if they're not already there
    #
    
    set hasBond 0
    set hasAngle 0
    set hasImproper 0
    set hasVdw 0
	
    foreach elem $potsToUse {
	if {[$elem instanceName] == "BOND"} {
	    set hasBond 1
	}
	if {[$elem instanceName] == "ANGL"} {
	    set hasAngle 1
	}
	if {[$elem instanceName] == "IMPR"} {
	    set hasImproper 1
	}
	if {[$elem instanceName] == "VDW"} {
	    set hasVdw 1
	}
    }

    if {! $hasBond} {
	set bondPot [create_XplorPot "BOND"]
	lappend potsToUse $bondPot
    }

    if {! $hasAngle} {
	set anglePot [create_XplorPot "ANGL"]
	lappend potsToUse $anglePot
    }

    if {! $hasImproper} {
	set improperPot [create_XplorPot "IMPR"]
	lappend potsToUse $improperPot
    }
    

    if {(! [flagExists $args -noVdw])  && (! $hasVdw)} {

	marvinPyth command \
             {protocol.initNBond(nbxmod=4,
                                 repel=float(repel),rcon=float(rcon))} \
	    [list [list repel $radius] [list rcon $kvdw]]
	

	set vdwPot  [create_XplorPot VDW]
	lappend potsToUse $vdwPot
    }


    potList removeAll
    foreach pot $potsToUse {
#	if {[catch {potList byName [$pot instanceName]}]} {
	    potList add [$pot cget -this]
#	}
    }

    runMinimizer \
	-numSteps $numSteps \
	-stepsize 0.01 

    potList removeAll
    
    if {[info exists bondPot]} {
	[$bondPot __deref__] -delete
	$bondPot -delete
    }
    if {[info exists anglePot]} {
	[$anglePot __deref__] -delete
	$anglePot -delete
    }
    if {[info exists improperPot]} {
	[$improperPot __deref__] -delete
	$improperPot -delete
    }
    if {[info exists vdwPot]} {
	[$vdwPot __deref__] -delete
	$vdwPot -delete
    }
}


proc evalCovalentGeom {} {
    
    XplorCommand "print threshold 0.01 bonds end"
    set curBondViols [XplorVariableNamed violations]
    
    XplorCommand "print threshold 5.0 angles end"
    set curAngleViols [XplorVariableNamed violations]
    
    XplorCommand "print threshold 5.0 impropers end"
    set curImprViols [XplorVariableNamed violations]

    return [list $curBondViols $curAngleViols $curImprViols]
}

    
proc cleanCovalentGeom args {

    set maxIters  [flagVal $args -maxIters 40] 

    XplorCommand "set message off echo off print off end"

    marvinPyth command {import protocol}
    if {[flagExists $args -useVdw]} {
	set useVdw 1
	marvinPyth command \
	    {protocol.initNBond(nbxmod=2,onlyCA=True,rcon=0.01,repel=0.9)}	    
	
      } else {
	set useVdw 0
    }

    if {[flagExists $args -useDynamics]} {
	set useMD 1
	XplorCommand "vector do (mass = 100) (all)"
	XplorCommand "vector do (fbeta = 10) (all)"
    } else {
	set useMD 0
    }

    set curViols [evalCovalentGeom]
    set curTotViols [expr [lindex $curViols 0] + [lindex $curViols 1] + [lindex $curViols 2]]

    set nIters 0
    while {($nIters < $maxIters) && ($curTotViols > 0)} {

	incr nIters

	#
	# randomly choose which minimization scheme to apply
	#

	if {$useVdw} {
	    set r [randomIntegerBetween 0 9]
	} else {
	    set r [randomIntegerBetween 0 7]
	}

	# FIX:  change to xplorPots
	
	switch $r {
	    
	    0 -
	    1 {XplorCommand "flags exclude * include bond end"}
	    2 -
	    3 {XplorCommand "flags exclude * include bond angle end"}
	    4 -
	    5 {XplorCommand "flags exclude * include bond impr end"}
	    6 -
	    7 {XplorCommand "flags exclude * include bond angle impr end"}
	    8 -
	    9 {XplorCommand "flags exclude * include bond angle impr vdw end"}
	}

	if {$useMD} {
	    set r [randomIntegerBetween 0 2]
	} else {
	    set r 0
	}

	switch $r {
	    
	    0 -
	    1 {XplorCommand "mini powell nstep 500 nprint 100 drop 10 end"}
	    2 {XplorCommand "dynamics verlet 
                               nstep 500 time 0.001 nprint 100 iprfrq 0
                               iasvel maxwell firsttemp 1000
                               tbath 1000
                            end"}
	}

	set curViols [evalCovalentGeom]
	set curTotViols [expr [lindex $curViols 0] + [lindex $curViols 1] + [lindex $curViols 2]]
    }

    if {$curTotViols > 0} {

	if {! $useVdw} {
	    set msg [format "%s\n(%d bond viols, %d angle viols, %d improper viols)\n%s" \
			 "Covalent geometry still violated after cleanCovalentGeom" \
			 [lindex $curViols 0] [lindex $curViols 1] [lindex $curViols 2] \
			 "Try increasing -maxIters flag or including -useVdw flag" ]
	} else {
	    set msg [format "%s\n(%d bond viols, %d angle viols, %d improper viols)\n%s" \
			 "Covalent geometry still violated after cleanCovalentGeom -useVdw" \
			 [lindex $curViols 0] [lindex $curViols 1] [lindex $curViols 2] \
			 "Try increasing -maxIters flag" ]
	}

	error $msg
    }

    XplorCommand "set message on echo on print OUTPUT end"
}


#
# equivalents of the original Marvin protocols
#

proc origPass1 args {

    global errorInfo

    #
    # grab the flags
    #

    set totStructs            [flagVal $args -numStructs 500]
    set outFilenameTemplate   [requiredFlagVal $args -outFilenameTemplate]
    set startingRemList       [flagVal $args -remarksList [list]]
    set potsToUseDuringSA     [requiredFlagVal $args -potsToUseDuringAnnealing]
    set potsToUseDuringRandom [requiredFlagVal $args -potsToUseDuringRandomization]
    set noePotList            [flagVal $args -noePot]
    set overallSeed           [flagVal $args -randomSeed 0]
    set startWithCurCoords    [flagExists $args -startWithCurrentCoords]
    
    #
    # save the starting coords
    #

    set origCoords [xplorSim atomPosArr]

    #
    # initial NOE behavior flags 
    #

    foreach noePot $noePotList {

	$noePot useSplitAssignmentBehavior
	$noePot useLinearPotential
	$noePot setViolationCutoff         1.0
	$noePot setInverseBound            4.0
	$noePot setInverseMethylCorrection 0.0
	$noePot setSwitchViolation         0.0
	$noePot noInversePotential

	marvinPyth command {import pasd}
	$noePot setLongRangePrimarySeqCutoff \
	    [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.longRangeResidCutoff} \
			   {} {out}] \
		      0] 1]
	$noePot updatePrimarySeqDists
	$noePot setMaximumMonteCarloAttempts 100
	$noePot useOriginalViolationScore
    }


    for {set count [firstStructNum $totStructs]} {$count < [lastStructNum  $totStructs]} {incr count} {
	
	puts "beginning of calculation of structure number $count"
	
	if {[catch {
	    
	    setMarvinRandomSeed -overall $overallSeed -pass 1 -count $count
	    
	    resetCoordsAndTemperature \
		-temperature 4000.0 \
		-coords $origCoords
	    
	    if {! $startWithCurCoords} {
		randomCollapse \
		    -bathTemp 4000 \
		    -potsToUse $potsToUseDuringRandom \
		    -mini \
		    -CaVDW \
		    -timeLength 10.0

		resetTemperature 4000.0
	    }
	    
	    #
	    # First high-temp stage:  optimize all possible peak & shift assigns -- nothing turned off
	    # 
	    # Activate all peakAssigns & shiftAssigns
	    #
	    
	    foreach noePot $noePotList {
		$noePot activateAllAssigns
	    } 
	    
	    # returned time here to 20ns
	    
	    XplorCommand "coor rgyr end"
	    puts [format "Rgyr = %s" [XplorVariableNamed rg]]

	    highTemp \
		-bathTemp 4000.0 \
		-timeLength 40.0 \
		-krama 0.1 \
		-knoe 1.0 \
		-kInverseNoe 0.0 \
		-kdih 200.0 \
		-eToleranceFactor 0.01 \
		-characteristicDistViol 9999999.9 \
		-characteristicNoeCompleteness 0.0 \
		-characteristicScatter 9999999.9 \
		-distanceViolationWeight  1.0 \
		-noeCompletenessWeight    1.0 \
		-scatterWeight            0.0 \
		-previousLikelihoodWeight 0.0 \
		-numNOEreshuffles 1 \
		-CaVDW \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList

	    XplorCommand "coor rgyr end"
	    puts [format "Rgyr = %s" [XplorVariableNamed rg]]


       	    #
	    # Second high-temp stage:  start cutting out assigns
	    #

	    highTemp \
		-bathTemp 4000 \
		-timeLength 40.0 \
		-krama 0.1 \
		-knoe 1.0 \
		-kInverseNoe 0.0 \
		-kdih 200.0 \
		-characteristicDistViol   5.5 \
		-characteristicNoeCompleteness 0.0 \
		-characteristicScatter    99999.9 \
		-distanceViolationWeight  1.0 \
		-previousLikelihoodWeight 0.0 \
		-noeCompletenessWeight    1.0 \
		-scatterWeight            0.0 \
		-characteristicDeltaDV    100.0 \
		-characteristicDeltaNC    0.001 \
		-characteristicDeltaPS    0.01 \
		-numNOEreshuffles 10 \
		-CaVDW \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-numMCsteps 1 \
		-characteristicDeltaFactor 1.0

	    #
	    # cool it, adding full vdw and dihedrals
	    #

	    XplorCommand "coor rgyr end"
	    puts [format "Rgyr = %s" [XplorVariableNamed rg]]

	    cooling \
		-distanceViolationWeight       {1.0 1.0} \
		-previousLikelihoodWeight      {0.0 0.0} \
		-noeCompletenessWeight         {1.0 1.0} \
		-scatterWeight                 {0.0 0.0} \
		-characteristicDistViol        {5.5 2.0} \
		-characteristicNoeCompleteness {0.0 0.0} \
		-characteristicScatter         {9999.9 9999.9} \
		-characteristicDeltaDV         {0.1 0.001} \
		-characteristicDeltaNC         {0.001 0.001} \
		-characteristicDeltaPS         {0.01 0.01} \
		-krama {0.1 10.0} \
		-knoe {1.0 30.0} \
		-kInverseNoe {0.0 0.0} \
		-radius {0.9 0.9} \
		-kdih {200.0 200.0} \
		-timeLength 250.0 \
		-bathTemp {4000.01 100.0} \
		-numNOEreshuffles 64 \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-numMCsteps 1 \
		-characteristicDeltaFactor 1.0


	    finalMinimization \
		-potsToUse $potsToUseDuringSA 
	    
	    set remList $startingRemList
	    
	    foreach line [energyReport $potsToUseDuringSA] { 
		lappend remList $line
	    }
	    lappend remList " "

	    foreach line [individualStructNoeReports $noePotList] {
		lappend remList $line
	    }
	    lappend remList " "
	    
	    lappend remList "Done with reworked original Marvin pass1 protocol"
	    
	    set outFileName [format $outFilenameTemplate $count]
	    writePDB -fileName $outFileName -remarks $remList
	    
	} errmsg ] } {
	    puts stderr "Note:  Error in calculating structure number $count in pass1.  Skipped. ($errmsg)"
	    puts stderr $errorInfo
	}
    }	
}

proc origPass2 args {

    global errorInfo

    #
    # grab the flags
    #

    set totStructs            [flagVal $args -numStructs 500]
    set outFilenameTemplate   [requiredFlagVal $args -outFilenameTemplate]
    set startingRemList       [flagVal $args -remarksList [list]]
    set potsToUseDuringSA     [requiredFlagVal $args -potsToUseDuringAnnealing]
    set potsToUseDuringRandom [requiredFlagVal $args -potsToUseDuringRandomization]
    set noePotList            [flagVal $args -noePot]
    set overallSeed           [flagVal $args -randomSeed 0]
    set startWithCurCoords    [flagExists $args -startWithCurrentCoords]

    #
    # save the starting coords
    #

    set origCoords [xplorSim atomPosArr]

    #
    # overall behavior flags 
    #

    foreach noePot $noePotList {

	$noePot useSplitAssignmentBehavior
	$noePot useLinearPotential
	$noePot setViolationCutoff         1.0
	$noePot setInverseBound            4.0
	$noePot setInverseMethylCorrection 0.0
	$noePot setSwitchViolation         0.0
	$noePot noInversePotential

	marvinPyth command {import pasd}
	$noePot setLongRangePrimarySeqCutoff \
	    [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.longRangeResidCutoff} \
			   {} {out}] \
		      0] 1]
	$noePot updatePrimarySeqDists
	$noePot setMaximumMonteCarloAttempts 100
	$noePot useOriginalViolationScore
	$noePot disallowShiftAssignmentInactivation
    }


    for {set count [firstStructNum $totStructs]} {$count < [lastStructNum  $totStructs]} {incr count} {
    
	puts "beginning of calculation of structure number $count"

	if {[catch {
	
	    setMarvinRandomSeed -overall $overallSeed -pass 2 -count $count
	    
	    resetCoordsAndTemperature \
		-temperature 4000.0 \
		-coords $origCoords

	    if {! $startWithCurCoords} {
		randomCollapse \
		    -bathTemp 4000 \
		    -potsToUse $potsToUseDuringRandom \
		    -mini \
		    -CaVDW \
		    -timeLength 10.0

		resetTemperature 4000.0
	    }
	    
	    foreach noePot $noePotList {
		$noePot activateAllAssigns
	    } 

	    #
	    # start cutting out assigns, based solely on the likelihoods from the previous pass
	    # to get us in the right neighborhood
	    #
	    
	    highTemp \
		-bathTemp 4000.0 \
		-timeLength 15.0 \
		-krama 0.1 \
		-knoe 1.0 \
		-kInverseNoe 0.0 \
		-kdih 200.0 \
		-characteristicDistViol 10.0 \
		-characteristicNoeCompleteness 0.0 \
		-characteristicScatter 9999999.9 \
		-distanceViolationWeight  0.0 \
		-noeCompletenessWeight    0.0 \
		-scatterWeight            0.0 \
		-previousLikelihoodWeight 1.0 \
		-numNOEreshuffles 10 \
		-CaVDW \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList

	    #
	    # now switch to cutting assigns based equally on likelihood & violation
	    #

	    highTemp \
		-bathTemp 4000.0 \
		-timeLength 40.0 \
		-krama 0.1 \
		-knoe 1.0 \
		-kInverseNoe 0.0 \
		-kdih 200.0 \
		-characteristicDistViol 5.5 \
		-characteristicNoeCompleteness 0.0 \
		-characteristicScatter    99999.9 \
		-distanceViolationWeight  0.5 \
		-noeCompletenessWeight    0.0 \
		-scatterWeight            0.0 \
		-previousLikelihoodWeight 0.5 \
		-characteristicDeltaDV    0.1 \
		-numNOEreshuffles 10 \
		-CaVDW \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList

	    #
	    # gradually switch over to cutting out on the basis of distance violation
	    #

	    cooling \
		-distanceViolationWeight  {0.5  1.0} \
		-previousLikelihoodWeight {0.5  0.0} \
		-noeCompletenessWeight    {0.0  0.0} \
		-scatterWeight            {0.0 0.0} \
		-characteristicDistViol        {5.5  2.0} \
		-characteristicNoeCompleteness {0.0 0.0} \
		-characteristicScatter         {9999.9 9999.9} \
		-characteristicDeltaDV         {0.1 0.001} \
		-characteristicDeltaNC         {0.01 0.01} \
		-characteristicDeltaPS         {0.01 0.01} \
		-krama {0.1 10.0} \
		-knoe {1.0 30.0} \
		-kInverseNoe {0.0 0.0} \
		-kdih {200.0 200.0} \
		-timeLength 250.0 \
		-bathTemp {4000.01 100.0} \
		-numNOEreshuffles 64 \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-numMCsteps 1 \
		-characteristicDeltaFactor 1.0

	    	    
	    finalMinimization \
		-potsToUse $potsToUseDuringSA 

	    set remList $startingRemList
	    
	    foreach line [energyReport $potsToUseDuringSA] { 
		lappend remList $line
	    }
	    lappend remList " "

	    foreach line [individualStructNoeReports $noePotList] {
		lappend remList $line
	    }
	    lappend remList " "

	    lappend remList "Done with Marvin pass2 protocol version 2"
	    
	    set outFileName [format $outFilenameTemplate $count]
	    writePDB -fileName $outFileName -remarks $remList
	    
	} errmsg ] } {
	    puts stderr "Note:  Error in calculating structure number $count in pass2.  Skipped. ($errmsg)"
	    puts stderr $errorInfo
	}
    }
}



proc origPass3 args {

    global errorInfo

    #
    # grab the flags
    #

    set totStructs            [flagVal $args -numStructs 500]
    set outFilenameTemplate   [requiredFlagVal $args -outFilenameTemplate]
    set startingRemList       [flagVal $args -remarksList [list]]
    set potsToUseDuringSA     [requiredFlagVal $args -potsToUseDuringAnnealing]
    set potsToUseDuringRandom [requiredFlagVal $args -potsToUseDuringRandomization]
    set noePotList            [flagVal $args -noePot]
    set overallSeed           [flagVal $args -randomSeed 0]
    set startWithCurCoords    [flagExists $args -startWithCurrentCoords]
    set ktalos                [flagVal $args -ktalos [list 200.0 200.0 200.0]]
    set krama                 [flagVal $args -krama  [list 0.2 0.1 10.0]]

    #
    # save the starting coords
    #

    set origCoords [xplorSim atomPosArr]

    #
    # overall behavior flags 
    #

    foreach noePot $noePotList {

	$noePot useSingleAssignmentBehavior
	$noePot useQuadraticPotential
	$noePot setViolationCutoff         1.0
	$noePot setInverseBound            4.0
	$noePot setInverseMethylCorrection 0.0
	$noePot noInversePotential

	marvinPyth command {import pasd}
	$noePot setLongRangePrimarySeqCutoff \
	    [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.longRangeResidCutoff} \
			   {} {out}] \
		      0] 1]
	$noePot updatePrimarySeqDists
	$noePot setMaximumMonteCarloAttempts 100
	$noePot useOriginalViolationScore
	$noePot disallowShiftAssignmentInactivation
    }


    for {set count [firstStructNum $totStructs]} {$count < [lastStructNum  $totStructs]} {incr count} {
		
	puts "beginning of calculation of structure number $count"

	if {[catch {
	
	    setMarvinRandomSeed -overall $overallSeed -pass 3 -count $count
	    
	    resetCoordsAndTemperature \
		-temperature 4000.0 \
		-coords $origCoords
    
	    if {! $startWithCurCoords} {
		randomCollapse \
		    -bathTemp 4000 \
		    -potsToUse $potsToUseDuringRandom \
		    -mini \
		    -CaVDW \
		    -timeLength 10.0

		resetTemperature 4000.0
	    }

	    foreach noePot $noePotList {
		$noePot activateAllAssigns
	    } 

	    #
	    # cut out restraints & decide assigns, based solely on the likelihoods from the previous pass
	    #
	    
	    highTemp \
		-bathTemp 4000.0 \
		-timeLength 50.0 \
		-krama [lindex $krama 0] \
		-knoe 3.0 \
		-kInverseNoe 0.0 \
		-kdih [lindex $ktalos 0] \
		-characteristicDistViol    9999.0 \
		-characteristicNoeCompleteness 0.0 \
		-characteristicScatter    9999.9 \
		-distanceViolationWeight  0.0 \
		-noeCompletenessWeight    0.0 \
		-scatterWeight            0.0 \
		-previousLikelihoodWeight 1.0 \
		-numNOEreshuffles 10 \
		-CaVDW \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList

	    #
	    # gradually switch over to cutting out restraints & deciding assigns based on distance viols
	    #

	    cooling \
		-distanceViolationWeight       {0.5 1.0} \
		-previousLikelihoodWeight      {0.5 0.0} \
		-noeCompletenessWeight         {0.0 0.0} \
		-scatterWeight                 {0.0 0.0} \
		-characteristicDistViol        {2.0 0.7} \
		-characteristicNoeCompleteness {0.0 0.0} \
		-characteristicScatter         {0.02 0.02} \
		-characteristicDeltaDV         {0.33 0.0033} \
		-characteristicDeltaNC         {0.01 0.01} \
		-characteristicDeltaPS         {0.01 0.01} \
		-krama [list [lindex $krama 1] [lindex $krama 2]] \
		-knoe {3.0 30.0} \
		-kInverseNoe {0.0 0.0} \
		-kdih [list [lindex $ktalos 1] [lindex $ktalos 2]] \
		-timeLength 250.0 \
		-bathTemp {4000.01 100.0} \
		-numNOEreshuffles 64 \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-numMCsteps 10 \
		-characteristicDeltaFactor 0.25
	    	    
	    finalMinimization \
		-potsToUse $potsToUseDuringSA 

	    set remList $startingRemList
	    
	    foreach line [energyReport $potsToUseDuringSA] { 
		lappend remList $line
	    }
	    lappend remList " "

	    foreach line [individualStructNoeReports $noePotList] {
		lappend remList $line
	    }
	    lappend remList " "

	    lappend remList "Done with original Marvin pass3 protocol"
	    
	    set outFileName [format $outFilenameTemplate $count]
	
	    writePDB -fileName $outFileName -remarks $remList
	    
	} errmsg ] } {
	    puts stderr "Note:  Error in calculating structure number $count in pass3.  Skipped. ($errmsg)"
	    puts stderr $errorInfo
	}
    }
}






#
# the main MARVIN protocol procedures
#


proc pass1 args {

    global errorInfo

    #
    # grab the flags
    #

    set totStructs            [flagVal $args -numStructs 500]
    set outFilenameTemplate   [requiredFlagVal $args -outFilenameTemplate]
    set startingRemList       [flagVal $args -remarksList [list]]
    set potsToUseDuringSA     [requiredFlagVal $args -potsToUseDuringAnnealing]
    set potsToUseDuringRandom [requiredFlagVal $args -potsToUseDuringRandomization]
    set noePotList            [flagVal $args -noePot]
    set completeNoePotList    [flagVal $args -completeNoePot]
    set overallSeed           [flagVal $args -randomSeed 0]
    set startWithCurCoords    [flagExists $args -startWithCurrentCoords]
    
    #
    # save the starting coords
    #

    set origCoords [xplorSim atomPosArr]

    #
    # initial NOE behavior flags 
    #

    foreach noePot $noePotList {

	$noePot useSplitAssignmentBehavior
	$noePot useLinearPotential
	$noePot setViolationCutoff         1.0
	$noePot setInverseForceConst       0.0
	$noePot setSwitchViolation         1.0
	$noePot noInversePotential

	marvinPyth command {import pasd}
	$noePot setLongRangePrimarySeqCutoff \
	    [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.longRangeResidCutoff} \
			   {} {out}] \
		      0] 1]
	$noePot updatePrimarySeqDists
	$noePot setMaximumMonteCarloAttempts 100
    }

    foreach completeNoePot $completeNoePotList {

	$completeNoePot useInversePotential
	$completeNoePot setInverseBound            4.0
	$completeNoePot setInverseMethylCorrection 0.0
   }


    for {set count [firstStructNum $totStructs]} {$count < [lastStructNum  $totStructs]} {incr count} {
	
	puts "beginning of calculation of structure number $count"
	
	if {[catch {
	    
	    setMarvinRandomSeed -overall $overallSeed -pass 1 -count $count
	    
	    resetCoordsAndTemperature \
		-temperature 4000.0 \
		-coords $origCoords
	    
	    if {! $startWithCurCoords} {
		randomCollapse \
		    -bathTemp 4000 \
		    -potsToUse $potsToUseDuringRandom \
		    -mini \
		    -CaVDW \
		    -timeLength 10.0

		resetTemperature 4000.0
	    }
	    
	    #
	    # First high-temp stage:  optimize all possible peak & shift assigns -- nothing turned off
	    # 
	    # Activate all peakAssigns & shiftAssigns
	    #
	    
	    foreach noePot $noePotList {
		$noePot activateAllAssigns
	    } 
	    
	    # returned time here to 20ns
	    
	    XplorCommand "coor rgyr end"
	    puts [format "Rgyr = %s" [XplorVariableNamed rg]]

	    highTemp \
		-bathTemp 4000.0 \
		-timeLength 20.0 \
		-krama 0.2 \
		-knoe 1.0 \
 		-kInverseNoe 5.0 \
		-kdih 200.0 \
		-eToleranceFactor 0.01 \
		-characteristicDistViol 9999999.9 \
		-characteristicNoeCompleteness 0.0 \
		-characteristicScatter 9999999.9 \
		-distanceViolationWeight  1.0 \
		-noeCompletenessWeight    1.0 \
		-scatterWeight            1.0 \
		-previousLikelihoodWeight 0.0 \
		-numNOEreshuffles 1 \
		-noVDW \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-completeNoePot $completeNoePotList

	    XplorCommand "coor rgyr end"
	    puts [format "Rgyr = %s" [XplorVariableNamed rg]]


       	    #
	    # Second high-temp stage:  start cutting out assigns
	    #

	    highTemp \
		-bathTemp 4000 \
		-timeLength 60.0 \
		-krama 0.2 \
		-knoe 1.0 \
		-kInverseNoe 5.0 \
		-kdih 200.0 \
		-characteristicDistViol   10.0 \
		-characteristicNoeCompleteness 0.10 \
		-characteristicScatter    1.0 \
		-distanceViolationWeight  1.0 \
		-previousLikelihoodWeight 0.0 \
		-noeCompletenessWeight    1.0 \
		-scatterWeight            0.0 \
		-characteristicDeltaDV    0.1 \
		-characteristicDeltaNC    0.001 \
		-characteristicDeltaPS    0.01 \
		-numNOEreshuffles 10 \
		-noVDW \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-completeNoePot $completeNoePotList \
		-numMCsteps 10 \
		-characteristicDeltaFactor 0.25

	    #
	    # cool it, adding full vdw and dihedrals
	    #

	    XplorCommand "coor rgyr end"
	    puts [format "Rgyr = %s" [XplorVariableNamed rg]]

	    cooling \
		-distanceViolationWeight       {1.0 1.0} \
		-previousLikelihoodWeight      {0.0 0.0} \
		-noeCompletenessWeight         {1.0 1.0} \
		-scatterWeight                 {0.0 0.0} \
		-characteristicDistViol        {10.0 2.0} \
		-characteristicNoeCompleteness {0.1 0.0} \
		-characteristicScatter         {1.0 1.0} \
		-characteristicDeltaDV         {0.1 0.01} \
		-characteristicDeltaNC         {0.001 0.001} \
		-characteristicDeltaPS         {0.01 0.01} \
		-krama {0.2 10.0} \
		-knoe {1.0 30.0} \
		-kInverseNoe {5.0 0.0} \
		-radius {0.8 0.9} \
		-kdih {200.0 200.0} \
		-timeLength 250.0 \
		-bathTemp {4000.01 100.0} \
		-numNOEreshuffles 64 \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-completeNoePot $completeNoePotList \
		-numMCsteps 10 \
		-characteristicDeltaFactor 0.25


	    finalMinimization \
		-potsToUse $potsToUseDuringSA 
	    
	    set remList $startingRemList
	    
	    foreach line [energyReport $potsToUseDuringSA] { 
		lappend remList $line
	    }
	    lappend remList " "

	    foreach line [individualStructNoeReports $noePotList] {
		lappend remList $line
	    }
	    lappend remList " "
	    
	    lappend remList "Done with Marvin pass1 protocol version 4--shiftAssignment likelihoods from noeCompl"
	    
	    set outFileName [format $outFilenameTemplate $count]
	    writePDB -fileName $outFileName -remarks $remList
	    
	} errmsg ] } {
	    puts stderr "Note:  Error in calculating structure number $count in pass1.  Skipped. ($errmsg)"
	    puts stderr $errorInfo
	}
    }	
}


proc pass2 args {

    global errorInfo

    #
    # grab the flags
    #

    set totStructs            [flagVal $args -numStructs 500]
    set outFilenameTemplate   [requiredFlagVal $args -outFilenameTemplate]
    set startingRemList       [flagVal $args -remarksList [list]]
    set potsToUseDuringSA     [requiredFlagVal $args -potsToUseDuringAnnealing]
    set potsToUseDuringRandom [requiredFlagVal $args -potsToUseDuringRandomization]
    set noePotList            [flagVal $args -noePot]
    set completeNoePotList    [flagVal $args -completeNoePot]
    set overallSeed           [flagVal $args -randomSeed 0]
    set startWithCurCoords    [flagExists $args -startWithCurrentCoords]
    set ktalos                [flagVal $args -ktalos [list 200.0 200.0 200.0]]
    set krama                 [flagVal $args -krama  [list 0.2 0.1 10.0]]

    #
    # save the starting coords
    #

    set origCoords [xplorSim atomPosArr]

    #
    # initial behavior flags 
    #

    foreach noePot $noePotList {
	
	$noePot useSplitAssignmentBehavior
	$noePot useLinearPotential
	$noePot setViolationCutoff         1.0
	$noePot setInverseForceConst       0.0
	$noePot setSwitchViolation         1.0
	$noePot noInversePotential
	
	marvinPyth command {import pasd}
	$noePot setLongRangePrimarySeqCutoff \
	    [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.longRangeResidCutoff} \
			   {} {out}] \
		      0] 1]
	$noePot updatePrimarySeqDists
	$noePot setMaximumMonteCarloAttempts 100
	$noePot useOriginalViolationScore
	$noePot disallowShiftAssignmentInactivation
    }
    
    foreach completeNoePot $completeNoePotList {
	
	$completeNoePot useInversePotential
	$completeNoePot setInverseBound            4.0
	$completeNoePot setInverseMethylCorrection 0.0
    }
    
    
    for {set count [firstStructNum $totStructs]} {$count < [lastStructNum  $totStructs]} {incr count} {
    
	puts "beginning of calculation of structure number $count"

	if {[catch {
	
	    setMarvinRandomSeed -overall $overallSeed -pass 2 -count $count
	    
	    resetCoordsAndTemperature \
		-temperature 4000.0 \
		-coords $origCoords

	    if {! $startWithCurCoords} {
		randomCollapse \
		    -bathTemp 4000 \
		    -potsToUse $potsToUseDuringRandom \
		    -mini \
		    -CaVDW \
		    -timeLength 10.0

		resetTemperature 4000.0
	    }

	    
	    
	    #
	    # start cutting out assigns, based solely on the likelihoods from the previous pass
	    # to get us in the right neighborhood
	    #

	    highTemp \
		-bathTemp 4000.0 \
		-timeLength 20.0 \
		-krama [lindex $krama 0] \
		-knoe 1.0 \
		-kInverseNoe 5.0 \
		-kdih [lindex $ktalos 0] \
		-eToleranceFactor 0.01 \
		-characteristicDistViol 10.0 \
		-characteristicNoeCompleteness 0.0 \
		-characteristicScatter 9999999.9 \
		-distanceViolationWeight  0.0 \
		-noeCompletenessWeight    0.0 \
		-scatterWeight            0.0 \
		-previousLikelihoodWeight 1.0 \
		-numNOEreshuffles 10 \
		-noVDW \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-completeNoePot $completeNoePotList

	    #
	    # now switch to cutting assigns based equally on likelihood & violation
	    #

	    highTemp \
		-bathTemp 4000 \
		-timeLength 60.0 \
		-krama [lindex $krama 0] \
		-knoe 1.0 \
		-kInverseNoe 5.0 \
		-kdih [lindex $ktalos 0] \
		-characteristicDistViol   10.0 \
		-characteristicNoeCompleteness 0.05 \
		-characteristicScatter    1.0 \
		-distanceViolationWeight  1.0 \
		-previousLikelihoodWeight 1.0 \
		-noeCompletenessWeight    1.0 \
		-scatterWeight            0.0 \
		-characteristicDeltaDV    0.1 \
		-characteristicDeltaNC    0.001 \
		-characteristicDeltaPS    0.01 \
		-numNOEreshuffles 10 \
		-noVDW \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-completeNoePot $completeNoePotList \
		-numMCsteps 10 \
		-characteristicDeltaFactor 0.25

	    #
	    # gradually switch over to cutting out on the basis of distance violation
	    #

	    cooling \
		-distanceViolationWeight  {1.0 1.0} \
		-previousLikelihoodWeight {1.0 0.0} \
		-noeCompletenessWeight    {1.0 1.0} \
		-scatterWeight            {0.0 0.0} \
		-characteristicDistViol        {10.0 2.0} \
		-characteristicNoeCompleteness {0.05 0.05} \
		-characteristicScatter         {1.0 1.0} \
		-characteristicDeltaDV         {0.01 0.01} \
		-characteristicDeltaNC         {0.01 0.01} \
		-characteristicDeltaPS         {0.01 0.01} \
		-krama [list [lindex $krama 1] [lindex $krama 2]] \
		-knoe {1.0 30.0} \
		-kInverseNoe {5.0 5.0} \
		-kdih [list [lindex $ktalos 1] [lindex $ktalos 2]] \
		-timeLength 250.0 \
		-bathTemp {4000.01 100.0} \
		-numNOEreshuffles 64 \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-completeNoePot $completeNoePotList \
		-numMCsteps 10 \
		-characteristicDeltaFactor 0.25


	    finalMinimization \
		-potsToUse $potsToUseDuringSA 

	    set remList $startingRemList
	    
	    foreach line [energyReport $potsToUseDuringSA] { 
		lappend remList $line
	    }
	    lappend remList " "

	    foreach line [individualStructNoeReports $noePotList] {
		lappend remList $line
	    }
	    lappend remList " "

	    lappend remList "Done with Marvin pass2 protocol version 2"
	    
	    set outFileName [format $outFilenameTemplate $count]
	    writePDB -fileName $outFileName -remarks $remList

	} errmsg ] } {
	    puts stderr "Note:  Error in calculating structure number $count in pass2.  Skipped. ($errmsg)"
	    puts stderr $errorInfo
	}
    }
}


proc pass3 args {

    global errorInfo

    #
    # grab the flags
    #

    set totStructs            [flagVal $args -numStructs 500]
    set outFilenameTemplate   [requiredFlagVal $args -outFilenameTemplate]
    set startingRemList       [flagVal $args -remarksList [list]]
    set potsToUseDuringSA     [requiredFlagVal $args -potsToUseDuringAnnealing]
    set potsToUseDuringRandom [requiredFlagVal $args -potsToUseDuringRandomization]
    set noePotList            [flagVal $args -noePot]
    set completeNoePotList    [flagVal $args -completeNoePot]
    set overallSeed           [flagVal $args -randomSeed 0]
    set startWithCurCoords    [flagExists $args -startWithCurrentCoords]

    #
    # save the starting coords
    #

    set origCoords [xplorSim atomPosArr]

    #
    # initial behavior flags 
    #

    foreach noePot $noePotList {

	$noePot useSingleAssignmentBehavior
	$noePot useQuadraticPotential
	$noePot setViolationCutoff         1.0
	$noePot setInverseForceConst       0.0
	$noePot setSwitchViolation         1.0
	$noePot noInversePotential

	marvinPyth command {import pasd}
	$noePot setLongRangePrimarySeqCutoff \
	    [lindex \
		 [lindex \
		      [marvinPyth command {out=pasd.longRangeResidCutoff} \
			   {} {out}] \
		      0] 1]
	$noePot updatePrimarySeqDists
	$noePot setMaximumMonteCarloAttempts 100
    }

    foreach completeNoePot $completeNoePotList {

	$completeNoePot useInversePotential
	$completeNoePot setInverseBound            4.0
	$completeNoePot setInverseMethylCorrection 0.0
    }


    for {set count [firstStructNum $totStructs]} {$count < [lastStructNum  $totStructs]} {incr count} {
		
	puts "beginning of calculation of structure number $count"

	if {[catch {
	
	    setMarvinRandomSeed -overall $overallSeed -pass 3 -count $count
	    
	    resetCoordsAndTemperature \
		-temperature 4000.0 \
		-coords $origCoords
    
	    if {! $startWithCurCoords} {
		randomCollapse \
		    -bathTemp 4000 \
		    -potsToUse $potsToUseDuringRandom \
		    -mini \
		    -CaVDW \
		    -timeLength 10.0

		resetTemperature 4000.0
	    }
	   
	    #
	    # cut out restraints & decide assigns, based solely on the likelihoods from the previous pass
	    #
	    
	    highTemp \
		-bathTemp 4000.0 \
		-timeLength 50.0 \
		-krama 0.1 \
		-knoe 3.0 \
		-kInverseNoe 1.0 \
		-kdih 200.0 \
		-characteristicDistViol 2.0 \
		-characteristicNoeCompleteness 0.0 \
		-characteristicScatter    0.02 \
		-distanceViolationWeight  0.0 \
		-noeCompletenessWeight    0.0 \
		-scatterWeight            0.0 \
		-previousLikelihoodWeight 1.0 \
		-numNOEreshuffles 10 \
		-noVDW \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-completeNoePot $completeNoePotList

	    #
	    # gradually switch over to cutting out restraints & deciding assigns based on distance viols
	    #

	    cooling \
		-distanceViolationWeight       {0.5 1.0} \
		-previousLikelihoodWeight      {0.5 0.0} \
		-noeCompletenessWeight         {0.5 1.0} \
		-scatterWeight                 {0.0 0.0} \
		-characteristicDistViol        {2.0 0.7} \
		-characteristicNoeCompleteness {0.1 0.1} \
		-characteristicScatter         {0.02 0.02} \
		-characteristicDeltaDV         {0.33 0.0033} \
		-characteristicDeltaNC         {0.01 0.01} \
		-characteristicDeltaPS         {0.01 0.01} \
		-krama {0.1 10.0} \
		-knoe {3.0 30.0} \
		-kInverseNoe {1.0 1.0} \
		-kdih {200.0 200.0} \
		-timeLength 250.0 \
		-bathTemp {4000.01 100.0} \
		-numNOEreshuffles 64 \
		-potsToUse $potsToUseDuringSA \
		-noePot $noePotList \
		-completeNoePot $completeNoePotList \
		-numMCsteps 10 \
		-characteristicDeltaFactor 0.25
	    	    
	    finalMinimization \
		-potsToUse $potsToUseDuringSA 

	    set remList $startingRemList
	    
	    foreach line [energyReport $potsToUseDuringSA] { 
		lappend remList $line
	    }
	    lappend remList " "

	    foreach line [individualStructNoeReports $noePotList] {
		lappend remList $line
	    }
	    lappend remList " "

	    lappend remList "Done with Marvin pass3 protocol version 3 (for use from beginning, inverse forces)"
	    
	    set outFileName [format $outFilenameTemplate $count]
	
	    writePDB -fileName $outFileName -remarks $remList
	    
	} errmsg ] } {
	    puts stderr "Note:  Error in calculating structure number $count in pass3.  Skipped. ($errmsg)"
	    puts stderr $errorInfo
	}
    }
}

