#
# sa_tools.tcl
#
# TCL procedures of use during annealing
#
# JJK 5/9/02
#

package provide marvin 1.0
package require pyinterp
package require simulationworld

#
# create a global (to Marvin) python interpreter object
# Note that this is used by procs defined in other 
# TCL files as well!
#

if {[catch {PyInterp marvinPyth} e]} {
    puts stderr [format "Error creating Python interpreter at startup:\n%s" $e]
    puts stdout [format "Error creating Python interpreter at startup:\n%s" $e]
    exit 1
}


#
# scaling functions
#

#
# exponential scaling mechanism:
#
# curVal = initVal * ((finalVal / initVal)^(1/totCycles))^curCycle
#
# Can be called with initVal finalVal totCycles curCycle or 
# {initVal finalVal} totCycles curCycle
#
# Note that cycle numbers need to start at zero for the
# scaled value to start at initVal.  
#
# This means that there will be no cycle called with curCycle == totCycles,
# since N cycles will be called with curCycle = 0, 1, ...N-1
#
# Therefore, I correct for this by subtracting one from totCycles, 
# so that in the final cycle, where curCycle == N-1, the scaled 
# value will be finalVal
#
# Note that if finalVal == 0, the scaled value will be zero for any 
# curCycle > 0
#

proc powerScale args {

    if {[llength $args] == 4} {
	set initVal   [lindex $args 0]
	set finalVal  [lindex $args 1]
	set totCycles [lindex $args 2]
	set curCycle  [lindex $args 3]
    } elseif {[llength $args] == 3} {
	set temp [lindex $args 0]
	if {[llength $temp] > 1} {
	    set initVal   [lindex $temp 0]
	    set finalVal  [lindex $temp 1]
	} else {
	    set initVal $temp
	    set finalVal $temp
	}
	set totCycles [lindex $args 1]
	set curCycle  [lindex $args 2]
    } else {
	error "powerScale not called with 3 or 4 arguments"
    }
    
    set totCycles [expr $totCycles - 1]

    if {$initVal == 0} {
	return $initVal
    } elseif {$totCycles < 1} {
        return $finalVal
    } else {
        return [expr $initVal * (pow( (double($finalVal) / double($initVal)), \
					  (double($curCycle)/double($totCycles))))]
    }
}

#
# linear scaling mechanism:
#
# curVal = ((curCycle / totCycles) * (finalVal - initVal)) + initVal
#
# Can be called with initVal finalVal totCycles curCycle or 
# {initVal finalVal} totCycles curCycle
#
# Note that cycle numbers need to start at zero for the
# scaled value to start at initVal.  
#
# This means that there will be no cycle called with curCycle == totCycles,
# since N cycles will be called with curCycle = 0, 1, ...N-1
#
# Therefore, I correct for this by subtracting one from totCycles, 
# so that in the final cycle, where curCycle == N-1, the scaled 
# value will be finalVal
#

proc linearScale args {
    
    if {[llength $args] == 4} {
	set initVal   [lindex $args 0]
	set finalVal  [lindex $args 1]
	set totCycles [lindex $args 2]
	set curCycle  [lindex $args 3]
    } elseif {[llength $args] == 3} {
	set temp [lindex $args 0]
	if {[llength $temp] > 1} {
	    set initVal   [lindex $temp 0]
	    set finalVal  [lindex $temp 1]
	} else {
	    set initVal $temp
	    set finalVal $temp
	}
	set totCycles [lindex $args 1]
	set curCycle  [lindex $args 2]
    } else {
	error "linearScale not called with 3 or 4 arguments"
    }
    
    set totCycles [expr $totCycles - 1]

    if {($totCycles < 1)} {
        return $finalVal
    } else {
        return [expr ((double($curCycle)/double($totCycles))* \
			  (double($finalVal)-double($initVal)))+double($initVal)]
    }
}

proc runMinimizer args {

    #
    # grab the name of the Python-side IVM object 
    # from the MarvinDefaults, which is defined in
    # setupIVM
    #
    
    namespace eval MarvinDefaults {
	variable IVMObjectName
    }
    set ivmObjectName $MarvinDefaults::IVMObjectName
    if {$ivmObjectName == {}} {
	error "no IVM object name defined in RunMinimizer.  Did you run setupIVM?"
    } 

    
    #
    # normal flags here.  The Python IVM object
    # maintains default values from the last execution
    # 
    
    set temp [flagVal $args -numSteps]
    if {$temp != {}} {
	marvinPyth command "$ivmObjectName.setNumSteps($temp)"
    }
    
    set temp [flagVal $args -printInterval]
    if {$temp != {}} {
	marvinPyth command "$ivmObjectName.setPrintInterval($temp)"
    }
    
    set temp [flagVal $args -stepsize] 
    if {$temp != {}} {
	marvinPyth command "$ivmObjectName.setStepsize($temp)"
    }

    #
    # Energy tolerance must be zero for minimization
    #

    marvinPyth command "$ivmObjectName.setETolerance(0)"

    # 
    # Since this is minimization, force stepType to be powell
    #
    
    set cmd [format "%s.setStepType(\"powell\")" $ivmObjectName]
    marvinPyth command $cmd
    
    # 
    # Run the minimization, trapping any Python exceptions
    #

    set cmd [format "fail=0\ntry:\n   %s.run()\nexcept SystemError:\n   fail=1; import traceback; traceback.print_exc()  " \
		 $ivmObjectName]

    set outVars [marvinPyth command $cmd {} {fail}]

    #
    # check the TCL variables I extracted for one
    # named 'fail' with a value != 0.  Raise a TCL
    # exception if it's there.
    #
    
    foreach item $outVars {
	
	set name [lindex $item 0]
	set val  [lindex $item 1]
	
	if {($name == "fail") && ($val != "0")} {
	    error [format "Python failure: %s" $val]
	}
    }    
}

    

#
# runDynamics
#
# A TCL wrapper around the Python calls to run a 
# dynamics trajectory.  Raises a TCL exception 
# if it crashes
#

proc runDynamics args {

    #
    # grab the name of the Python-side IVM object 
    # from the MarvinDefaults, which is defined in
    # setupIVM.  Also grab the EToleranceFactor and 
    # MaxDeltaEFactor
    #
    
    namespace eval MarvinDefaults {
	variable IVMObjectName
	variable EToleranceFactor
	variable MaxDeltaEFactor
    }
    set ivmObjectName $MarvinDefaults::IVMObjectName
    if {$ivmObjectName == {}} {
	error "no IVM object name defined in RunMinimizer.  Did you run setupIVM?"
    } 

    
    #
    # normal flags here.  The Python IVM object
    # maintains default values from the last execution
    # 
    
       
    set temp [flagVal $args -timeLength]
    if {$temp != {}} {
	marvinPyth command "$ivmObjectName.setFinalTime($temp)"
    }
    
    set temp [flagVal $args -printInterval]
    if {$temp != {}} {
	marvinPyth command "$ivmObjectName.setPrintInterval($temp)"
    }
    
    set temp [flagVal $args -resetCMInterval]
    if {$temp != {}} {
	marvinPyth command "$ivmObjectName.setResetCMInterval($temp)"
    }
    
    set temp [flagVal $args -bathTemp]
    if {$temp != {}} {
	marvinPyth command "$ivmObjectName.setBathTemp($temp)"
    }
    
    set temp [flagVal $args -responseTime]
    if {$temp != {}} {
	marvinPyth command "$ivmObjectName.setResponseTime($temp)"
    }
    
    set temp [flagVal $args -stepsize] 
    if {$temp != {}} {
	marvinPyth command "$ivmObjectName.setStepsize($temp)"
    }

    #
    # energy tolerance stuff:
    #
    # Values for eTolerance and MaxDeltaE depend on both the number of
    # atoms and on the bath temperature.
    #

    set outVars [marvinPyth command "bt=$ivmObjectName.bathTemp()" {} {bt}]
    
    foreach item $outVars {
	
	set name [lindex $item 0]
	set val  [lindex $item 1]
	
	if {$name == "bt"} {
	    set bathTemp $val
	}
    }

    set et  [expr $MarvinDefaults::EToleranceFactor * [xplorSim numAtoms] * $bathTemp]
    set mde [expr $MarvinDefaults::MaxDeltaEFactor  * [xplorSim numAtoms] * $bathTemp]

    # overwrite my values with calc from GMC script

    set eToleranceFactor [flagVal $args -energyToleranceFactor 0.001]
    set et [expr $bathTemp * $eToleranceFactor]
    set mde 100000.0
    
    marvinPyth command "$ivmObjectName.setETolerance($et)"
    marvinPyth command "$ivmObjectName.setMaxDeltaE($mde)"

    # 
    # Since this is MD, force stepType to be pc6, with adjustable
    # stepsizes and no fixed number of them
    #
    
    set cmd [format "%s.setStepType(\"pc6\")" $ivmObjectName]
    marvinPyth command $cmd
    
    marvinPyth command "$ivmObjectName.setNumSteps(0)"

    if {[flagExists $args -constantStepsize]} {
	marvinPyth command "$ivmObjectName.setAdjustStepsize(0)"
    } else {
	marvinPyth command "$ivmObjectName.setAdjustStepsize(1)"
    }

    #
    # randomize the velocities around the bath temperature
    #

    marvinPyth command "randomizeVelocities($bathTemp)"

    # debugging
    # marvinPyth command "$ivmObjectName.setPrintInterval(1)"
    # marvinPyth command "$ivmObjectName.setNumSteps(269)"

    # 
    # Run the dynamics, trapping any Python exceptions
    #


    set cmd [format "fail=0\ntry:\n   %s.run()\nexcept SystemError:\n   fail=1; import traceback; traceback.print_exc()  " \
		 $ivmObjectName]

    set outVars [marvinPyth command $cmd {} {fail}]

    #
    # check the TCL variables I extracted for one
    # named 'fail' with a value != 0.  Raise a TCL
    # exception if it's there.
    #
    
    foreach item $outVars {
	
	set name [lindex $item 0]
	set val  [lindex $item 1]
	
	if {($name == "fail") && ($val != "0")} {
	    error [format "Python failure: %s" $val]
	}
    }  

    # debugging
    #   writePDB -fileName tmp.pdb
    #   exit
}


proc grabPDBFiles args {
    
    #
    # a simple routine to assemble lists of PDB files for structure-summary
    # tasks, with checking for average PDB files along the way
    #
    # Reads each one, and returns a list of <filename, atomPosArr pointer> pairs
    #

    set fnames     [requiredFlagVal $args -fileNames]
    set modelNums  [flagVal $args -modelNumbers]
    set avgPDBfile [flagVal $args -avgPDBfile]
    
    marvinPyth command "from selectTools import correctSymmetricSidechains"

    #
    # make sure fnames doesn't include the averagePDBfile
    #

    if {$avgPDBfile != ""} {
	set loc [lsearch -exact $fnames $avgPDBfile]
	if {$loc != -1} {
	    puts stderr "Removing file $avgPDBfile from analysis because it is also the name given for the avgPDBfile"
	    set fnames [lreplace $fnames $loc $loc]
	}
    }

    #
    # make sure fnames doesn't include anything that looks like
    # an averaged PDB file
    #

    if {! [flagExists $args -noAvgChecking]} {

	set done 0
	while {! $done} {
	    set loc [lsearch -glob $fnames *_avg.*]
	    if {$loc == -1} {
		set done 1
	    } else {
		puts stderr [format "Removing file %s because it looks like an averaged PDB file." \
				 [lindex $fnames $loc]]
		puts stderr "If this is in error, re-run grabPDBFiles with -noAvgChecking flag."
		set fnames [lreplace $fnames $loc $loc]
	    }
	}
    }

    #
    # make sure fnames doesn't include anything that looks like 
    # a converged link to a PDB file
    #

    if {! [flagExists $args -noConvergedChecking]} {


	set done 0
	while {! $done} {
	    set loc [lsearch -glob $fnames *_converged_*]
	    if {$loc == -1} {
		set done 1
	    } else {
		puts stderr [format "Removing file %s because it looks like a converged PDB file link." \
				 [lindex $fnames $loc]]
		puts stderr "If this is in error, re-run grabPDBFiles with -noConvergedChecking flag."
		set fnames [lreplace $fnames $loc $loc]
	    }
	}
    }


    #
    # read each PDB file, grab its coords from the simulation, 
    # and return them in a list with their filenames
    #

    set retVal [list]
    set fileCount 0

    foreach curName $fnames curModNum $modelNums {

	if {$curModNum == ""} {
	    updateUser [format "reading %s (%d of %d)\r" $curName [incr fileCount] [llength $fnames]]
	    readPDB -fileName $curName
	    marvinPyth command "correctSymmetricSidechains()"
	    set curAPA [xplorSim atomPosArr]
	    lappend retVal [list $curName $curAPA]
	} else {
	    updateUser [format "reading %s model %s (%d of %d)\r" $curName $curModNum [incr fileCount] [llength $fnames]]
	    readPDB -fileName $curName -model $curModNum
	    marvinPyth command "correctSymmetricSidechains()"
	    set curAPA [xplorSim atomPosArr]
	    lappend retVal [list [format "%s:%s" $curName $curModNum] $curAPA]
	}

    }

    return $retVal
}


proc meanStruct args {

    set inFileNames  [flagVal $args -fileNames]
    set files        [flagVal $args -files]
    set sel          [flagVal $args -selection [AtomSel -args "(name ca or name c or name n)"]]

    if {($inFileNames == "") && ($files == "")} {
	error "Neither -fileNames nor -files defined in call to meanStruct"
    }

    set selString [$sel string]

    if {[llength $inFileNames] > 0} {
	set files [grabPDBFiles -fileNames $inFileNames -noAvgChecking -noConvergedChecking]
    }

    set firstFile [lindex $files 0]

    xplorSim setAtomPosArr [lindex $firstFile 1]

    XplorCommand "coor copy end"
    XplorCommand "vector do (store1 = 0) (all)"
    XplorCommand "vector do (store2 = 0) (all)"
    XplorCommand "vector do (store3 = 0) (all)"

    foreach curFile $files {

	xplorSim setAtomPosArr [lindex $curFile 1]
	XplorCommand "coor fit sele ($selString) end"
	XplorCommand "vector do (store1 = store1 + x) (all)"
	XplorCommand "vector do (store2 = store2 + y) (all)"
	XplorCommand "vector do (store3 = store3 + z) (all)"
    }
    
    set nFiles [llength $files]

    XplorCommand "vector do (x = store1 / $nFiles) (all)"
    XplorCommand "vector do (y = store2 / $nFiles) (all)"
    XplorCommand "vector do (z = store3 / $nFiles) (all)"

    if {! [flagExists $args -noCleanCovalentGeom]} {
	if {[catch {cleanCovalentGeom -maxIters 100 -useVdw} msg]} {
	    puts [format "WARNING::meanStruct: %s" $msg]
	}
    }
}

proc compareToReference args {

    set refFileName [requiredFlagVal $args -referenceFileName]
    set sel         [flagVal $args -selection [AtomSel -args "(name ca or name c or name n)"]]

    set selString [$sel string]

    XplorCommand "coor copy end"
    readPDB -fileName $refFileName
    XplorCommand "coor fit sele ($selString) end"
    XplorCommand "coor rms sele ($selString) end"
    set rmsd [XplorVariableNamed result]
    XplorCommand "coor swap end"

    return $rmsd
}

proc structPrecision args {

    set inFileNames  [flagVal $args -fileNames]
    set files        [flagVal $args -files]
    set bbnSel       [flagVal $args -backboneSelection [AtomSel -args "(name ca or name c or name n)"]]
    set heavySel     [flagVal $args -heavyatomSelection [AtomSel -args "(not name h*)"]]

    XplorCommand "coor copy end"

    if {($inFileNames == "") && ($files == "")} {
	error "Neither -fileNames nor -files defined in call to structPrecision"
    }

    if {[llength $inFileNames] > 0} {
	set files [grabPDBFiles -fileNames $inFileNames -noAvgChecking -noConvergedChecking]
    }

    set bbnString   [$bbnSel string]
    set heavyString [$heavySel string]

    if {[flagExists $args -calcMean]} {
	meanStruct -files $files -selection $bbnSel
	XplorCommand "coor copy end"
    }

    set backboneRMSDs  [list]
    set heavyatomRMSDs [list]


    set retVal "Precision report:"

    foreach curFile $files {
	
	xplorSim setAtomPosArr [lindex $curFile 1]
	XplorCommand "coor fit sele ($bbnString) end"
	XplorCommand "coor rms sele ($bbnString) end"
	set curBackboneRMSD [XplorVariableNamed result]
	lappend backboneRMSDs $curBackboneRMSD

	XplorCommand "coor fit sele ($heavyString) end"
	XplorCommand "coor rms sele ($heavyString) end"
	set curHeavyatomRMSD [XplorVariableNamed result]
	lappend heavyatomRMSDs $curHeavyatomRMSD


	if {[flagExists $args -verbose]} {
	    set retVal [format "%s\n   %s %f %f" $retVal [lindex $curFile 0] $curBackboneRMSD $curHeavyatomRMSD]
	}
    }
    
    XplorCommand "coor swap end"

    set retVal [format "%s\n   mean backbone  RMSD to current struct is %f +/- %f" \
		    $retVal [mean $backboneRMSDs] [standardDeviation $backboneRMSDs]]
    
    set retVal [format "%s\n   mean heavyatom RMSD to current struct is %f +/- %f" \
		    $retVal [mean $heavyatomRMSDs] [standardDeviation $heavyatomRMSDs]]
    
    return $retVal
}


proc createViewScript args {

    set scriptFileName [flagVal $args -scriptFileName "view.inp"]
    set psfFileName    [requiredFlagVal $args -psfFileName]
    set inFileNames    [requiredFlagVal $args -pdbFileNames]
    set fitToFile      [flagVal $args -fitToFileName]
    set sel            [flagVal $args -selection [AtomSel -args "(name ca or name c or name n)"]]

    if {[catch {set outUnit [open $scriptFileName w]}]} {
	
	error "createViewScript::Error opening output file $scriptFileName"
    }

    puts $outUnit [format "struct @%s end\n" $psfFileName]

    if {$fitToFile != ""} {
	puts $outUnit [format "coor @%s" $fitToFile]
	puts $outUnit "coor copy end"
    }

    puts $outUnit [format "for \$infile in ("]

    foreach inFileName $inFileNames {
	puts $outUnit [format "   \"%s\" " $inFileName]
    }

    puts $outUnit "   ) loop f" 
    puts $outUnit [format "   coor @@\$infile"]
    
    if {$fitToFile != ""} {
	puts $outUnit [format "   coor fit sele %s end" [$sel string]]
    }

    puts $outUnit [format "   ctcl \"set xplor_result \[string trimleft \[XplorVariableNamed infile \] ./ \] \" "]
    puts $outUnit [format "   ps define \$result bonds %s color name blue end end end" [$sel string]]
    puts $outUnit "end loop f"

    close $outUnit
}
    

proc setupDelphicTorsions args {

    upvar env env

    # note that I need curly braces here for $krama

    XplorCommand {eval ($krama = 1.0)}

    XplorCommand \
	"set echo off message off end
         rama 
            nres 10000
            @GAUSSIANS:shortrange_gaussians.tbl
            @GAUSSIANS:new_shortrange_force.tbl
         end
         set echo on message on end"

    XplorCommand "@GAUSSIANS:newshortrange_setup.tbl"

}



proc setupIVM args {

    set fixedSelects  [flagVal $args -fixedSelections]
    set rigidSelects  [flagVal $args -rigidSelections]
    set tolMult       [flagVal $args -toleranceMultiplier 1.0]
   
    marvinPyth command "import ivm"
    marvinPyth command "import selectTools"
    marvinPyth command "from xplor import select"
    marvinPyth command "from selectTools import IVM_groupRigidSidechain"
    marvinPyth command "from selectTools import IVM_groupRigidBackbone"
    marvinPyth command "from selectTools import IVM_breakRiboses"
    marvinPyth command "from selectTools import IVM_breakProlines"
    marvinPyth command "from selectTools import IVM_breakDisulfides"
    marvinPyth command "from atomAction import *"
    marvinPyth command "import xplor"
    
    marvinPyth command "sim=xplor.simulation"
    marvinPyth command "d = ivm.IVM(sim)"
    marvinPyth command "d.setVerbose(0)"

    #
    # EToleranceFactor and MaxDeltaEFactor are used to set 
    # values of eTolerance and maxDeltaE.  Actual values of these params also 
    # depend on the number of atoms in the system and on
    # the bath temperature.  They're calculated for each call to
    # runDynamics
    #

    namespace eval MarvinDefaults {
	variable IVMObjectName
	variable EToleranceFactor
	variable MaxDeltaEFactor
    }

    set MarvinDefaults::IVMObjectName "d"

    #
    # raw values are designed to give CDS's default values for 
    # eTolerance and maxDeltaE in the protein G case @ 4000K
    #

    set MarvinDefaults::EToleranceFactor [expr $tolMult * 1.54e-6]
    set MarvinDefaults::MaxDeltaEFactor  [expr $tolMult * 3.85e-5]

    #
    # generate all the fixed and rigid groups
    #

    foreach curSel $fixedSelects {

	set tempString [$curSel string]
	marvinPyth command "d.fix( '$tempString' )"
    }

    foreach curSel $rigidSelects {

	set tempString [$curSel string]
	marvinPyth command "d.group( '$tempString' )"
    }


    marvinPyth command "IVM_groupRigidSidechain(d)"
    marvinPyth command "IVM_groupRigidBackbone(d)"
    marvinPyth command "IVM_breakProlines(d)"
    marvinPyth command "IVM_breakRiboses(d)"
    marvinPyth command "IVM_breakDisulfides(d)"

    marvinPyth command "d.autoTorsion()"

    marvinPyth command "d.setPrintInterval(10)"
    marvinPyth command "d.setResetCMInterval(10)"
    marvinPyth command "d.setResponseTime(5)"
    marvinPyth command "d.setAdjustStepsize(1)"

    #
    # grab the IVM's potList and create a TCL version of it
    #

    marvinPyth command {from pyInterp import portableStringRep}
    set tmp [marvinPyth command {ret=portableStringRep(d.potList())} {} \
	     {ret}]
    set potListPtr [lindex [lindex $tmp 0] 1]
    rcPotList potList -this $potListPtr
}


proc firstStructNum {tot} {

    upvar \#0 env env
    return [expr int(($env(XPLOR_PROCESS) * $tot) / $env(XPLOR_NUM_PROCESSES))]
}

proc lastStructNum {tot} { 

    upvar \#0 env env
    return [expr int((($env(XPLOR_PROCESS) + 1) * $tot) / $env(XPLOR_NUM_PROCESSES))]
}

proc initParamPsfPdb args {

    global env
    
    set paramFileNames [flagVal $args -paramFileName $env(TOPPAR)/protein.par]
    set psfFileNames   [requiredFlagVal $args -psfFileName]
    set pdbFileNames   [flagVal $args -pdbFileName ""]
    set makeRandom     [flagExists $args -randomCoords]
    
    foreach elem $paramFileNames {
	if {! [file exists $elem]} {
	    error "Parameter file $elem does not exist"
	}
	
	XplorCommand "param @$elem end"
    }

    foreach elem $psfFileNames {
	if {! [file exists $elem]} {
	    error "PSF file $elem does not exist"
	}
	
	XplorCommand "struct @$elem end"
    }

    if {$makeRandom} {
 	set tempSel [AtomSel -args "all"]
 	makeRandomCoords -selection $tempSel
 	rename $tempSel ""
    }

    foreach elem $pdbFileNames {
	if {! [file exists $elem]} {
	    error "PDB file $elem does not exist"
	}

	readPDB -fileName $elem
    }

    XplorCommand "vector do (mass  = 100.0) (all)" 
}


proc setupRgyr args {

    set restraintList [flagVal $args -restraints]
    set correction [flagVal $args -correction -1.0]

    foreach item $restraintList {
	
	if {[llength $item] == 1} {
	    set selPtr [lindex $item 0]
	    set nResidues [numResiduesInSelection $selPtr]
	    set rgTarget [expr (2.2 * pow($nResidues,0.38)) + $correction]
	} else {
	    set selPtr [lindex $item 0]
	    set rgTarget [lindex $item 1]
	}

	AtomSel tempSel -this $selPtr
	set selString [tempSel string]
	rename tempSel ""

	#
	# create an Rg restraint 
	# 

	
	XplorCommand \
	    "collapse 
                assign ($selString) 100.0 $rgTarget 
                scale 1.0 
             end"
	
    }
}


proc resetTemperature {newTemp} {

    # FIX:  change to an AtomSelAction
    for {set i 0} {$i < [xplorSim numAtoms]} {incr i} {

	set dev [expr sqrt($newTemp * [kBoltzmann] / [xplorSim atomMass $i] )]
	set vx  [expr [nextGaussianRandom] * $dev]
	set vy  [expr [nextGaussianRandom] * $dev]
	set vz  [expr [nextGaussianRandom] * $dev]
	
	xplorSim setAtomVel $i [list $vx $vy $vz]
    }
}
    

proc resetCoordsAndTemperature args {

    set newTemp    [requiredFlagVal $args -temperature]
    set coordArray [requiredFlagVal $args -coords]

    xplorSim setAtomPosArr $coordArray
    resetTemperature $newTemp

    marvinPyth command "d.init()"
}


proc distanceMovedFrom {anAtomPosArray} {

    #
    # returns the RMSD relative to the given atomPosArray
    # Doesn't mess up the main or comparison coordinate sets
    #

    #
    # Note that if the comparison set is undefined, coor swap spits up a COOR-ERR 
    # statement.  This should shut it up.
    #

    XplorCommand "set message off echo off print off end"

    XplorCommand "vector show elem (x) (known)"
    if {[XplorVariableNamed SELECT] > 0} {
	
	XplorCommand "coor swap end"
    } else {
	XplorCommand "coor copy end"
    }

    set curCompAPA [xplorSim atomPosArr]

    xplorSim setAtomPosArr $anAtomPosArray

    XplorCommand "coor fit end"
    XplorCommand "coor rms end"
    set rmsd [XplorVariableNamed result]
    
    xplorSim setAtomPosArr $curCompAPA
    XplorCommand "coor swap end"
   
    XplorCommand "set message on echo on print on end"

    return $rmsd
}


#
# ensemble cluster analysis, bassed on Kelley et al., Prot. Eng. 9, 1063-1065 (1996)
#    

proc interClusterDistance {aCluster bCluster pairwiseRMSDArrayName} {

    #
    #  Kelley et al., first eqn
    #

    upvar $pairwiseRMSDArrayName pairwiseRMSD

    set sumDist 0
  
    foreach aName $aCluster {
	foreach bName $bCluster {
	    set sumDist [expr $sumDist + $pairwiseRMSD([join [list $aName $bName]])]
	}
    }

    set retVal [expr $sumDist / double([llength $aCluster] * [llength $bCluster])]
}


proc clusterSpread {aCluster pairwiseRMSDArrayName} {

    upvar $pairwiseRMSDArrayName pairwiseRMSD

    #
    # Given a list of PDB files that are a cluster, and the name of the pairwise RMSD array, 
    # calculate the average RMSD of all pairs of structures in the cluster 
    # (ignoring pairing a structure with itself)
    #

    set sumDist 0

    foreach iName $aCluster {
	foreach kName $aCluster {
	    if {$iName != $kName} {
		set sumDist [expr $sumDist + $pairwiseRMSD([join [list $iName $kName]])]
	    }
	}
    }

    set n [llength $aCluster]
    if {$n > 0} {
	set retVal [expr $sumDist / double(($n * $n) - $n)]
    } else {
	set retVal 0
    }
    return $retVal
}


proc avgSpread {clusterList pairwiseRMSDArrayName} {

    #
    # the paper says not to include clusters of size 1, but I don't see how you could
    # calculate anything for the first stage of clustering without them
    #


    set nGoodClusters 0

    foreach curCluster $clusterList {
	incr nGoodClusters
	set sumSpread [expr $sumSpread + [clusterSpread $curCluster $pairwiseRMSDArrayName]]
    }
    
    set retVal [expr $sumSpread / $nGoodClusters]
    return $retVal
}

proc cluster args {

    #
    # This follows Kelley et al, Prot Eng 9, 1063-1065 (1996)
    #

    set files [requiredFlagVal $args -files]
    set selString [flagVal $args -selString "(name ca or name c or name n)"]

    #
    # eval RMSD between all pairs of structures
    #

    foreach fileA $files {
	set nameA [lindex $fileA 0]
	xplorSim setAtomPosArr [lindex $fileA 1]
	XplorCommand "coor copy end"

	foreach fileB $files {
	    set nameB [lindex $fileB 0]
	    xplorSim setAtomPosArr [lindex $fileB 1]

	    XplorCommand [format "coor fit sele %s end" $selString]
	    XplorCommand [format "coor rms sele %s end" $selString]
	    set pairwiseRMSD([join [list $nameA $nameB]]) [XplorVariableNamed "RMS"]
	}
    }

    
    #
    # start with each structure as its own cluster
    #

    set curClusters [list]
    foreach file $files {
	lappend curClusters [lindex $file 0]
    }

    
    #
    # 
}
