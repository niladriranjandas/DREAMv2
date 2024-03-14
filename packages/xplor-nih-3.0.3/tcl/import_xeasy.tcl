#
# tools for importing XEASY formatted files
#
#
# XEASY table headers look like 
#
#
# Number of dimensions 3
# FORMAT xeasy3D
#INAME 1 H
#INAME 2 C
#INAME 3 HC
#CYANAFORMAT hCH


package require marvin
package provide xeasy 0.1

proc readXeasyPeakTable args {

    set fname [requiredFlagVal $args -fileName]
    set colNames  [requiredFlagVal $args -cols]

    global errorInfo

    set retVal [list]

    # open file

    if {[catch {set inUnit [open $fname r]}]} {
	error "Error opening input file $fname"
    }

    #
    # read the number of dimensions line
    #

    while {1} {

	if {[catch {set l [nextLine $inUnit]}]} {
	    error "Can't find Number of dimensions line in Xeasy file $fname"
	}

	if {[string match "*Number of dimensions*" $l]} {
	    set nDimensions [lindex $l end]
	    break
	}
    }


    set syntheticVarsLine [list "VARS" "PkID"]
    for {set x 1} {$x <= $nDimensions} {incr x} {
	lappend syntheticVarsLine [format "SHIFT%d" $x]
    }

    lappend syntheticVarsLine "JUNK"
    lappend syntheticVarsLine "JUNK"
    lappend syntheticVarsLine "Intensity"

    
    #
    # read the INAME comment lines, which define the order
    # of the shift position columns, and replace the SHIFT* 
    # entries in the synthetic VARS line
    #

    while {1} {

	if {[catch {set l [nextLine $inUnit]}]} {
	    error "EOF during search for headers in Xeasy file $fname"
	}


	if {[string first "\#" $l] != 0} {
	    break
	}

	set firstWord [lindex $l 0]

	if {$firstWord == "\#INAME"} {

	    set oldColName [format "SHIFT%d" [lindex $l 1]]
	    set newColName [lindex $l 2]
	    
	    set p [lsearch -exact $syntheticVarsLine $oldColName]
	    if {$p != -1} {
		set syntheticVarsLine [lreplace $syntheticVarsLine $p $p $newColName]
	    }
	}
    }
    


    set columnPosns [list]

    foreach colName $colNames {
	lappend columnPosns [expr [lsearch -exact $syntheticVarsLine $colName] - 1]
    }

    # make sure none of these are missing

    foreach pos $columnPosns colName $colNames {
	if {$pos == -2} {
	    error "Can't find column $colName in XEASY-style header of file $fname" 
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
	# skip commented lines
	#

	if {[lindex $l 0] == "\#"} {
	    continue
	}

	set curVal [list]

	#
	# try extracting out the columns' data
	#
	
	if {[catch {
	    
	    foreach pos $columnPosns {
		lappend curVal [lindex $l $pos]
	    }
	    
	}]} {
	    
	    error "cant parse line: $l from XEASY peak file $fname"
	}

	lappend retVal $curVal
    }

    puts stdout [format "%d shifts were read from XEASY formatted file %s" \
	   [llength $retVal] $fname]

    close $inUnit
    return $retVal
}



proc process3dXeasyPeakTable args {

    set fname         [requiredFlagVal $args -fileName]
    set pot           [requiredFlagVal $args -pot]
    set fromProtName  [requiredFlagVal $args -fromProtonColumnName]
    set fromHeavyName [requiredFlagVal $args -fromHeavyatomColumnName]
    set toProtName    [requiredFlagVal $args -toProtonColumnName]
    set namePrefix    [requiredFlagVal $args -namePrefix]
    set remVar        [flagVal $args -remarksVariableName ""]

    #
    # The column names for peak ID and intensity are hardwired, 
    # because they aren't actually defined in the XEASY header.  
    # Rather, they're synthesized in readXeasyPeakTable, to create something
    # that looks like an nmrPipe table, which is then read normally
    #

    set peaks [readXeasyPeakTable -fileName $fname -cols [list "PkID" "Intensity" $fromProtName $toProtName $fromHeavyName]]

    foreach peak $peaks {

	set curID   [lindex $peak 0]
	set curName [format "%s%s" $namePrefix $curID]
	set tempPeak [Peak -args $curName]
	$tempPeak setIntensity [lindex $peak 1]
	$tempPeak setFromProtonShift [lindex $peak 2]
	$tempPeak setToProtonShift [lindex $peak 3]
	$tempPeak setFromHeavyatomShift [lindex $peak 4]
	$tempPeak appendToNote [format "from XEASY file %s, peak %s" $fname $curID]
	
	$pot addPeak [$tempPeak cget -this]
	$tempPeak -disown
	rename $tempPeak ""
    }
    
    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "Created %d peaks from XEASY file %s" [llength $peaks] $fname]
    }
    puts stdout [format "Created %d peaks from XEASY file %s" \
		     [llength $peaks] $fname]

    return ""
}

proc process2dXeasyPeakTable args {

    set fname         [requiredFlagVal $args -fileName]
    set pot           [requiredFlagVal $args -pot]
    set fromProtName  [requiredFlagVal $args -fromProtonColumnName]
    set toProtName    [requiredFlagVal $args -toProtonColumnName]
    set namePrefix    [requiredFlagVal $args -namePrefix]
    set remVar        [flagVal $args -remarksVariableName ""]

    #
    # The column names for peak ID and intensity are hardwired, 
    # because they aren't actually defined in the XEASY header.  
    # Rather, they're synthesized in readXeasyPeakTable, to create something
    # that looks like an nmrPipe table, which is then read normally
    #

    set peaks [readXeasyPeakTable -fileName $fname \
		   -cols [list "PkID" "Intensity" $fromProtName $toProtName]]

    foreach peak $peaks {

	set curID   [lindex $peak 0]
	set curName [format "%s%s" $namePrefix $curID]
	set tempPeak [Peak -args $curName]
	$tempPeak setIntensity [lindex $peak 1]
	$tempPeak setFromProtonShift [lindex $peak 2]
	$tempPeak setToProtonShift [lindex $peak 3]
	$tempPeak appendToNote [format "from XEASY file %s, peak %s" $fname $curID]
	
	$pot addPeak [$tempPeak cget -this]
	$tempPeak -disown
	rename $tempPeak ""
    }
    
    if {$remVar != ""} {
	upvar $remVar tempRem
	lappend tempRem [format "Created %d peaks from XEASY file %s" [llength $peaks] $fname]
    }
    puts stdout [format "Created %d peaks from XEASY file %s" \
		     [llength $peaks] $fname]

    return ""
}
