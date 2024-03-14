#
# TCL code for reading NOE peak tables in 
# nmrDraw format
#

package provide nmrdraw 1.0
package require marvin
package require pipp

#
# nmrDraw peak tables are virtually identical to PIPP PCK tables,
# so I just have these wrappers to change the name
#

proc process3dCnmrDrawPeakTable args {

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

proc process3dNnmrDrawPeakTable args {

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

proc process3dnmrDrawPeakTable args {

    set remVar [flagVal $args -remarksVariableName ""]
    if {$remVar != ""} {
	upvar $remVar $remVar
    }

    eval process3dPippPeakTable $args
    return ""
}
