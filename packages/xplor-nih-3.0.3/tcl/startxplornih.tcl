#	  
# startxplornih.tcl
#	   
# System-wide TCL commands for NIH xplor.  
# This file is sourced during startup of the TCL interp.
#	 
# Not intended for users to play with!
#	  
# John Kuszewski 5/1/02
#

#
# check if this is running in a tclsh process, or an xplor -tcl process
# set in xplorNIH.i
#

#clear errors from startup
set errorInfo ""

#proc unknown args {  }
#proc history args {  }

if {[info exists xplorNIHflag]} {
    unset xplorNIHflag
} else {

    # do something to prevent package require xplornih from doing anything!
}
	   
#	  
# load in the standard packages
#	  
	   
package require simulationworld
package require cdsvector
package require potlist
package require vec3
package require atom
package require atomsel
package require atomselaction
package require pdbtool
package require derivlist
	   
#	  
# create the interface to the XplorWrap::shell and XplorWrap::command 
# methods 
#	  
	   
proc XplorShell {} {
	   
   global wrapPtr
	   
   return [::xplorwrap::XplorWrap_shell $wrapPtr]
}	  
	   
proc XplorCommand {commandString {varList {}}} {
	   
   global wrapPtr
	   
   return [::xplorwrap::XplorWrap_command $wrapPtr $commandString $varList]
}	  

#	  
# create the global Simulation instance
# that's connected to the (static) XplorSimulation
#	  
# Note that in reality, xplorSimulationPtr is actually
# a ptr to a (static) XplorSimulation instance, 
# which in turn points to the (static) XplorSimulation
# instance
#	  
	   
XplorSimulation xplorSim -this $xplorSimulationPtr


#
# get the value of a given xplor-side variable or variables
#

proc XplorVariableNamed {l} {

    return [XplorCommand "" $l]
}
