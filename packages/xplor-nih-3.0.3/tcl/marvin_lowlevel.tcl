#
# low level stuff for MARVIN
#
# JJK 3/14/02
#

package provide marvin 1.0

package require pasdpot
package require marvinassignment

#
# used for silencing procs that echo lots of output to the screen--
# just call verboseProc .... ; noOutput
#

proc noOutput {} {
}

#
# return either the next character from fileID, 
# or an error if it's at EOF
#

proc nextChar {fileID} {

    if {[eof $fileID]} {
	error "at EOF in nextChar"
    } else {
	return [read $fileID 1]
    }
}

proc nextWord {fileID} {

    global errorInfo

    #
    # eat whitespace from the current position in this file
    # until we find non-whitespace or the EOF
    #

    # if we hit EOF before the end of the whitespace, 
    # the untrapped error from nextChar quits us out

    while {1} {

	set curChar [nextChar $fileID]
	
	if {![isWhitespace $curChar]} {
	    break
	}
    }

    #
    # continue reading chars and add them to curWord
    # until we hit a whitespace or EOF
    #
    
    set curWord $curChar
        
    while {1} {
	
	if {[catch {set curChar [nextChar $fileID]}]} {
	    set errorInfo ""
	    break
	}
	
	if {[isWhitespace $curChar]} {
	    
	    break
	} 

	set curWord "$curWord$curChar"
    }
    
    #
    # return the word
    #
    
    return $curWord
}

proc nextLine {fileID} {

    global errorInfo

    set curLine ""
    set done 0

    # make sure we have some text to return

    if {[eof $fileID]} {
	error "at EOF in nextLine"
    }

    set curChar ""

    while {1} {

	set prevChar $curChar

	if {[catch {set curChar [nextChar $fileID]}]} {
	    set errorInfo ""
	    break
	}

	if {($curChar == "\n") && ($prevChar != "\\")} {

	    break
	}

	if {($curChar == "\n") && ($prevChar == "\\")} {
	    set curChar " "
	}
	
	set curLine "$curLine$curChar"
    }

    return $curLine
}   

proc isNumber {aWord} {

    global errorInfo
    
    if {[catch {expr double($aWord)}]} {
	set errorInfo ""
	return 0
    } else {
	return 1
    }
}

proc isInteger {aWord} {

    if {![isNumber $aWord]} {
	return 0
    } else {

	if {$aWord == [expr int($aWord)]} {
	    return 1
	} else {
	    return 0
	}
    }
}

	    
proc isPosNumber {aWord} {

    if {![isNumber $aWord]} {
	return 0
    } elseif {$aWord < 0} {
	return 0
    } else {
	return 1
    }
}

proc isZeroToOneNumber {aWord} {

    if {![isNumber $aWord]} {

	return 0

    } elseif {$aWord < 0} {

	return 0

    } elseif {$aWord > 1} {

	return 0

    } else {

	return 1

    }
}

proc isWhitespace {aWord} {

    if {[string trim $aWord] == ""} {
	return 1
    } else {
	return 0
    }
}


proc nextSel {fileID} {

    #
    # read the next selection from the given file
    #

    #
    # eat leading text
    #

    set curChar ""

    while {1} {

	if {[catch {set curChar [nextChar $fileID]}]} {
	    error "Unexpected EOF searching for initial ( in nextSel"
	}
	
	if {$curChar == "("} {
	    break
	}
    }

    set selectionText $curChar
    set numOpenParens 1

    #
    # keep eating characters until the count of
    # open parens reaches zero
    #

    while {($numOpenParens > 0)} {

	if {[catch {set curChar [nextChar $fileID]}]} {
	    error "Unexpected EOF in nextSel"
	}
	    
	set selectionText [format "%s%s" $selectionText $curChar]
	    
	if {$curChar == "("} {
		
	    incr numOpenParens
	}

	if {$curChar == ")"} {
		
	    incr numOpenParens -1
	}
    }

    return $selectionText

}


proc isAFlag { target } {

    # a flag is any string that starts with - 
    # but is not a number and is not a list
    
    if { ([llength $target] == 1) && ([regexp {^\-} $target]) && (! [string is double -strict $target]) } {
	
	return 1
	
    } else {

	return 0
    }
}


proc flagExists {flagList targetList} {

    foreach target $targetList {
	if {[lsearch -exact $flagList $target] != -1} {
	    return 1
	} 
    }

    return 0
}

proc requiredFlagVal args {

    set flagList [lindex $args 0]
    set target   [lindex $args 1]

    if {![flagExists $flagList $target]} {
	error "Missing required flag $target"
    } 
    
    return [flagVal $flagList $target]
}

proc flagVal args {

    if {[llength $args] == 2} {

	set flagList   [lindex $args 0]
	set targetList [lindex $args 1]
	set defVal     {}

    } elseif {[llength $args] == 3} {

	set flagList   [lindex $args 0]
	set targetList [lindex $args 1]
	set defVal     [lindex $args 2]

    } else {

	error "flagVal must be called with either 2 or 3 parameters"
    }


    set ret ""

    # get the target flag's position

    foreach target $targetList {
	set loc [lsearch -exact $flagList $target]
	if {$loc != -1} {
	    break
	}
    }

    if {$loc != -1} {

	# starting at the next position, 
	# grab items from flagList until
	# we get to another flag or the end
	
	# FIX?  Not sure this is a good behavior, since it turns 
	# -flag x y into a return value of xy.  Better to just return
	# a value of x?  Note that this doesn't affect structures like 
	# -flag "a string value" or -flag [list 1 2 3].  In those cases, 
	# the text following the flag is a single entity in the flag list

	incr loc

	while {(! [isAFlag [lindex $flagList $loc]]) && ($loc < [llength $flagList])} {

	    append ret [lindex $flagList $loc]
	    incr loc
	}
    }

    if {$ret == {}} {
	set ret $defVal
    }
    
    return $ret
}


proc chooseRandomItem {aList} {

    set numItems [llength $aList]
    set randomIndex [randomIntegerBetween 0 [expr $numItems - 1]]

    return [lindex $aList $randomIndex]
}


proc origNextXplorChar {fileID} {

    # 
    # returns the next non-xplor-comment character
    #
    # if we hit EOF inside of either of this proc's while loops,
    # the untrapped error from nextChar quits us out
    #

    set curChar [nextChar $fileID]
    
    #
    # rules for the first xplor comment type: ! .. \n
    #

    if {$curChar == "!"} {
	while {$curChar != "\n"} {
	    set curChar [nextChar $fileID]
	} 

	#
	# at this point, curChar is the newline.
	# Eat it and return the next non-xplor-comment char
	#
	set curChar [origNextXplorChar $fileID]
    }

    #
    # rules for the second xplor comment type: { ... }
    #

    if {$curChar == "\{"} {
	while {$curChar != "\}"} {
	    set curChar [nextChar $fileID]
	}

	#
	# at this point, curChar is the close-brace.
	# Eat it and return the next non-xplor-comment char
	#
	set curChar [origNextXplorChar $fileID]
    }

    return $curChar
}


proc nextXplorChar {fileID} {

    #
    # returns the next non-xplor-comment character
    #
    # if we hit EOF inside of this proc, the
    # untrapped error from nextChar quits us out
    #
    # Re-written to avoid nested calls to nextXplorChar, 
    # which can cause crashes in skipping over very long
    # groups of comments
    #

    while {1} {

	set curChar [nextChar $fileID]
	set curCharOK 1

	#
	# handle the first xplor comment type: ! .. \n
	#
	
	if {$curChar == "!"} {
	    while {$curChar != "\n"} {
		set curChar [nextChar $fileID]
	    } 
	    
	    # at this point, curChar is the newline
	    set curCharOK 0
	}

	#
	# handle the second xplor comment type: { ... }
	#
	
	if {$curChar == "\{"} {
	    while {$curChar != "\}"} {
		set curChar [nextChar $fileID]
	    }
	    
	    # at this point, curChar is the close-brace
	    set curCharOK 0
	}

	if {$curCharOK} {
	    break
	}
    }

    return $curChar
}


proc nextXplorWord {fileID} {

    global errorInfo

    #
    # eat whitespace from the current position in this file
    # until we find non-whitespace or the EOF
    #

    # if we hit EOF before the end of the whitespace, 
    # the untrapped error from nextChar quits us out

    while {1} {

	set curChar [nextXplorChar $fileID]
	
	if {![isWhitespace $curChar]} {
	    break
	}
    }

    #
    # continue reading chars and add them to curWord
    # until we hit a whitespace or EOF
    #
    
    set curWord $curChar
        
    while {1} {
	
	if {[catch {set curChar [nextXplorChar $fileID]}]} {
	    set errorInfo ""
	    break
	}
	
	if {[isWhitespace $curChar]} {
	    
	    break
	} 

	set curWord "$curWord$curChar"
    }
    
    #
    # return the word
    #
    
    return $curWord
}


proc eatXplorTextUpTo {fileID targetWords} {

    set foundTarget 0

    while {$foundTarget == 0} {

        set curWord [nextXplorWord $fileID]

	if {[lsearch -exact $targetWords $curWord] != -1} {
            set foundTarget 1
        } 

        if {[eof $fileID]} {
	    error "unexpected EOF while searching for xplor text $targetWords"
        }
    }
}
	

proc eatNormalTextUpTo {fileID targetWords} {

    set foundTarget 0

    while {$foundTarget == 0} {

        set curWord [nextWord $fileID]

	if {[lsearch -exact $targetWords $curWord] != -1} {
            set foundTarget 1
        } 

        if {[eof $fileID]} {
	    error "unexpected EOF while searching for $targetWords"
        }
    }
}

proc flattenList {a} {

    set res [list]
    
    foreach item $a {

	if {[llength $item] == 1} {
	    lappend res $item
	} else {
	    set flatItem [flattenList $item] 
	    foreach b $flatItem {
		lappend res $b
	    }
	}
    }
    return $res
}

proc max args {

    if {[llength $args] == 0} {
	error "Attempt to calculate the max of an empty list"
    }

    set vals [flattenList $args]

    set curMax [lindex $vals 0]
    
    foreach val $vals {

	if {$val > $curMax} {
	    set curMax $val
	}
    }

    return $curMax
}

proc min args {

    if {[llength $args] == 0} {
	error "Attempt to calculate the min of an empty list"
    }

    set vals [flattenList $args]

    set curMin [lindex $vals 0]
    
    foreach val $vals {

	if {$val < $curMin} {
	    set curMin $val
	}
    }

    return $curMin
}

proc sum {aList} {

    set sum 0

    foreach val $aList {
	set sum [expr $sum + $val]
    }
    return $sum
}

proc mean {aList} {

    if {[llength $aList] == 0} {
	error "Attempt to calculate the mean of an empty list"
    }

    set sum [sum $aList]
    return [expr $sum / double([llength $aList])]
}

proc standardDeviation {aList} {

    if {[llength $aList] == 0} {
	error "Attempt to calculate the std dev of an empty list"
    }

    if {[llength $aList] == 1} {
	return 0
    }

    set m [mean $aList]
    set sum 0
    
    foreach val $aList {
	
	set cur [expr pow(($val - $m),2)]
	set sum [expr $sum + $cur]
    }
    
    set variance [expr $sum / double([llength $aList] - 1)]
    set dev [expr sqrt($variance)]
    
    return $dev
}

proc median {aList} {
        
    if {[llength $aList] == 0} {
	error "Attempt to calculate the median of an empty list"
    }

    set l [lsort -real $aList]

    if {[expr fmod([llength $l], 2)] == 0} {
	
	
	set a [lindex $l [expr [llength $l] / 2]]
	set b [lindex $l [expr ([llength $l] / 2) + 1]]
	set med [expr ($a + $b) / 2]
	       
    } else {
	
	set med [lindex $l [expr int(ceil([llength $l] / double(2)))]]

    }
    
    return $med
}

proc firstItem {l} {

    if {[llength $l] == 0} {
	error "Attempt to get the first item of an empty list"
    }

    return [lindex $l 0]
}

proc lastItem {l} {

    if {[llength $l] == 0} {
	error "Attempt to get the last item of an empty list"
    }

    return [lindex $l [expr [llength $l] - 1]]
}

proc nextGaussianRandom args {

    set targetMean [flagVal $args -mean 0.0]
    set targetDev  [flagVal $args -dev  1.0]

    # from Numerical Recipes, section 7.2

    # uniformRandom is a call to the xplor random number generator,
    # defined in SimulationWorld.i

    set r 0

    while {($r >= 1.0) || ($r == 0.0)} {

	set v1 [expr 2.0 * [uniformRandom] - 1.0]
	set v2 [expr 2.0 * [uniformRandom] - 1.0]
	set r  [expr ($v1 * $v1) + ($v2 * $v2)]
    }

    set fac [expr sqrt(-2.0 * log($r) / $r)]

    set gaussRandNum [expr $v2 * $fac]

    return [expr $targetMean + ($gaussRandNum * $targetDev)]
}

proc randomIntegerBetween {a b} {

    #
    # returns a random integer between integers a and b, inclusive
    # flat distribution
    #

    if {$a < $b} {
	set begin $a
	set finish $b
    } else {
	set begin $b
	set finish $a
    }

    set range [expr ($finish - $begin) + 1]
    set randVal [expr int([uniformRandom] * $range) + $begin]
    if {$randVal > $finish} {
	set randVal $finish
    }
    
    return $randVal
}


proc shuffledList {aList} {

    #
    # Given a list, return it with its elements in random order
    #

    set nElements [llength $aList]

    for {set curPos 0} {$curPos < $nElements} {incr curPos} {

	# choose a position between curPos and the end of the list

	set otherPos [randomIntegerBetween $curPos [expr $nElements - 1]]

	# swap contents of curPos and otherPos

	set temp [lindex $aList $curPos]
	set aList [lreplace $aList $curPos $curPos [lindex $aList $otherPos]]
	set aList [lreplace $aList $otherPos $otherPos $temp]
    }

    return $aList
}


#
# support for the histogram proc
#

proc whichBin {minVal maxVal nBins aVal} {


    #
    # defend against zero binWidth arising from maxVal == minVal
    #

    if {$maxVal == $minVal} {
	set binWidth 1
    } else {
	set binWidth [expr ($maxVal - $minVal) / double($nBins)]
    }

    set aBin [expr floor(($aVal - $minVal) / $binWidth)]

    #
    # FIX:  Make these range enforcements optional, or move them elsewhere
    #

    if {$aBin < 0} {
	set aBin 0
    }

    if {$aBin >= $nBins} {
	set aBin [expr $nBins -1]
    }

    return [expr int($aBin)]
}


proc binCenter {minVal maxVal nBins whichBin} {

    #
    # Given a range minVal .. maxVal, 
    # divided into nBins bins 
    # which are numbered 0 .. nBins - 1, 
    # Return the value of the center of bin #whichBin
    #

    #
    # defend against maxVal == minVal
    #

    if {$maxVal == $minVal} {
	return $maxVal
    } 

    set binWidth [expr ($maxVal - $minVal) / double($nBins)]

    set center [expr $minVal + ($whichBin * $binWidth) + (0.5 * $binWidth)]
    return $center
}


proc histogram args {

    set data                [requiredFlagVal $args -data]
    set minVal              [flagVal $args -min [min $data]]
    set maxVal              [flagVal $args -max [max $data]]
    set nBins               [flagVal $args -nBins 20]
    set nRows               [flagVal $args -nRows 20]
    set title               [flagVal $args -title "Histogram"]
    set nBlankLinesToFollow [flagVal $args -numBlankLinesToFollow 2]

    if {[flagExists $args -integerBins]} {
	set nBins [expr int(ceil($maxVal) - floor($minVal)) + 1]
    }
    
    for {set binCount 0} {$binCount < $nBins} {incr binCount} {
	set numInBin($binCount) 0
	set fracInBin($binCount) 0
    }

    foreach datum $data {
	incr numInBin([whichBin $minVal $maxVal $nBins $datum])
    }

    set totData [llength $data]

    if { $totData==0 } { return "" }

    for {set binCount 0} {$binCount < $nBins} {incr binCount} {
	set fracInBin($binCount) [expr $numInBin($binCount) / double($totData)]
    }

    set fracPerRow [expr 1.0 / double($nRows)]

    for {set binCount 0} {$binCount < $nBins} {incr binCount} {
	set nCols($binCount) [expr round($fracInBin($binCount) / $fracPerRow)]
    }


    #
    # draw the bottom axis
    #

    set curRow "0.0 +"
    for {set binCount 0} {$binCount < $nBins} {incr binCount} {
	set curRow [format "%s-" $curRow]
    }

    set histPic $curRow

    #
    # now draw the histo
    #

    for {set rowCount 0} {$rowCount < $nRows} {incr rowCount} {

	#
	# draw the vertical axis, with a label at the top
	#

	if {$rowCount < [expr $nRows - 1]} {
	    set curRow "    |"
	} else {
	    set curRow "1.0 |"
	}

	#
	# draw the current line
	#

	for {set binCount 0} {$binCount < $nBins} {incr binCount} {
	    if {$nCols($binCount) > $rowCount} {
		set curRow [format "%sX" $curRow]
	    } else {
		set curRow [format "%s " $curRow]
	    }
	} 
	set histPic [format "%s\n%s" $curRow $histPic]
    }
    
    #
    # add the vert axis label
    #

    set histPic [format "fraction of\n%d datapoints\n%s" $totData $histPic]

    #
    # add the title
    #

    set nSpaces [expr ($nBins - [llength $title]) / 2]
    for {set count 0} {$count < $nSpaces} {incr count} {
	set title [format " %s" $title]
    }

    set histPic [format "%s\n\n%s" $title $histPic]

    #
    # add labels to the bottom axis
    # Integer values need both the check and the 
    # rounding to avoid problems with eg., minVal = 4.0
    #

    if {[isInteger $minVal]} {
	set minAsString [format "%d" [expr round($minVal)]]
    } else {
	set minAsString [format "%f" $minVal]
    }

    if {[isInteger $maxVal]} {
	set maxAsString [format "%d" [expr round($maxVal)]]
    } else {
	set maxAsString [format "%f" $maxVal]
    }

    set nSpaces [expr $nBins - [string length $minAsString] - [string length $maxAsString]]
    set nSpaces [max [list 1 $nSpaces]]

    set botLine [format "     %s" $minAsString]
    for {set count 0} {$count < $nSpaces} {incr count} {
	set botLine [format "%s " $botLine]
    }
    set botLine [format "%s%s" $botLine $maxAsString]

    set histPic [format "%s\n%s" $histPic $botLine]

    for {set count 0} {$count < $nBlankLinesToFollow} {incr count} {
	set histPic [format "%s\n" $histPic]
    }

    return $histPic
}

proc removeIfExists {fname} {
    if {[file exists $fname]} {
	exec rm $fname
    }
}

proc listIncludes {aList anItem} {

    if {[lsearch -exact $aList $anItem] == -1} {
        return 0
    } else {
        return 1
    }
}


#
# this is used for PIPP and TALOS filereading
#

# FIX:  what about DNA residues?
	
proc oneToThreeLetterResName {aName} {

    set aName [string tolower $aName]

    switch $aName {

	g {return gly}
	a {return ala}
	v {return val}
	i {return ile}
	l {return leu}
	f {return phe}
	p {return pro}
	m {return met}
	w {return trp}
	c {return cys}
	s {return ser}
	t {return thr}
	n {return asn}
	q {return gln}
	y {return tyr}
	h {return his}
	d {return asp}
	e {return glu}
	k {return lys}
	r {return arg}
	default {error "Unrecognized 1-letter residue name: $aName"}
    }
}


#
# A convenient place to put an algorithm for generating
# random seeds as a function of the Marvin pass, the 
# structure number, and an overall (user set) seed
#

proc setMarvinRandomSeed args {

    if {[llength $args] == 1} {

	setRandomSeed [lindex $args 0]

    } else {

	set overall [flagVal $args -overall 0]
	set pass    [flagVal $args -pass 1]
	set count   [flagVal $args -count 0]

	# NOTE: this algorithm limits the number of structs per pass!

	set maxStructsPerPass 10000
	set maxPasses 10

	set temp [expr ($overall * $maxPasses * $maxStructsPerPass) \
		  + ($pass * $maxStructsPerPass) \
		  + $count]
	
	setRandomSeed $temp
    }
}


#
# Since I'm so often building up long strings in this way, 
# I might as well speed it up
#

proc appendLineToString {stringVarName newLine} {

    upvar $stringVarName s 

    if {$s == ""} {
	set s $newLine
    } else {
	set s [format "%s\n%s" $s $newLine]
    }
}


proc sign {x} {

    if {$x > 0} {
	return 1
    } else {
	return -1
    }
}

 
#
# A general update-user routine.
# Off by default.  Turned on via the xplor
# command line flag -verbose.
# By default, doesn't attach a newline after each
# message, so updates can just overwrite previous ones
# by ending in a \r.  
#

proc updateUser {msg} {

    global hasVerboseFlag
    global argv
    global updateUserNewline

    if {! [info exists hasVerboseFlag]} {
	set hasVerboseFlag [expr [flagExists $argv -verbose] || [flagExists $argv -verb]] 
    }

    if { !$hasVerboseFlag} {
	return
    }

    if {! [info exists updateUserNewline]} {
	set updateUserNewline 0
    }

    if { $updateUserNewline } {
	puts stdout $msg
    } else {
	puts -nonewline stdout $msg
    }
}


#
# Given a list of lists, return a list of every permutation generated 
# by choosing one item of each element of the input list
#

proc listPermutations {l} {

    proc addPermutation {startPermutations newElems} {

	#
	# Given a starting list of permutations, and a list of 
	# new elements, generate a new list, each element of which
	# is an element of the starting list of permutations fused 
	# with one element of the list of new elements.
	#
	# If the list of starting permutations is empty, just
	# return the list of new elements.
	#
	
	if {[llength $startPermutations] == 0} {
	    return $newElems
	} else {

	    set retVal [list]
	    
	    foreach startPerm $startPermutations {
		foreach newElem $newElems {
		    set newPerm [concat $startPerm $newElem]
		    lappend retVal $newPerm
		}
	    }
	    return $retVal
	}
    }

    set allPermutations [list]

    foreach elem $l {
	set allPermutations [addPermutation $allPermutations $elem]
    }

    return $allPermutations
}
   
