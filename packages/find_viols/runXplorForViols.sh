#!/bin/bash

xplorseeconstrainfile="$1"
protname="$2"

xplorcmd="../packages/xplor-nih-3.0.3/bin/xplor -py"

if [ ! -f $xplorseeconstrainfile ]; then
	echo "$xplorseeconstrainfile NOT found"
	exit -1
fi

#--
	folder=`dirname $xplorseeconstrainfile`
	opfilename=$folder"/"$protname".viols"
	$xplorcmd $xplorseeconstrainfile > $opfilename

	opfilename2=$folder"/"$protname".onlyviols"
	awk '/--------/,EOF { getline; print $0; }' $opfilename > $opfilename2
