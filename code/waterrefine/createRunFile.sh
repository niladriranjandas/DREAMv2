#!/bin/bash

protname="$1"
pdbfile="$2"
noefile="$3"
dihedfile="$4" # this file is optional

opfilename=$protname"_wref.py"
oppattern=$protname"_wref"

################################################
if [ ! -f $pdbfile ]; then
	echo "No pdb file found"
	exit -1
fi

if [ ! -f $noefile ]; then
	echo "No noe file found"
	exit -1
fi

if [ ! -z $dihedfile ]; then
	if [ ! -f $dihedfile ]; then
		echo "No dihed file found"
		exit -1
	fi
fi
################################################
	if [ -z $dihedfile ]; then
		considerTemplate="template_upl.py"
	else
		considerTemplate="template_upl_dihed.py"
	fi

	pdbnametmp=${pdbfile//\//\\/}
	sed "s/PATH_TO_PDB/$pdbnametmp/g" $considerTemplate > $opfilename
	
	noefiletmp=${noefile//\//\\/}	
	sed -i "s/PATH_TO_NOE_TBL/$noefiletmp/g" $opfilename
	
	if [ ! -z $dihedfile ]; then
		dihedfiletmp=${dihedfile//\//\\/}
		sed -i "s/PATH_TO_DIHED_TBL/$dihedfiletmp/g" $opfilename
	fi
	
	sed -i "s/OP_PREFIX/$oppattern/g" $opfilename
