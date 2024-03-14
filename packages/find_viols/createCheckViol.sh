#!/bin/bash

pdbfile="$1"
noefile="$2"
opfilename="$3"
templatefile="../packages/find_viols/template_see_constraints.py"


if [ ! -f $pdbfile ]; then
	echo "$pdbfile NOT found"
	exit -1
fi
if [ ! -f $noefile ]; then
	echo "$noefile NOT found"
	exit -1
fi

# ---
        pdbnametmp=${pdbfile//\//\\/}
        sed "s/PDB_NAME/$pdbnametmp/g" $templatefile > $opfilename

        noefiletmp=${noefile//\//\\/}
        sed -i "s/NOE_CONSTRAINS/$noefiletmp/g" $opfilename

#echo "s/PDB_NAME/$pdbnametmp/g $templatefile > $opfilename"
