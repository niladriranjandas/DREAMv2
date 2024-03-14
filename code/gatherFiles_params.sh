#!/bin/bash
protname="$1"

if [ -z $protname ]; then
	echo "arg-1 needed"
	exit -1
fi

alreadyrun="../proteinsparams/listoffiles.txt"
waterref_folder="waterrefine"
proteinfolder="../protein"

# ---------------- #
protfolder="$protname"

# make folder in protein folder #
finalfolder=$proteinfolder/$protfolder/finalpdbfiles
mkdir $finalfolder

for pdbs in `ls $waterref_folder/$protfolder*/$protfolder*wref`
do
	newname=$pdbs".pdb"
	pdbname=`basename $pdbs`
	echo "$pdbname"
	#echo "grep -v TIP3 $pdbs > $newname"
	grep -v TIP3 "$pdbs" > "$newname"
	
	# make folder in protein folder #
	cp "$newname" $finalfolder/.
done
# ---------------- #