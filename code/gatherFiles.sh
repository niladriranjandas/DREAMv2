#!/bin/bash
protname="$1"

if [ -z $protname ]; then
	echo "arg-1 needed"
	exit -1
fi

alreadyrun="../proteinsparams/listoffiles.txt"
waterref_folder="waterrefine"
proteinfolder="../protein"

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
	protfolder=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`
	
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
else
	#echo "code to be done"
	waterref_whichfiles=$protname"_whichwaterref_mat.txt"
	waterref_which_folder="gap_correct/RESIDUAL_FILES/$protname"
	
	waterref_folder="waterrefine"
	MAXFILE=5
	
	protfolder=$protname"_allparam"
	# make folder in protein folder #
	finalfolder=$proteinfolder/$protfolder/finalpdbfiles
	mkdir $finalfolder	
	
	for whichfiles in `head -$MAXFILE $waterref_which_folder/$waterref_whichfiles`
	do
		foldername=`basename $whichfiles`
		foldername=${foldername::-4}
		
		waterref_soln=`ls $waterref_folder/$foldername/*wref`
		newname=$waterref_soln".pdb"
		grep -v TIP3 $waterref_soln > "$newname"
		
		# make folder in protein folder #
		cp "$newname" $finalfolder/.				
	done
fi
