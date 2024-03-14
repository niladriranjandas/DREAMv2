#!/bin/bash
protname="$1"

CREATE_CONSTRAINT_CMD="../packages/find_viols/createCheckViol.sh"
FIND_VIOL_CMD="../packages/find_viols/runXplorForViols.sh"

alreadyrun="../proteinsparams/listoffiles.txt"
pdb_folder="finalpdbfiles"
datafolder="../protein"
violfiles="distviols"

xplorfolder="waterrefine"
# ---------------------------- #

protfolder="$protname"
pdbfilesfolder=$datafolder/$protfolder/$pdb_folder	
paramfile=$datafolder"/"$protfolder"/"$protfolder".xml"	

if [ -f $paramfile ]; then
	
        xplornoe_=`xmllint --xpath "string(//upl_xplor)" "$paramfile"`
	xplordihed_=`xmllint --xpath "string(//aco_xplor)" "$paramfile"`			
	if [ -z $xplornoe_ ]; then		
		# assumed that makeXplorUplAco.sh has already ran
		xplornoe=`ls $datafolder/$protfolder/*noe.tbl`
		xplordihed=`ls $datafolder/$protfolder/*aco.tbl`
	else
		xplornoe=$datafolder/$protfolder/$xplornoe_
		if [ ! -z $xplordihed_ ]; then
			xplordihed=$datafolder/$protfolder/$xplordihed_
		else
			xplordihed=""
		fi
	fi
	
	if [ ! -z $xplornoe ]; then
		if [ ! -d $datafolder/$protfolder/$violfiles ]; then
			mkdir $datafolder/$protfolder/$violfiles
			#echo "mkdir $datafolder/$protfolder/$violfiles"
		fi
		for pdbs in `ls $pdbfilesfolder/*pdb`
		do
			filename=`basename $pdbs`
			filename=${filename::-4}
			constraintrunfile="$datafolder/$protfolder/$violfiles/$filename.py"
		
			filename_folder=${filename::-5} 				  # removing "_wref"
                        noefolder="$datafolder/$protfolder/$xplorfolder/$filename_folder"
			xplornoe=`ls $noefolder/*noe.tbl`
			if [ ! -z $xplornoe ]; then
				$CREATE_CONSTRAINT_CMD $pdbs $xplornoe $constraintrunfile
				$FIND_VIOL_CMD $constraintrunfile $filename
				#echo "$CREATE_CONSTRAINT_CMD $pdbs $xplornoe $constraintrunfile"
				#echo "$FIND_VIOL_CMD $constraintrunfile $filename"
			fi
		done		
	else
		echo "XPLOR_NOE not found"
	fi	
else
	echo "$paramfile NOT found"
	exit -1	
fi

# ---------------------------- #
