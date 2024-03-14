#!/bin/bash
protname="$1"

CREATE_CONSTRAINT_CMD="../packages/find_viols/createCheckViol.sh"
FIND_VIOL_CMD="../packages/find_viols/runXplorForViols.sh"

alreadyrun="../proteinsparams/listoffiles.txt"
pdb_folder="finalpdbfiles"
datafolder="../protein"
violfiles="distviols"

xplorfolder="waterrefine"

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
	protfolder=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`	
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
			mkdir $datafolder/$protfolder/$violfiles
			#echo "mkdir $datafolder/$protfolder/$violfiles"
			for pdbs in `ls $pdbfilesfolder/*pdb`
			do
				filename=`basename $pdbs`
				filename=${filename::-4}
				constraintrunfile="$datafolder/$protfolder/$violfiles/$filename.py"
			
				$CREATE_CONSTRAINT_CMD $pdbs $xplornoe $constraintrunfile
				$FIND_VIOL_CMD $constraintrunfile $filename
				#echo "$CREATE_CONSTRAINT_CMD $pdbs $xplornoe $constraintrunfile"
				#echo "$FIND_VIOL_CMD $constraintrunfile $filename"
			done		
		else
			echo "XPLOR_NOE not found"
		fi	
	else
		echo "$paramfile NOT found"
		exit -1	
	fi
else
	protfolder=$protname"_allparam"	
	pdbfilesfolder=$datafolder/$protfolder/$pdb_folder	
	paramfile=$datafolder"/"$protname"/"$protname".xml"
	
	if [ -f $paramfile ]; then
		
	        xplornoe_=`xmllint --xpath "string(//upl_xplor)" "$paramfile"`
		xplordihed_=`xmllint --xpath "string(//aco_xplor)" "$paramfile"`			
		if [ -z $xplornoe_ ]; then		
			# assumed that makeXplorUplAco.sh has already ran
			xplornoe=`ls ../protein/$protname/*noe.tbl`
			xplordihed=`ls ../protein/$protname/*aco.tbl`
		else
			xplornoe=$datafolder/$protname/$xplornoe_
			if [ ! -z $xplordihed_ ]; then
				xplordihed=$datafolder/$protname/$xplordihed_
			else
				xplordihed=""
			fi
		fi
		
		if [ ! -z $xplornoe ]; then		
			mkdir $datafolder/$protfolder/$violfiles
			#echo "mkdir $datafolder/$protfolder/$violfiles"
			for pdbs in `ls $pdbfilesfolder/*pdb`
			do
				filename=`basename $pdbs`
				filename=${filename::-4}
				constraintrunfile="$datafolder/$protfolder/$violfiles/$filename.py"
				
				filename_folder=${filename::-5}
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
			echo "XPLOR NOE file not found"
		fi	
	else
		echo "$paramfile NOT found"
		exit -1	
	fi	
fi
