#!/bin/bash

protname="$1"
email="$2"

#SENDMAIL_CMD="python front_end/sendMail.py"
SENDMAIL_CMD="python front_end/sendMail_v2.py"
JOINPDB_CMD="pymol -ca ../packages/make_ensemble/alignNsaveEnsemble_v2.py -- "

alreadyrun="../proteinsparams/listoffiles.txt"
datafolder="../protein"
violfiles="distviols"
pdbfolder="finalpdbfiles"

progressfile="progress_msgs/sendfiles_list.runstatus.pdf"

emailfolder="sendfiles"
emailfileslist="sendfiles_list.txt"

# ------------- #

protfolder="$protname"
violsfolder=$datafolder/$protfolder/$violfiles
pdbnramafolder=$datafolder/$protfolder/$pdbfolder

emailfolderpath=$datafolder/$protfolder/$emailfolder
mkdir $emailfolderpath	
if [ -f $emailfolderpath/$emailfileslist ]; then
	rm $emailfolderpath/$emailfileslist
fi
#ls $pdbnramafolder/*pdb
pdbs_all_commasep_=`ls $pdbnramafolder/*pdb | tr '\n' ','`
pdbs_all_commasep=${pdbs_all_commasep_::-1}
ensemble_name=$pdbnramafolder"/"$protname"_ensemble.pdb"
$JOINPDB_CMD $pdbs_all_commasep $ensemble_name
if [ -f $ensemble_name ]; then
        echo "$ensemble_name file generated"
else
        echo "ERROR: $ensemble_name file not found"
fi

for pdbs in `ls $pdbnramafolder/*pdb`
do		
	file_pdb=`basename $pdbs`
	file=${file_pdb::-4}
	
	ramafile=$pdbnramafolder/$file_pdb".png"
	distviolfile=$violsfolder/$file".png"
	
	echo "$pdbs" >> $emailfolderpath/$emailfileslist
	if [ -f $ramafile ]; then
		if [ -f $distviolfile ]; then
			newfilename=$file".rama.viol.pdf"
			convert $ramafile $distviolfile $emailfolderpath/$newfilename
			echo "$emailfolderpath/$newfilename" >> $emailfolderpath/$emailfileslist
		else
			echo "$distviolfile NOT found"		
		fi
	else
		echo "$ramafile NOT found"
	fi
done

if [ -f $datafolder/$protfolder/$progressfile ]; then
        echo "$datafolder/$protfolder/$progressfile" >> $emailfolderpath/$emailfileslist
fi

if [ ! -z "$email" ]; then
	$SENDMAIL_CMD "$email" "$emailfolderpath/$emailfileslist" "$protname"
fi
# ------------- #
