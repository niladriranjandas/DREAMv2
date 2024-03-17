#!/bin/bash

#
#Copyright (c) 2024 niladriranjandas
#send emails to end user

protname="$1"
email="$2"

SENDMAIL_CMD="python front_end/sendMail.py"


if [ -z $protname ]; then
	echo "arg-1 needed"
	exit -1
fi

alreadyrun="../proteinsparams/listoffiles.txt"
protpath="../protein"
foldername="finalpdbfiles"

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
	protfolder=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`
        
	waterref_folder=$protpath"/"$protfolder"/"$foldername
	echo $waterref_folder
	oppfilename=$waterref_folder"/"$protname"_pdbfiles.txt"
	if [ -f "$oppfilename" ]; then
		rm -f "$oppfilename"
	fi
	for pdbs in `ls $waterref_folder/*pdb`
	do
		curr_path=`pwd`
		full_path=$curr_path"/"$pdbs
		
		echo "$full_path" >> "$oppfilename"
	done
	
	for pdbs in `ls $waterref_folder/*png`
	do
		curr_path=`pwd`
		full_path=$curr_path"/"$pdbs
		
		echo "$full_path" >> "$oppfilename"
	done	
        echo "$SENDMAIL_CMD" "$email" "$oppfilename" "$protname"	
	$SENDMAIL_CMD "$email" "$oppfilename" "$protname" 
else
	#echo "code it"
	protfolder=$protname"_allparam"
        
	waterref_folder=$protpath"/"$protfolder"/"$foldername
	echo $waterref_folder
	oppfilename=$waterref_folder"/"$protname"_pdbfiles.txt"
	if [ -f "$oppfilename" ]; then
		rm -f "$oppfilename"
	fi
	for pdbs in `ls $waterref_folder/*pdb`
	do
		curr_path=`pwd`
		full_path=$curr_path"/"$pdbs
		
		echo "$full_path" >> "$oppfilename"
	done
	
	for pdbs in `ls $waterref_folder/*png`
	do
		curr_path=`pwd`
		full_path=$curr_path"/"$pdbs
		
		echo "$full_path" >> "$oppfilename"
	done	
        echo "$SENDMAIL_CMD" "$email" "$oppfilename" "$protname"	
	$SENDMAIL_CMD "$email" "$oppfilename" "$protname" 

fi
