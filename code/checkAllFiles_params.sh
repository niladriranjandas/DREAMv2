#!/bin/bash

# Copyright (c) 2024 niladriranjandas
# check if all the files in xml are present


isFileEmpty () {
	filename="$1"
	
	if [ ! -f $filename ]; then
		return 0
	else
		if [ ! -s $filename ]; then
			return 0
		else
			return 1
		fi
	fi
}

protname="$1"
alreadyrun="../proteinsparams/listoffiles.txt"

PROTFOLDER="../protein"
xmlprotip=$PROTFOLDER"/"$protname"/"$protname".xml"
ERR_FLAG=1

		seq_cyana_file=`xmllint --xpath "string(//seq_file)" "$xmlprotip"`

		upl_cyana_file=`xmllint --xpath "string(//upl_file)" "$xmlprotip"`
		upl_xplor_file=`xmllint --xpath "string(//upl_xplor)" "$xmlprotip"`

		hbond_cyana_file=`xmllint --xpath "string(//upl_xplor)" "$xmlprotip"`

		aco_cyana_file=`xmllint --xpath "string(//ang_file)" "$xmlprotip"`
		aco_xplor_file=`xmllint --xpath "string(//aco_xplor)" "$xmlprotip"`

		 if [ ! -z $seq_cyana_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$seq_cyana_file"
        	        if [ "$?" == 0 ]; then
        	        	echo "ERROR: $seq_cyana_file NOT found"
        	        	ERR_FLAG=0
        	        	#exit
        	        fi
		 fi
		 if [ ! -z $upl_cyana_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$upl_cyana_file"
        	        if [ "$?" == 0 ]; then
	        	        echo "ERROR: $upl_cyana_file NOT found"
	        	        ERR_FLAG=0
        	        	#exit
        	        fi
		 fi
		 if [ ! -z $upl_xplor_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$upl_xplor_file"
        	        if [ "$?" == 0 ]; then
	        	        echo "ERROR: $upl_xplor_file NOT found"
	        	        ERR_FLAG=0
        	        	#exit
        	        fi
		 fi
		 if [ ! -z $hbond_cyana_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$hbond_cyana_file"
        	        if [ "$?" == 0 ]; then
	        	        echo "ERROR: $hbond_cyana_file NOT found"
	        	        ERR_FLAG=0
        	        	#exit
        	        fi
		 fi	
		 if [ ! -z $aco_cyana_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$aco_cyana_file"
        	        if [ "$?" == 0 ]; then
	        	        echo "ERROR: $aco_cyana_file NOT found"
	        	        ERR_FLAG=0
        	        	#exit
        	        fi
		 fi		 	 		 		 
		 if [ ! -z $aco_xplor_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$aco_xplor_file"
        	        if [ "$?" == 0 ]; then
	        	        echo "ERROR: $aco_xplor_file NOT found"
	        	        ERR_FLAG=0
        	        	#exit
        	        fi
		 fi
		  
		 if [ $ERR_FLAG == 0 ]; then
		 	echo "ERROR: exiting due to missing/empty file"
		 	exit
		 fi
		 
		 echo "--BEGIN--"
