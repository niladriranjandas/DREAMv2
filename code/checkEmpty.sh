#!/bin/bash

# Copyright (c) 2024 niladriranjandas
# check if input files in xml are empty or not

protname="$1"

alreadyrun="../proteinsparams/listoffiles.txt"
PROTFOLDER="../protein"

	found=`egrep ^$protname, $alreadyrun`
	if [ ! -z "$found" ]; then
        	paramfile_=`egrep ^$protname $alreadyrun | awk -F',' '{print $2}'`
        	paramfile=$PROTFOLDER/$paramfile_
	else
		paramfile=$PROTFOLDER"/"$protname"/"$protname".xml"		
	fi

#### check if files in empty ####

seq_file=`xmllint --xpath "string(//seq_file)" "$paramfile"`

upl_cyana_file=`xmllint --xpath "string(//upl_file)" "$paramfile"`
upl_xplor_file=`xmllint --xpath "string(//upl_file)" "$paramfile"`

aco_cyana_file=`xmllint --xpath "string(//ang_file)" "$xmlprotip"`
aco_xplor_file=`xmllint --xpath "string(//aco_xplor)" "$xmlprotip"`

        ##################### seq_file ############################
	if [ ! -z "$seq_file" ]; then
		if [ ! -f "$seq_file" ]; then
			echo "$seq_file NOT found"
			echo "Mail error"
			exit			
			if [ ! -s "$seq_file" ]; then
				echo "$seq_file is empty"
				echo "Mail error"
				exit				
			fi
		fi
	else
		echo "$seq_file not found"
		echo "Mail error"
		exit
	fi
	
        ################# upl_cyana_file #########################	
	if [ ! -z "$upl_cyana_file" ]; then
		if [ ! -f "$upl_cyana_file" ]; then
			echo "$upl_cyana_file NOT found"
			echo "Mail error"
			exit			
			if [ ! -s "$upl_cyana_file" ]; then
				echo "$upl_cyana_file is empty"
				echo "Mail error"
				exit				
			fi
		fi
	else
		echo "$upl_cyana_file not found"
		echo "Mail error"
		exit
	fi
	
        ##################### upl_xplor_file ############################
	if [ ! -z "$upl_xplor_file" ]; then
		if [ ! -f "$upl_xplor_file" ]; then
			echo "$upl_xplor_file NOT found"
			echo "Mail error"
			exit			
			if [ ! -s "$upl_xplor_file" ]; then
				echo "$upl_xplor_file is empty"
				echo "Mail error"
				exit				
			fi
		fi
	else
		echo "$upl_xplor_file not found"
		echo "Mail error"
		exit
	fi
	
        ##################### aco_cyana_file ############################
	if [ ! -z "$aco_cyana_file" ]; then
		if [ ! -f "$aco_cyana_file" ]; then
			echo "$aco_cyana_file NOT found"
			echo "Mail error"
			exit			
			if [ ! -s "$aco_cyana_file" ]; then
				echo "$aco_cyana_file is empty"
				echo "Mail error"
				exit				
			fi
		fi
	else
		echo "$aco_cyana_file not found"
		echo "Mail error"
		exit
	fi
	
        ##################### aco_xplor_file ############################
	if [ ! -z "$aco_xplor_file" ]; then
		if [ ! -f "$aco_xplor_file" ]; then
			echo "$aco_xplor_file NOT found"
			echo "Mail error"
			exit			
			if [ ! -s "$aco_xplor_file" ]; then
				echo "$aco_xplor_file is empty"
				echo "Mail error"
				exit				
			fi
		fi
	else
		echo "$aco_xplor_file not found"
		echo "Mail error"
		exit
	fi	
				

