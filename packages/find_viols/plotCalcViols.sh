#!/bin/bash

protname="$1"

alreadyrun="../proteinsparams/listoffiles.txt"
datafolder="../protein"
violfiles="distviols"

#PYTHON_CODE="../packages/find_viols/findViolations_s_m_w.py"
#PYTHON_CODE="../packages/find_viols/findViolations_s_m_w_v2.py"
PYTHON_CODE="../packages/find_viols/findViolations_s_m_w_v3.py"

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
        protfolder=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`
	violsfolder=$datafolder/$protfolder/$violfiles
	
	for viols in `ls $violsfolder/*onlyviols`
	do
		filename=`basename $viols`
		filename_noext=${filename::-10}
		#echo "python3 $PYTHON_CODE $viols $filename_noext"
		python3 $PYTHON_CODE $viols $filename_noext
	done 
else
	protfolder=$protname"_allparam"
	violsfolder=$datafolder/$protfolder/$violfiles
	
	for viols in `ls $violsfolder/*onlyviols`
	do
		filename=`basename $viols`
		filename_noext=${filename::-10}
		#echo "python3 $viols $filename_noext"
		python3 $PYTHON_CODE $viols $filename_noext
	done 	
fi
