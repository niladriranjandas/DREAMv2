#!/bin/bash
#
#Copyright (c) 2024 niladriranjandas
#convert the dihedral angle bound from CYANA to XPLOR format for "normal mode"

protname="$1"

aco_convert_cmd="python3 ../packages/PdbStat510_20130909/acoCyanaToXplor.py"

alreadyrun="../proteinsparams/listoffiles.txt"

 PROTFOLDER="../protein"
 xmlprotip=$PROTFOLDER"/"$protname"/"$protname".xml"

 aco_cyana_file=`xmllint --xpath "string(//ang_file)" "$xmlprotip"`
 aco_xplor_file=`xmllint --xpath "string(//aco_xplor)" "$xmlprotip"`

 if [ ! -f $aco_cyana_file ]; then
		#echo " $aco_convert_cmd $PROTFOLDER/$protname/$aco_cyana_file > $PROTFOLDER/$protname/$aco_xplor_file"
 		$aco_convert_cmd "$PROTFOLDER/$protname/$aco_cyana_file" > "$PROTFOLDER/$protname/$aco_xplor_file"
 fi
