#!/bin/bash

#from DSSP/dssp-3.0.0/doc/mkdssp.1
#H;Alpha Helix
#B;Beta Bridge
#E;Strand
#G;Helix\-3
#I;Helix\-5
#T;Turn
#S;Bend

usage(){
        echo "usage : ./percentSeconday.sh <pdb_file_name>"
}


if [ $# -ne 1 ]; then
        echo "ERROR"
        usage
        exit 1
fi

#ignore='*,S,T,G'
#include='H,B,E,G,I'       # changed on 3-3-21
include='H,B,E,G,I,T,S'    # changed on 3-3-21
include='H,B,E,G,I,T'       # changed on 7-5-21
pdb_file="$1"

if [ ! -f $pdb_file ]; then
        echo "ERROR: $pdb_file doesnot exist."
        exit 1
fi

name_file=`echo $pdb_file | awk -F'/' '{print $NF}'`
protein_name=`echo $name_file | awk -F'.pdb' '{print $1}'`

### run mkdssp to get secondary structures ####
        dssp_file=$protein_name"_dssp.txt"
        #mkdssp -i $pdb_file |  awk '/#  RESIDUE/,/^$/' | awk '{print $1 " " substr($0,17,1)}' | sed -n '1!p' | sed 's/  / */g'> $dssp_file   # changed on 3-3-21
        #mkdssp -i $pdb_file |  awk '/#  RESIDUE/,/^$/' | awk '{print $2 " " substr($0,17,1)}' | sed -n '1!p' | sed 's/  / */g'> $dssp_file  # changed on 3-3-21
	dssp -i $pdb_file |  awk '/#  RESIDUE/,/^$/' | awk '{print $2 " " substr($0,17,1)}' | sed -n '1!p' | sed 's/  / */g'> $dssp_file  # changed on 4-1-23

        if [ $? -ne 0 ]; then
                echo "ERROR: dssp run error"
                exit 1
        fi

###############################################
# which to ignore. blank is already considered
#IFS=', ' read -r -a array <<< "$ignore"

# which to ignore. 
IFS=', ' read -r -a array <<< "$include"

while IFS= read -r line || [[ -n "$line" ]]; do
    if [ ! -z "$line" ]; then
         sec=`echo $line | awk '{print $2}'`
         resi=`echo $line | awk '{print $1}'`        
    	 if [[ " ${array[@]} " =~ " ${sec} " ]]; then
	    printf "%s," "$resi"
	 fi
    fi
done < "$dssp_file"

