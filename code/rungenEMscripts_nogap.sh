#!/bin/bash

#max=5
max=1

pdb_name="$1"
pdb_folder="$2"
upl_file="$3"

#pattern=$pdb_name"_[0-9]_ancloc_refine_[0-9]+.pdb"
#files_=`ls $pdb_folder | egrep $pattern | tr '\n' ','`
#files=${files_::-1}

c=1
files=""
while [[ $c -le $max ]]
do
	#pattern=$pdb_name"_"$c"_ancloc_refine_[0-9]+.pdb"
	pattern="grp_registration_"$pdb_name".pdb"
        #echo $pattern
	file_=`ls $pdb_folder | egrep $pattern | tr '\n' '|' | awk -F'|' '{print $(NF-1)}'`
	echo $file_
	if [ ! -z $file_ ]; then
		files_tmp=$files","$pdb_folder"/"$file_        #files_tmp=$files","$file_ 
        	files=$files_tmp
	fi
        c=$((c+1))
done

final_files=${files:1}

echo "./genEMscripts.sh $final_files $upl_file $pdb_name"
