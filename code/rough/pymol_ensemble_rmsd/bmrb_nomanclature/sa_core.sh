#!/bin/bash

prots="2m4k,2kok,5x1x,1pbu,1xxe,2lav,2l7b"
dirs="2m4k_anchor_pseudo,2kok_anchor_pseudo,5x1x_anchor_pseudo,1pbu_anchor_pseudo,1xxe_anchor_pseudo,2lav_anchor_pseudo,2l7b_anchor_pseudo"

IFS="," read -r -a proteins <<< $prots
IFS="," read -r -a dir_arr <<< $dirs

target_prefix="/home/niladri/Downloads/niladri"

i=0
for dirs_i in "${dir_arr[@]}"
do
	pdb_name="grp_"${proteins[i]}".pdb"
	resis=`./getResiFromPDB.sh /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/"${proteins[i]}"/output_pdb/"$pdb_name"`
	for pdbs in `ls $target_prefix/$dirs_i/????_cluster_????.pdb`
	do
		name_=`echo $pdbs | awk -F'/' '{print $NF}'`
		name=`echo $name_ | awk -F'.pdb' '{print $1}'`
		name_mod=$name"_core".pdb
		target_file=$target_prefix"/"$dirs_i"/"$name_mod
		if [ ! -z $resis ]; then
			./writePDBforRange.sh $pdbs $resis > $target_file
		fi		
	done
	i=$((i+1))
done
