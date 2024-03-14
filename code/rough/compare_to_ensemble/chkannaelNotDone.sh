#!/bin/bash

#check the archive/.. protein names, search Experiment_list/NMR.csv and fill no. residues
#check archive/????/*nosol.pdb and fill names if present
# usage: ./chkannaelNotDone /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive /home/niladri/Documents/Disco_etc_all_in_1/our_algo/Experiment_list/NMR.csv
#         ./chkannaelNotDone

if [ -z $1 ]; then
	archive_dir='/home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive'
else
	archive_dir=$1	
fi

if [ -z $2 ]; then
	exp_file='/home/niladri/Documents/Disco_etc_all_in_1/our_algo/Experiment_list/NMR.csv'
else
	exp_file=$2
fi

pattern='nosol.pdb$'
#--------------------------------------------------------------------------------
	for i in `ls -d $archive_dir/????`
	do
            protein_name=`echo $i | awk -F'/' '{print $NF}'`
               #get the resi lenght from downloaded file
		resi_len_=`egrep -i $protein_name $exp_file`
		if [ ! -z "$resi_len_" ]; then
			resi_len=`echo $resi_len_ | awk -F',' '{print $2}'`
		else
			resi_len="NOT FOUND"
		fi

               #get the details about generated file
		protein_lowercase=`echo $protein_name | tr '[:upper:]' '[:lower:]'`
		gen_file=`ls $archive_dir/$protein_lowercase/gromacs_stuff | egrep $pattern`
		if [ -z "$gen_file" ]; then
			gen_file="NOT FOUND"
		fi

             printf "%s,%s,%s\n" $protein_name "$resi_len" "$gen_file"  
	done


