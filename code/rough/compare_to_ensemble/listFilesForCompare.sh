#!/bin/bash

#e.g: ./listFilesForCompare.sh ~/Documents/Disco_etc_all_in_1/our_algo/archive/ ~/Documents/Disco_etc_all_in_1/SPROS_try/proteins/

if [ $# != 2 ]; then
	echo "usage: ./listFilesForCompare <path_to_gen_file_archive> <path_to_ground_truth dir>"
	exit 1
fi

if [ ! -d $1 ]; then
	echo "err: generated file archive dir not found."
	exit 1
fi

if [ ! -d $2 ]; then
	echo "err: ground truth dir not found."
fi

##########################################################################################################
pattern_gen="nosol.pdb"

##########################################################################################################
 curr_dir=`pwd`
	for i in `ls ${1}`
	do
		#----locate generated file----#
		gen_file_=`find $1 -iname "*$i*$pattern_gen"`
                for j in `echo $gen_file_`
                do
			if [ ! -z $j ]; then
				gen_file_name_=`echo ${j} | awk -F"/" '{print $NF}'`
				gen_file_name=`echo ${gen_file_name_} | awk -F".pdb" '{print $1}'`
				gen_file_path=`echo ${j} | awk -F"/$gen_file_name_" '{print $1}'`
				printf "%s,%s," $gen_file_name $gen_file_path

				#--------locate the ground truth -------------#
				truth_file_=`find $2 -iname "${i}.pdb"`
				if [ ! -z $truth_file_ ]; then
					truth_file_name_=`echo ${truth_file_} | awk -F"/" '{print $NF}'`
					truth_file_name=`echo ${truth_file_name_} | awk -F".pdb" '{print $1}'`
					truth_file_path=`echo ${truth_file_} | awk -F"/$truth_file_name_" '{print $1}'`
					printf "%s,%s\n" $truth_file_name $truth_file_path
				else
					printf "FILE NOT FOUND\n"
				fi
			else
				printf "%s,FILE NOT FOUND\n" $i
			fi
    	        done
	done


				
