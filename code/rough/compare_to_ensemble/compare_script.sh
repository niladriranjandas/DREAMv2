#!/bin/bash

#usage: /compare_script.sh compare_files_all.txt all_log.txt best_rmsd_log.txt

if [ $# != 3 ]; then
   echo "usage: ./compare_script.sh <compare_file> <detailed_log_file> <best_rmsd_log_file>"
   exit 1
fi

if [ ! -f $1 ]; then
   echo "compare file not found $1"
   exit 1
fi

#--------- read the compare file -----------
	while read -r line
        do
           gen_name=`echo $line | awk -F',' '{print $1}'`
           gen_file_path=`echo $line | awk -F',' '{print $2}'`
           gen_file_full=$gen_file_path/$gen_name'.pdb'

           truth_name=`echo $line | awk -F',' '{print $3}'`
           truth_file_path=`echo $line | awk -F',' '{print $4}'`
           truth_file_full=$truth_file_path/$truth_name'.pdb'

           if [ ! -f $gen_file_full ]; then
                echo "${gen_file_full} not found"
           elif [ ! -f $truth_file_full ]; then
                echo "${truth_file_full} not found"
           else
                pymol -cq compare_truth.py -- $gen_name $gen_file_path $truth_name $truth_file_path >> $2
               #echo "pymol -cq compare_truth.py -- $gen_name $gen_file_path $truth_name $truth_file_path "
                if [ $? -eq 0 ]; then
                	echo "${gen_name} : ${truth_name} success"
                        pymol -cq compare_truth_printbest.py -- $gen_name $gen_file_path $truth_name $truth_file_path >> $3
			if [ $? != 0 ]; then
				echo "-------best rmsd file write error------"
			fi
                else
                        echo "${gen_name} : ${truth_name} error"
                fi
           fi
        done < $1
                

#-------------------
grep -v -i '^pymol' $3 > tmp.txt
mv tmp.txt $3
