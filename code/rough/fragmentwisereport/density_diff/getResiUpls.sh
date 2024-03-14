#!/bin/bash

#----------------------------
# ./getResiUpls.sh <uplfile1,uplfile2,...> <resi1, resi2, ...>
#----------------------------

usage () {
	echo " ./getResiUpls <uplfile1,uplfile2,...> <resi1, resi2, ...>"
}

containsElement () {
	array_a="$1"
	srch_i="$2"

	for ele in "${array_a[@]}"
	do
		if [ "$srch_i" -eq "$ele" ]; then
			echo 1
			return
		fi
	done
	echo 0
}

upls="$1"
resi_list="$2"

if [ $# -ne 2 ]; then
	echo "ERROR "
	usage
fi

##########################################
strong=2.5
medium=5
weak=7

IFS="," read -r -a files <<< "$upls"
IFS="," read -r -a resi_list <<< "$resi_list"

#########################################

	for resis in "${resi_list[@]}"
	do
		strong_count=0
		medium_count=0
		weak_count=0
		for file_i in "${files[@]}"
		do
			resi_i=`cat $file_i | awk '{print $1}'`
			resi_j=`cat $file_i | awk '{print $4}'`
			dij_all=`cat $file_i | awk '{print $7}'`	

			#echo $resis $resi_i $resi_j $dij_all
			

			echo $resi_i | grep -q "$resis"
			flag_i=$?
			echo $resi_j | grep -q "$resis"				
			flag_j=$?
			
			flag_strong=0
			flag_medium=0
			flag_weak=0

			if [[ $flag_i -eq 0 || $flag_j -eq 0 ]]; then
				dij_all=`cat $file_i | awk '{print $7}'`
				echo $dij_all			
				IFS="" read -r -a dij_arr <<< "$dij_all"
				for dij in "${dij_arr[@]}"
				do				
					flag_strong=`awk -v hi="$strong" -v dij="$dij" 'BEGIN{if (dij <= hi) print 1}'`
					flag_medium=`awk -v hi="$medium" -v lo="$strong" -v dij="$dij" 'BEGIN{if (dij > lo && dij <= hi) print 1}'`
					flag_weak=`awk -v hi="$weak" -v lo="$meduim" -v dij="$dij" 'BEGIN{if (dij > lo && dij <= hi) print 1}'`

					if [[ $flag_strong -eq 1 ]]; then
						strong_count=$((strong_count+1))
					elif [[ $flag_medium -eq 1 ]]; then
						medium_count=$((medium_count+1))
					elif [[ $flag_weak -eq 1 ]]; then
						weak_count=$((weak_count+1))
					fi
				done
			fi
		done
		#printf "\n %s\t%d\t%d\t%d" "$resis" "$strong_count" "$medium_count" "$weak_count"
	done	
