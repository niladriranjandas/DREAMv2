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
strong=2.7
medium=3.5
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
			while IFS='' read -r line || [[ -n "$line" ]]; do
			    resi_i=`echo "$line" | awk '{print $1}'`
			    resi_j=`echo "$line" | awk '{print $4}'`
			    dij=`echo "$line" | awk '{print $7}'`	
			    atm_i=`echo "$line" | awk '{print $3}'`
			    atm_j=`echo "$line" | awk '{print $6}'`
			    if [[ "$resi_i" -eq "$resis" || "$resi_j" -eq "$resis" ]]; then
				flag_strong=0
				flag_medium=0
				flag_weak=0

				flag_strong=`awk -v hi="$strong" -v dij="$dij" 'BEGIN{if (dij > 0 && dij <= hi) print 1}'`
				flag_medium=`awk -v hi="$medium" -v lo="$strong" -v dij="$dij" 'BEGIN{if (dij > lo && dij <= hi) print 1}'`
				flag_weak=`awk -v hi="$weak" -v lo="$meduim" -v dij="$dij" 'BEGIN{if (dij > lo && dij <= hi) print 1}'`

                                flag_gt_1=`awk -v resii="$resi_i" -v resij="$resi_j" 'BEGIN{if ((resii-resij)>1 || (resij-resii)>1) print 1; else print 0;}'`
				if [[ "$flag_gt_1" -eq 1 ]]; then
					if [[ "$flag_strong" -eq 1 ]]; then
						#printf "\n %s\t%s\t%s\t%s\t%f" "$resi_i" "$atm_i" "$resi_j" "$atm_j" "$dij"
						#printf "\n %s,%s,%s,%s,%f" "$resi_i" "$atm_i" "$resi_j" "$atm_j" "$dij"
						strong_count=$((strong_count+1))
					elif [[ "$flag_medium" -eq 1 ]]; then
						#printf "\n %s\t%s\t%s\t%s\t%f" "$resi_i" "$atm_i" "$resi_j" "$atm_j" "$dij"
						printf "\n %s,%s,%s,%s,%f" "$resi_i" "$atm_i" "$resi_j" "$atm_j" "$dij"
						medium_count=$((medium_count+1))
					elif [[ "$flag_weak" -eq 1 ]]; then
						weak_count=$((weak_count+1))
					fi
				fi
			    fi
			done < "$file_i"
		done
		#printf "\n %s\t%d\t%d\t%d" "$resis" "$strong_count" "$medium_count" "$weak_count"
	done
echo

