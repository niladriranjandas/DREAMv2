#!/bin/bash

#*** write the upl's corr to a range of residues****
# usage: ./writeUplForRange.sh <upl_file> <range=lo1,up1:lo2,upl2>
##########################################

usage(){
	echo "usage: /writeUplForRange.sh <upl_file> <lo1,up1:lo2,upl2>"
}

if [[ $# -ne 2 ]]; then
	usage
	exit 1
fi

upl_file="$1"
range="$2"
if [ ! -f $upl_file ]; then
	echo "upl file $upl_file not found"
	exit 1
fi 

######################################################
	IFS=":" read -r -a range_arr <<< $range
	
	while IFS='' read -r line || [[ -n "$line" ]]; do
		resi_i=`echo $line | awk '{print $1}'`
		resi_j=`echo $line | awk '{print $4}'`

		for range_i in "${range_arr[@]}"
		do
			lo1=`echo $range_i | awk -F',' '{print $1}'`
			up1=`echo $range_i | awk -F',' '{print $2}'`

			if [[ $lo1 -le $up1 ]]; then
				if [[ $resi_i -ge $lo1 && $resi_i -le $up1 ]]; then
					echo "$line"
				elif [[ $resi_j -ge $lo1 && $resi_j -le $up1 ]]; then
					echo "$line"
				fi
			fi
		done		
	done < $upl_file
#####################################################


