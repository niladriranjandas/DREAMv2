#!/bin/bash

# code to get the pdb for certain residues supplied in a list separated by comma 
#

usage(){
	echo "usage: ./writePDBforRange.sh <pdb-file> <list separated by comma>"
	echo "       e.g., ./writePDBforRange.sh grp_1xxe.pdb 2,3,4,5"
}

if [[ $# != 2 ]]; then
	echo "error"
	usage
	exit 1
fi

pdbfile="$1"
resilist=",$2,"

if [ ! -f $pdbfile ]; then
	echo "No pdb file found in $pdbfile"
	exit 1
fi

if [ -z $resilist ]; then
	echo "residue list cannot be empty: $resilist"
fi


	while IFS='' read -r line || [[ -n "$line" ]]; 
	do
		echo $line | egrep -qi ^atom
		if [[ $? -eq 0 ]]; then
			resi_i=`echo "$line" | awk '{print $5}'`
			if [ -z $resi_i ]; then
				echo "$line"
				echo "error: residue not found in pdb file"
				exit 1
			else
				echo $resilist | grep -q ",$resi_i,"
				if [[ $? -eq 0 ]]; then
					echo "$line"
				fi
			fi
		else
			echo "$line"
		fi
	done < $pdbfile


		
