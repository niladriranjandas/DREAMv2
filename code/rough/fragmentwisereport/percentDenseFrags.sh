#!/bin/bash

usage(){
	echo "./percentDenseFrags <reg_pdb_file> <full_pdb_file>"
}

chcknoModels(){
	pdbfile="$1"
	count=`egrep -c ^MODEL $pdbfile`
	if [ $count -gt 1 ]; then
		echo 0
	else
		echo 1
	fi
}

if [ $# -ne 2 ]; then
	echo "Error"
	usage
	exit 1
fi

if [ ! -f "$1" ]; then
	echo "Error: $1 file not present"
	exit 1
fi

if [ ! -f "$2" ]; then
	echo "Error: $2 file not present"
	exit 1
fi

reg_pdb="$1"
org_pdb="$2"

#chcknoModels $reg_pdb

if [ ! $(chcknoModels $reg_pdb) ]; then
	echo "Error: cannot handle multiple models"
	exit 1
fi
if [ ! $(chcknoModels $org_pdb) ]; then
	echo "Error: cannot handle multiple models"
	exit 1
fi

n_seq_reg=`egrep ^ATOM $reg_pdb | awk '{print $5}' | uniq | wc -l`   #no chain ID therefor $5 in awk
n_seq_org=`egrep ^ATOM $org_pdb | awk '{print $6}' | uniq | wc -l`

percent=`awk -v n_reg=$n_seq_reg -v n_org=$n_seq_org 'BEGIN{print n_reg*100/n_org}'`

echo $percent





