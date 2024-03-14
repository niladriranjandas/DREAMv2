#!/bin/bash


#print the residue number from the pdb file

usage(){
	echo "./getResiFromPDB <pdb_file_name>"
}

if [[ $# -ne 1 ]]; then
	echo "error: not enough input"
	usage
	exit 1
fi

pdbfile="$1"

if [ ! -f $pdbfile ]; then
	echo "error: file not found $pdbfile"
	exit 1
fi	

startln=`echo $pdbfile | grep -in atom $pdbfile | head -1 | awk -F':' '{print $1}'`
endln=`echo $pdbfile | grep -in atom $pdbfile | tail -1 | awk -F':' '{print $1}'`

resi_tmp=`awk -v sl=$startln -v ndl=$endln 'FNR>=sl && FNR<=ndl' $pdbfile | awk '{print $5}' | uniq | awk 'BEGIN{ORS=","} {print $0}'`
resi="${resi_tmp::-1}"

echo $resi
