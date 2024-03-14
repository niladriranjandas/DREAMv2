#!/bin/bash

usage(){
	echo "convertBMRBatom <resi_3_letter> <atom_name>"
}

threeToOne(){
	resi3="$1"
	aminoacids="$2"
	grep -i $resi3 $aminoacids| awk '{print $2}'
}


nomanclature_file="nomanclature_edit.csv" 		# for Sa and ancloc_refine.pdb
#nomanclature_file="nomanclature_edit_cyanacol5.csv"	# for CYANA
aminoacids="aminoacid.txt"
resi3="$1"
atomname="$2"

if [ $# -ne 2 ]; then
	usage
	exit 1
fi

if [ ! -f $nomanclature_file ]; then
	echo "nomanclature file not found"
	exit 1
fi

resi1=`threeToOne $resi3 $aminoacids`
if [ -z $resi1 ]; then
	echo "$resi3 not found"
	exit 1
fi
#echo $zzz
#resi1=`treeToOne $resi3`
#if [ $resi1 == "X" ]; then
#	echo "atom name could not be found"
#	exit 1
#fi
#resi1="R"

lines=`grep "^$resi1" $nomanclature_file | awk '{print $2}' | grep -n  -w "^$atomname" | awk -F':' '{print $1}'`
if [ ! -z $lines ]; then
	lines_=$((lines-1))

	pdb=`grep "^$resi1" $nomanclature_file | awk '{print $4}'`      # fpr our algo after SA
	#pdb=`grep "^$resi1" $nomanclature_file | awk '{print $2}'`     # for ancloc_refine.pdb
	#pdb=`grep "^$resi1" $nomanclature_file | awk '{print $5}'`     # for CYANA

	IFS=$' ' read -d '' -r -a pdb_atoms <<< $pdb
	pdb_atom_name=${pdb_atoms[$lines_]}

	echo $pdb_atom_name
else
	echo XX
	#echo $atomname
fi

