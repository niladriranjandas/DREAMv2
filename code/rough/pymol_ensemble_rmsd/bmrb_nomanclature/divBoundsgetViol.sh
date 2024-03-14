#!/bin/bash
# ./divBoundsgetViol.sh <protein_ensemble> <upl=upl_file> <lo=lo_file>
#   protein_ensemble must be present in the current folder

usage(){
	echo "usage: /divBoundsgetViol.sh <protein_ensemble> <upl_file> [<lo_file>]"
}

if [[ "$#" -lt 2 ]]; then
	echo $#
	usage
	exit 1
fi

protein_ensemble="$1"
if [[ ! -f $protein_ensemble ]]; then
	echo "$protein_ensemble file not present"
	exit 1
fi
protein_ensemble_name=`echo $protein_ensemble | awk -F'.pdb' '{print $1}'`
if [[ $# == 2 ]]; then
	upl="$2"
	if [[ ! -f $upl ]]; then
		echo "$upl file not present"
		exit 1
	fi
fi

if [[ $# == 3 ]]; then
	lo="$3"
	if [[ ! -f $lo ]]; then
		echo "$lo file not present"
		exit 1
	fi
fi

## set bounds for strong/medium/weak bounds
l1=0		#strong
u1=2.7

l2=2.7		#medium
u2=3.5

l3=3.5		#weak
u3=6

## create a temporary folder to save the ensemble
if [[ ! -d divBoundsgetViol_tmp ]]; then
	mkdir divBoundsgetViol_tmp
fi
cp $protein_ensemble divBoundsgetViol_tmp/.

pymol -cq callSplitnSave.py -- $protein_ensemble_name divBoundsgetViol_tmp
if [[ $? -ne 0 ]]; then
	echo "Error in callSplitnSave.py"
	exit 1
fi

pat=$protein_ensemble_name"_????.pdb"

files_=`ls divBoundsgetViol_tmp/$pat`
if [[ -z $files_ ]]; then
	echo "No files found after spliting"
	#exit 1
	pat=$protein_ensemble_name".pdb"
	files_=`ls divBoundsgetViol_tmp/$pat`
	echo "Trying to proceed with PDB input"
	if [[ -z $files_ ]]; then
		echo "No files found with/without splitting"
		exit 1
	fi
fi
cd divBoundsgetViol_tmp
files_nodir=`ls $pat`
cd ..
files_tmp="${files_nodir//.pdb/,}"
files_tmp2="${files_tmp//[[:space:]]/}"
files="${files_tmp2::-1}"


## split the upls into categories
./writeUplForDistRange.sh "$upl" $l1,$u1 > divBoundsgetViol_tmp/$protein_ensemble_name"_strong.upl"
./writeUplForDistRange.sh "$upl" $l2,$u2 > divBoundsgetViol_tmp/$protein_ensemble_name"_medium.upl"
./writeUplForDistRange.sh "$upl" $l3,$u3 > divBoundsgetViol_tmp/$protein_ensemble_name"_lo.upl"

## call the getViol module

#pymol -cq getViolations.py -- "/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/gromacs_pos_restraint/2m4k_trr_extract/script_1" 2m4k_1_clusters_0001 /home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/pymol_ensemble_rmsd/2m4k.dist.1.upl_cpy

echo "------------distance viol: strong--------------"
pymol -cq getViolations.py -- divBoundsgetViol_tmp $files divBoundsgetViol_tmp/$protein_ensemble_name"_strong.upl"
echo "pymol -cq getViolations.py -- divBoundsgetViol_tmp $files divBoundsgetViol_tmp/$protein_ensemble_name _strong.upl"
echo "------------distance viol: medium--------------"
pymol -cq getViolations.py -- divBoundsgetViol_tmp $files divBoundsgetViol_tmp/$protein_ensemble_name"_medium.upl"
echo "------------distance viol: weak----------------"
pymol -cq getViolations.py -- divBoundsgetViol_tmp $files divBoundsgetViol_tmp/$protein_ensemble_name"_lo.upl"

## clean up
rm  $files_
rm divBoundsgetViol_tmp/$protein_ensemble_name".pdb"

