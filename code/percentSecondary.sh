#!/bin/bash

usage(){
	echo "usage : ./percentSeconday.sh <pdb_file_name>"
}

givSecndryFrm1letter(){
	#from DSSP/dssp-3.0.0/doc/$DSSP.1
	#H;Alpha Helix
	#B;Beta Bridge
	#E;Strand
	#G;Helix\-3
	#I;Helix\-5
	#T;Turn
	#S;Bend
	
	abbr="#H;Alpha_Helix #B;Beta_Bridge #E;Strand #G;Helix\-3 #I;Helix\-5 #T;Turn #S;Bend"

	alphabet="$1"
	
	secondary_strct=`echo $abbr | tr ' ' '\n' | egrep ^#$alphabet | awk -F';' '{print $2}'`
	echo $secondary_strct
}

if [ $# -ne 1 ]; then
	echo "ERROR"
	usage
	exit 1
fi

pdb_file="$1"

if [ ! -f $pdb_file ]; then
	echo "ERROR: $pdb_file doesnot exist."
	exit 1
fi

name_file=`echo $pdb_file | awk -F'/' '{print $NF}'`
protein_name=`echo $name_file | awk -F'.pdb' '{print $1}'`
 
### run $DSSP to get secondary structures ####
	dssp_file=$protein_name"_dssp.txt"
	$DSSP -i $pdb_file |  awk '/#  RESIDUE/,/^$/' | awk '{print $1 " " substr($0,17,1)}' | sed -n '1!p' > $dssp_file

	if [ $? -ne 0 ]; then
		echo "ERROR: dssp run error"
		exit 1
	fi

### calc the percent ####
	if [ ! -f $dssp_file ]; then
		echo "ERROR: dssp file not found"
		exit 1
	fi
	
	n_resi=`tail -1 $dssp_file | awk '{print $1}'`

	printf "%s\t%s" "$name_file" $n_resi
	printf "\n___________________________________"	

	for secndry in `cat $dssp_file | awk '{print $2}' | sed '/^$/d' | sort -u`
	do
		secndry_strct=`givSecndryFrm1letter $secndry`
		

		count=`grep -c $secndry $dssp_file`
		percnt=`awk -v n_resi=$n_resi -v count=$count 'BEGIN{print count*100/n_resi}'`
		#printf "\n%4s\t%15s\t%4s\t%6s" $secndry "$secndry_strct" $count $percnt
		printf "\n%-15s\t%4s\t%6s" "$secndry_strct" $count $percnt

		count=""
		percnt=""
	done

echo
