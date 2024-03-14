#!/bin/bash

pdb_file="$1"
aco_file="$2"

./pdbtorsions $pdb_file > tmp_torsion.txt

diff=0
while IFS='' read -r line || [[ -n "$line" ]]; do
	resi_no=`echo $line | awk '{print $1}'`
	resi_nm=`echo $line | awk '{print $2}'`
	resi_nm_cap=`echo $resi_nm | tr '[:lower:]' '[:upper:]'` 
	phi_psi=`echo $line | awk '{print $3}'`
	lo=`echo $line | awk '{print $4}'`
	hi=`echo $line | awk '{print $5}'`
	
	
	srch_str=$resi_nm_cap$resi_no
	line=`grep $srch_str tmp_torsion.txt`
	phi=`echo $line | awk '{print $2}'`
	psi=`echo $line | awk '{print $3}'`

	if [ $phi_psi == "PHI" ]; then
		if [ 1 -eq "$(echo "${phi} < ${lo}" | bc)" ]; then
			diff=$(echo "$lo - $phi"|bc)
		elif [ 1 -eq "$(echo "${phi} > ${hi}" | bc)" ]; then
			diff=$(echo "$phi -$hi"|bc)
		fi
	elif [ $phi_psi == "PSI" ]; then
		if [ 1 -eq "$(echo "${psi} < ${lo}" | bc)" ]; then
			diff=$(($lo - $psi))
		elif [ 1 -eq "$(echo "${psi} > ${hi}" | bc)" ]; then
			diff=$(($psi - $hi))
		fi
	else
		display "error"
	fi

done < $aco_file
