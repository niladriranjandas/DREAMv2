#!/bin/bash

function getpseudo(){
	resi_nm="$1"
	atm_nm="$2"
	filename="bmrb_nomanclature/pseudo.csv"

	atoms_=`grep -i $resi_nm[[:space:]]*$atm_nm $filename`
	atom=`echo $atoms_ | awk '{print $NF}'`

	ret_atm=""
	IFS="," read -r -a atom_arr <<< "$atom"
	for atom_i in "${atom_arr[@]}"
	do
		ret_atm_=`./getGROfromBMRB_edit.sh $resi_nm $atom_i`
		ret_atm="$ret_atm","$ret_atm_"
	done
	echo ${ret_atm:1}
}

gro_file="/root/Downloads/gmx_22nov17/lisa_10_64_25_37/2m4k/script1/2m4k_1_ancloc_refine_npt.gro"
upl_file="/home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/2m4k/2m4k.dist.1.upl"
pseudo_file="/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/pymol_ensemble_rmsd/bmrb_nomanclature/pseudo.csv"

write_line_file="lines_for_posre.txt"
posre_file_name="posre_H.itp"

if [[ -f line.txt ]]; then
	echo "removing temp file line.txt"
	sudo rm line.txt
fi

while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "========================================================================"		
	i_resi_no=`echo $line | awk '{print $1}'`
	i_resi_nm=`echo $line | awk '{print $2}'`	
	i_atom_nm_=`echo $line | awk '{print $3}'`
	echo $i_atom_nm_
	if [[ "$i_atom_nm_" == *Q* ]]; then
		echo "here"
		#getpseudo $i_resi_nm $i_atom_nm_
		i_atom_nm=`getpseudo $i_resi_nm $i_atom_nm_`
		#i_atom_nm="XX"
	else
		i_atom_nm=`./getGROfromBMRB_edit.sh $i_resi_nm $i_atom_nm_`
	fi

	j_resi_no=`echo $line | awk '{print $4}'`
	j_resi_nm=`echo $line | awk '{print $5}'`	
	j_atom_nm_=`echo $line | awk '{print $6}'`
	echo $j_atom_nm_
	if [[ "$j_atom_nm_" == *Q* ]]; then
		echo "here"
		#getpseudo $i_resi_nm $i_atom_nm_
		j_atom_nm=`getpseudo $j_resi_nm $j_atom_nm_`
		#j_atom_nm="XX"
	else
		j_atom_nm=`./getGROfromBMRB_edit.sh $j_resi_nm $j_atom_nm_`
	fi

		
	IFS="," read -r -a atoms_of_i <<< $i_atom_nm	
	for atoms_in_i in "${atoms_of_i[@]}"
	do
		#srch_pat_i=^[[:space:]]*"$i_resi_no$i_resi_nm"[[:space:]]*"$i_atom_nm"[[:space:]]
		srch_pat_i=^[[:space:]]*"$i_resi_no$i_resi_nm"[[:space:]]*"$atoms_in_i"[[:space:]]
		gro_line_i=`grep $srch_pat_i $gro_file`
		echo "grep $srch_pat_i $gro_file"
		echo "$gro_line_i" 
		echo "$gro_line_i" | awk '{print $3}' >> line.txt
	done

	IFS="," read -r -a atoms_of_j <<< $j_atom_nm		
	for atoms_in_j in "${atoms_of_j[@]}"
	do
		#srch_pat_j=^[[:space:]]*"$j_resi_no$j_resi_nm"[[:space:]]*"$j_atom_nm"[[:space:]]
		srch_pat_j=^[[:space:]]*"$j_resi_no$j_resi_nm"[[:space:]]*"$atoms_in_j"[[:space:]]
		gro_line_j=`grep $srch_pat_j $gro_file`
		echo "grep $srch_pat_j $gro_file"
		echo "$gro_line_j" 
		echo "$gro_line_j" | awk '{print $3}' >> line.txt
	done
done < "$upl_file"

cat line.txt | sort -n | uniq > $write_line_file

line_1="[ position_restraints ]"
line_2="; atom  type      fx      fy      fz"

echo $line_1 > $posre_file_name
echo $line_2 >> $posre_file_name

while IFS='' read -r line_atm || [[ -n "$line_atm" ]]; do
	if [[ ! -z $line_atm ]]; then
		printf "\n%6s%6s  %s  %s  %s" $line_atm 1 1000 1000 1000 >> $posre_file_name
	fi
done < "$write_line_file"


