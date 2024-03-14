#!/bin/bash

usage{}{
	echo "usage: ./sa_md_template.sh 2m4k_1"
	exit 1
}

checkErr(){
	if [ "$1" -ne 0 ]; then
		echo "Error in $2"
		exit 1
	fi
}

protein_name="2m4k_1"
posre_filename="posre.itp"

gmx pdb2gmx -f "$protein_name"_ancloc_refine.pdb -water tip3p -o "$protein_name"_process.gro -ignh <<<8
checkErr "$?" "pdb2gmx" 
gmx editconf -f "$protein_name"_process.gro -o "$protein_name"_newbox.gro -c -d 1.0 -bt cubic
checkErr "$?" "editconf"

gmx solvate -cp "$protein_name"_newbox.gro -cs spc216.gro -o "$protein_name"_solv.gro -p topol.top 
checkErr "$?" "solvate"
gmx grompp -f ions.mdp -c "$protein_name"_solv.gro -p topol.top -o ions.tpr
checkErr "$?" "grompp_ions"
gmx genion -s ions.tpr -o "$protein_name"_solv_ions.gro -p topol.top -pname NA -nname CL -neutral <<< 13
checkErr "$?" "genion"

gmx grompp -f minim.mdp -c "$protein_name"_solv_ions.gro -p topol.top -o "$protein_name"_em.tpr
checkErr "$?" "grompp_vaccum"
gmx mdrun -v -deffnm "$protein_name"_em
checkErr "$?" "mdrun_vaccume"

gmx grompp -f nvt.mdp -c "$protein_name"_em.gro -p topol.top -o "$protein_name"_nvt.tpr
checkErr "$?" "grompp_nvt"
gmx mdrun -v -deffnm "$protein_name"_nvt
checkErr "$?" "mdrun_nvt"

gmx grompp -f npt.mdp -c "$protein_name"_nvt.gro -t "$protein_name"_nvt.cpt -p topol.top -o "$protein_name"_npt.tpr
checkErr "$?" "grompp_npt"
gmx mdrun -deffnm "$protein_name"_npt
checkErr "$?" "mdrun_npt"

cp posre.itp posre.itp_bkp
#(delete gaps resi)
./delLoopsPosre.sh gaps.txt "$protein_name"_npt.gro
checkErr "$?" "delLoopsPosre"

gmx grompp -f md_edit.mdp -c "$protein_name"_npt.gro -t "$protein_name"_npt.cpt -p topol.top -o "$protein_name"_md.tpr
checkErr "$?" "grompp_md"
gmx mdrun -v -deffnm "$protein_name"_md
checkErr "$?" "mdrun_md"
