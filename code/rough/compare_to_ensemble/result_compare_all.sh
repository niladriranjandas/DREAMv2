#!/bin/bash

#-----------------------------------------------#
# generates the list to be read for comparision 
# txt file created with entries:                
# generated_pdb_name PDB_downloaded_file_path
#-----------------------------------------------#

result_folder="$1"
org_pdb_folder="$2"

log_rmsd="$3"
rama_plot_folder="$4"
errat_folder="$5"


    for protein in `ls -d $result_folder/*/`
    do
        protein=`echo $protein | awk -F'/' '{print $((NF-1))}'`
	protein_name=`echo $protein | awk -F".pdb" '{print $1}'`
        protein_upper=`echo $protein_name | tr [:lower:] [:upper:]`
        org_protein=`find "$org_pdb_folder" -name "$protein_upper.pdb"`
        
        for gen_pdbs in `ls *.pdb`
        do
            # gen_pdb gen_pdb_path              org_pdb       org_protein_path
            # gen_pdb $result_folder/$protein   $protein_upper $org_protein
	    gen_pdbs_name=`echo $gen_pdbs | awk -F".pdb" '{print $1}'`            

printf "\n %s %s/%s %s %s" $gen_pdb $result_folder $protein $protein_upper $org_protein	

            # 1) best rmsd 
            #pymol -cq compare_truth_printbest.py -- $gen_pdbs_name $result_folder/$protein $protein_upper $org_protein >> $log_rmsd
            
            # 2) ramachandran map
            filename_rama="$gen_pdbs".png
            #./drawRamaFromPDB.sh $gen_pdbs $rama_plot_folder/$filename_rama
            
            # 3) errat
            #./errat 
            
            # 4) residue-residue contact map
        done
    done
    



