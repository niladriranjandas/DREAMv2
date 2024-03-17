#!/bin/bash

# Copyright (c) 2024 niladriranjandas
# clean up files created and transfer them to the proteins folder

protname="$1"

alreadyrun="../proteinsparams/listoffiles.txt"
em_folder="md_related_v2"

backup_folder="../protein/$protname"

### folders to backup ###
md_beforegapcor=md_related_v2/
md_aftergapcor=md_related_gapcorrect_v2/
mod20_gapcor=gap_correct/MOD20/
modmulti_gapcor=gap_correct/MODMULTI/
residual_gapcor=gap_correct/RESIDUAL_FILES/
waterref_folder=waterrefine/
localize_tmp_folder=localize_tmp/


# -------------- #

 protfolder="$protname"
 backup_folder="../protein/$protfolder"
 
 count=0
 mkdir "$backup_folder"/gap_correct
 for pdbfile in `ls $em_folder/$protfolder/*/*em_nosol.pdb`
 do      
   	count=$((count+1))
         
   	pdbname=`basename $pdbfile`
   	echo "#####  $pdbname ####"
   	echo " "
   	pdbname_noext=${pdbname::-4}
 	## folder waterrefine ##
 	if [[ $count -eq 1 ]]; then
 		mkdir "$backup_folder"/waterrefine
 	fi
 	mv "$waterref_folder"/"$pdbname_noext" "$backup_folder"/waterrefine/.
 			
 	## folder gap_correct/MOD20 ##
 	if [[ $count -eq 1 ]]; then
 		mkdir "$backup_folder"/gap_correct/MOD20
 	fi
 	mv "$mod20_gapcor"/"$pdbname_noext" "$backup_folder"/gap_correct/MOD20/.
 	
 	## folder gap_correct/MODMULTI ##
 	if [[ $count -eq 1 ]]; then
 		mkdir "$backup_folder"/gap_correct/MODMULTI
 	fi
 	mv "$modmulti_gapcor"/"$pdbname_noext" "$backup_folder"/gap_correct/MODMULTI/.
 	
 	## folder gap_correct/RESIDUAL_FILES ##
 	if [[ $count -eq 1 ]]; then
 		mkdir "$backup_folder"/gap_correct/RESIDUAL_FILES
 	fi
 	mv "$residual_gapcor"/"$pdbname_noext" "$backup_folder"/gap_correct/RESIDUAL_FILES/.
 	
 	## folder md_related_gapcorrect_v2 ##
 	if [[ $count -eq 1 ]]; then
 		mkdir "$backup_folder"/md_related_gapcorrect_v2
 	fi
 	mv "$md_aftergapcor"/"$pdbname_noext" "$backup_folder"/md_related_gapcorrect_v2/.
 	
 done
 ## folder md before gap_correct
 mkdir "$backup_folder"/md_related_v2
 mv "$md_beforegapcor"/"$protfolder" "$backup_folder"/md_related_v2/.	
 ## folder localize_tmp/
 mkdir "$backup_folder"/localize_tmp
 mv $localize_tmp_folder/$protfolder "$backup_folder"/localize_tmp/.
 ## move stage1_<protname>_preloc 
 stage1preloc="stage1_"$protfolder"_preloc.m"
 [ -f "$stage1preloc" ] &&  mv "$stage1preloc" "$backup_folder"/.
# -------------- #
