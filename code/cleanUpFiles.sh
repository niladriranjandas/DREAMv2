#!/bin/bash
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

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
		
        protfolder=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`        
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
	[ -f "$stage1preloc" ] && mv "$stage1preloc" "$backup_folder"/.
else
	#echo "code it"
	mdbefore_whichfiles=$protname"_findwhichtorun.txt"
	
	backupfolder_new=$backup_folder"_allparam"
	
	# md_related_v2
	mkdir "$backupfolder_new"/md_related_v2
	while IFS= read -r line || [[ -n "$line" ]]; do
		mv $md_beforegapcor/$line "$backupfolder_new"/md_related_v2/.
	done < "$mdbefore_whichfiles"
	
	waterref_whichfiles=$protname"_whichwaterref_mat.txt"
	waterref_which_folder="gap_correct/RESIDUAL_FILES/$protname"
	MAXFILE=5
	# md_related_gapcorrect_v2
	mkdir "$backupfolder_new"/md_related_gapcorrect_v2	
	for whichfiles in `head -$MAXFILE $waterref_which_folder/$waterref_whichfiles`
	do
		foldername=`basename $whichfiles`
		foldername=${foldername::-4}		
		mv $md_aftergapcor/$foldername "$backupfolder_new"/md_related_gapcorrect_v2/.
	done
	
	# waterrefine
	mkdir "$backupfolder_new"/waterrefine
	for whichfiles in `head -$MAXFILE $waterref_which_folder/$waterref_whichfiles`
	do
		foldername=`basename $whichfiles`
		foldername=${foldername::-4}		
		mv $waterref_folder/$foldername "$backupfolder_new"/waterrefine/.
	done	
	
	## gap correct
	mkdir "$backupfolder_new"/gap_correct
	# gap_correct/MOD20
	mkdir "$backupfolder_new"/gap_correct/MOD20
	for whichfiles in `head -$MAXFILE $waterref_which_folder/$waterref_whichfiles`
	do
		foldername=`basename $whichfiles`
		foldername=${foldername::-4}		
		mv $mod20_gapcor/$foldername "$backupfolder_new"/gap_correct/MOD20/.
	done	
	# gap_correct/MODMULTI
	mkdir "$backupfolder_new"/gap_correct/MODMULTI
	for whichfiles in `head -$MAXFILE $waterref_which_folder/$waterref_whichfiles`
	do
		foldername=`basename $whichfiles`
		foldername=${foldername::-4}		
		mv $modmulti_gapcor/$foldername "$backupfolder_new"/gap_correct/MODMULTI/.
	done
	# gap_correct/RESIDUAL_FILES
	mkdir "$backupfolder_new"/gap_correct/RESIDUAL_FILES
	for whichfiles in `head -$MAXFILE $waterref_which_folder/$waterref_whichfiles`
	do
		foldername=`basename $whichfiles`
		foldername=${foldername::-4}		
		mv $residual_gapcor/$foldername "$backupfolder_new"/gap_correct/RESIDUAL_FILES/.
	done
	
	## gap_correct/RESIDUAL_FILES/$protname
	mv gap_correct/RESIDUAL_FILES/$protname "$backupfolder_new"/gap_correct/RESIDUAL_FILES/.
	
	## finally transfer the txt files ##
	mv $mdbefore_whichfiles "$backupfolder_new"/.
	
	## move _log_all.log.txt ##
        log_file=$protname"_log_all.log.txt"        
	mv paramwisetrial/$log_file "$backupfolder_new"/.

	### move files not there in $mdbefore_whichfiles
	mkdir "$backupfolder_new"/extra_files
	# modeller_stuff
	mkdir "$backupfolder_new"/extra_files/MODELLER_STUFF
	mv modeller_stuff/*_$protname_* "$backupfolder_new"/extra_files/MODELLER_STUFF/.
	# output_pdb
	mkdir "$backupfolder_new"/extra_files/output_pdb
	mv ../output_pdb/*_$protname_* "$backupfolder_new"/extra_files/output_pdb/.
	# log
	mkdir "$backupfolder_new"/extra_files/log
	mv ../log/*_$protname_* "$backupfolder_new"/extra_files/log/.

        ## move stage1_<protname>_preloc 
        #stage1preloc="stage1_"$protname"_preloc.m"
        #mv "$stage1preloc" "$backup_folder"/.
fi	
