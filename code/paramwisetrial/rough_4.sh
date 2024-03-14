#!/bin/bash

filename_forall=stage1_post_reggraph_2o4e.txt
logfile="abc.txt"
flagfile_pred="$1_stage1_solvedflag.txt"

flagfile="$1_tmpflagfile.txt"
echo "NO" > $flagfile
old_filename="";
while IFS= read -r line || [[ -n "$line" ]]; do
	 printf "\n============================================\n"  >> $logfile
	 echo $line  >> $logfile
	 cat $flagfile  >> $logfile	 
	 filename="stage1_$line.m"
	 echo "cd .." > $filename
	 echo "addpath(genpath(pwd));" >> $filename
	 echo "cd code;" >> $filename

	 echo "	fid_flagfile=fopen('$flagfile','w') " >> $filename
	 # ------ localize --------------------------
	 echo "load('../protein/$line/stage1_$line.mat')" >> $filename
	 echo "localize_after_div" >> $filename
	 # ------ global register ---------------------------
	 echo "curr_path = pwd;" >> $filename
	 #echo "cd /home/niladri/Documents/matlab_extras/cvx/" >> $filename
	 echo "cd ../packages/cvx/" >> $filename
	 echo "cvx_setup" >> $filename     #included to run in nmr machines
	 echo "cvx_startup" >> $filename
	 echo "cd (curr_path)" >> $filename
	 echo "consolidate" >> $filename
	 echo "chk_valid_structure">> $filename    # set flag is_valid_structure 
	 echo "if is_valid_structure" >> $filename
	 echo "    structure_done=1;" >> $filename
	 echo "    fprintf(fid_flagfile,'DONE')" >> $filename
	 echo "else" >> $filename
	 echo "    structure_done=0;" >> $filename
	 echo "    fprintf(fid_flagfile,'NO')" >> $filename
	 echo "end" >> $filename
	 echo "fclose(fid_flagfile)" >> $filename
	 # ------- make scripts for stage-2 initilization ------
	 echo "%-----prepare scripts for modeller --------" >> $filename
	 echo "prepForFillGap" >> $filename

	 #-------------- save the workspace ---------------
	 echo "%--------- save workspace -----------------" >> $filename
	 echo "save('../protein/$line/stage1_$line.mat')" >> $filename
	 echo "exit" >> $filename

	 # ----------- run matlab for stage-I --------------------
	 runornot=`awk '$0=="DONE"' $flagfile`
	 if [ -z "$runornot" ]; then
		if [ -s "$flagfile_pred" ]; then
			printf "\n \t\t\t Skipping Running:%s since predicted case ran\n" "$filename" >> $logfile
			echo "-------------------Prev run result" >> $logfile		
			cat $flagfile >> $logfile
			
		else		
			printf "\n \t\t\t Running:%s\n" "$filename" >> $logfile
			echo "-------------------Prev run result" >> $logfile		
			cat $flagfile >> $logfile
			echo "matlab -nodesktop -nosplash < $filename"
		fi
	 fi
	 # ----------- move results for stage-I ------------------
         cp $filename ../protein/$line/.
done < "$filename_forall"

