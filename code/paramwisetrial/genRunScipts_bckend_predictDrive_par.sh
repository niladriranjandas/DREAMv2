#!/bin/bash

function usage(){
	echo "./genRunScripts_bckend_predictDrive <protein-name> <protein-name_param1,protein-name_param2,...>"
}


### Global params ###
logfile="../log/$1_allparam_log.txt"
#protein_path="../protein/$1"
LOG_MASTER=$1"_allparam_minimal_log.txt"

### call the precursor to consolidate for all param ###
files=""
IFS=', ' read -r -a array <<< "$2"

echo "---- begin $1 ----" > $logfile
for filei in "${array[@]}"
do
	printf "\n\t\t %s stage1precursor: div conquer only predicted" "$filei" >> $logfile
	#./genRunScipts_stage1precursor.sh $filei
	paramwisetrial/genRunScipts_stage1precursor_par.sh $filei
        if [ $? -ne 0 ]; then
		exit -1
	fi
	matlab -nodesktop -nosplash < "stage1_precursor_$filei.m" &
	echo "stage1 precursor for stae1_precursor_$filei.m: $?" >> $LOG_MASTER
	files=$files",$filei/stage1_$filei.mat"
	printf "\n\t\t %s/stage1_%s.mat" "$filei" "$filei" >> $logfile
done
wait

for filei in "${array[@]}"
do
	# --- coppy the files ---
	#cp "stage1_precursor_$filei.m" ../protein/$filei/.
	mv "stage1_precursor_$filei.m" ../protein/$filei/.
done

files="${files:1}"

### do the dedrogram thing ###
# 1. generate dendrogram
# 2. for each dendrogram not solved find leaves to run 
# 3. Start with longest array of leaf .. mark solved or not as we go
# 4. save the result for each
filename_forall="stage1_post_reggraph_$1.txt"
flagfile_pred="$1_stage1_solvedflag.txt"

filename_dendrogramrun="dendrogram_$1.m"
printf "\n\n \t\t %s doDendogram: i/p:%s o/p:%s" "$filename_dendrogramrun" "$files" "$filename_forall" >> $logfile
echo "cd .." > $filename_dendrogramrun
echo "addpath(genpath(pwd));" >> $filename_dendrogramrun
echo "cd code;" >> $filename_dendrogramrun
echo "doDendrogram('$files','$filename_forall','$flagfile_pred')" >> $filename_dendrogramrun # write the leaves in a file
echo "exit" >> $filename_dendrogramrun
matlab -nodesktop -nosplash < "$filename_dendrogramrun"
echo "do dengrogram $filename_dendrogramrun: $?" >> $LOG_MASTER
cat $filename_forall >> $logfile

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
			matlab -nodesktop -nosplash < $filename &
			echo "modelling fragments and consolidating $filename" >> $LOG_MASTER
		fi
	 fi
done < "$filename_forall"
wait

while IFS= read -r line || [[ -n "$line" ]]; do
         filename="stage1_$line.m"
	 # ----------- move results for stage-I ------------------
         #cp $filename ../protein/$line/.
	 mv $filename ../protein/$line/.
done < "$filename_forall"

# ---- gather the params which led to solns -------
filetmp="$1_findwhichtorun.txt"
file_whichrun="$1_stage_intermediate_1_2.m"
echo "cd .." > $file_whichrun
echo "addpath(genpath(pwd));" >> $file_whichrun
echo "cd code;" >> $file_whichrun
#echo "runSecondStage2('$files','$filetmp')" >> $file_whichrun
echo "runSecondStage2_v2('$files','$filetmp')" >> $file_whichrun
echo "exit" >> $file_whichrun
matlab -nodesktop -nosplash < $file_whichrun
echo "---------------------Consolidate pdb generated for files:" >> $logfile
cat $filetmp >> $logfile

# ---------------------------------------------------------
while IFS= read -r line || [[ -n "$line" ]]; do
	#paramwisetrial/run_stage2_forparams.sh "$line" "$filetmp" "$logfile" &
        paramwisetrial/run_stage2_forparams.sh "$line" "$filetmp" "$logfile"
	echo "running gap correct for $line: $?" >> $LOG_MASTER
done < "$filetmp"
#wait

# make a seperate folder for extra files
mkdir ../protein/$1"_allparam"
#cp $logfile ../protein/$1"_allparam"/.
#cp $filename_dendrogramrun ../protein/$1"_allparam"/.
#cp $filename_forall ../protein/$1"_allparam"/.
cp $filetmp ../protein/$1"_allparam"/.
#cp $file_whichrun ../protein/$1"_allparam"/.
#cp $flagfile_pred ../protein/$1"_allparam"/.
mv $logfile ../protein/$1"_allparam"/.
mv $filename_dendrogramrun ../protein/$1"_allparam"/.
mv $filename_forall ../protein/$1"_allparam"/.
#mv $filetmp ../protein/$1"_allparam"/.
mv $file_whichrun ../protein/$1"_allparam"/.
mv $flagfile_pred ../protein/$1"_allparam"/.
mv $flagfile ../protein/$1"_allparam"/.
