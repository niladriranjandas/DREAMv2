#!/bin/bash

function usage(){
	echo "./genRunScripts_bckend_predictDrive <protein-name> <protein-name_param1,protein-name_param2,...>"
}


### Global params ###
logfile="../log/$1_allparam_log.txt"
#protein_path="../protein/$1"


### call the precursor to consolidate for all param ###
files=""
IFS=', ' read -r -a array <<< "$2"

echo "---- begin $1 ----" > $logfile
for filei in "${array[@]}"
do
	printf "\n\t\t %s stage1precursor: div conquer only predicted" "$filei" >> $logfile
	#./genRunScipts_stage1precursor.sh $filei
	paramwisetrial/genRunScipts_stage1precursor.sh $filei
        if [ $? -ne 0 ]; then
		exit -1
	fi
	matlab -nodesktop -nosplash < "stage1_precursor_$filei.m"		
	files=$files",$filei/stage1_$filei.mat"
	printf "\n\t\t %s/stage1_%s.mat" "$filei" "$filei" >> $logfile
done

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
			matlab -nodesktop -nosplash < $filename
		fi
	 fi
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
	  param_file=../protein/$line/$line.xml
 	 # ----------- check to see if grp register file generated --
	 if [ ! -f  ../output_pdb/grp_$line.pdb ]; then
     		echo "Stage-1 group registration file not generated in ../output_pdb/grp_"$line".pdb"
		exit 1
	 else
	        #cp ../output_pdb/grp_$line.pdb modeller_stuff/.
	        mv ../output_pdb/grp_$line.pdb modeller_stuff/.
	 fi

 
 
	 if [ ! -f modeller_stuff/alignment_"$line".ali ]; then
	     echo "Alignment file not generated for modeller in modeller_stuff/alignment_$line.ali"
	     exit 1
	 fi

	 if [ ! -f modeller_stuff/do_loop_dntmove_"$line".py ]; then
	      echo "Python script not generated for modeller in modeller_stuff/do_loop_dntmove_$line.py" 
	      exit 1
	 fi

	 # ---------- do Modeller run ---------------------------------
	 printf "\n \t\t\t Run Modeller for %s " "$line" >> $logfile
	 echo "==============Running Modeller==========================="
	 cd modeller_stuff
	 #mod9.19 do_loop_dntmove.py
	 #mod9.20 do_loop_dntmove_"$line".py
         #mod9.23 do_loop_dntmove_"$line".py
         mod10.0 do_loop_dntmove_"$line".py
	 cd ..

	 # --------- gather the structures from modeller run reducer ----
	files=`ls modeller_stuff/ | egrep ^grp_${line}_fill\.B[0-9]+\.pdb$`
	if [ `echo $files|wc -c` = 1 ]; then
             echo "ERROR: modeller files not found"
             exit 1
	fi
	cat $files >> $logfile

	filename="stage2_$line.m"
        if [ -f $filename ];then
      	     echo "DELETING $filename"
             rm -f $filename
        fi
	echo "%% ---------------- Do anchored localization -----------" > $filename
	c=0
	for i in `echo $files`
	do
	   c=$((c+1))
	   echo $i
	   wc -c $i
	   # reducer_fill_H -FLIP modeller_stuff/$i > ../output_pdb/grp_$line_${c}_H.pdb
	   #../packages/Reduce_fill_H_protein/reduce.3.23.130521.linuxi386/reducer_fill_H -FLIP  modeller_stuff/$i > ../output_pdb/grp_$line_${c}_H.pdb
           ../packages/Reduce_fill_H_protein/reduce.3.23.130521.linuxi386/reducer_fill_H -FLIP  modeller_stuff/$i > ../output_pdb/grp_"$line"_"${c}"_H.pdb
	   printf "\n \t\t\t Reducer ../output_pdb/grp_%s_%s_H.pdb" "$line" "${c}">> $logfile
	   echo "cd .." >> $filename
	   echo "addpath(genpath(pwd));" >> $filename
	   echo "cd code;" >> $filename
	   echo "load('../protein/$line/stage1_$line.mat')" >> $filename
	   #echo "OUTPUT.modeller_n_H  = '../output_pdb/grp_$line_${c}_H.pdb';" >> $filename
           hfile=grp_"$line"_"${c}"_H.pdb
	   echo "OUTPUT.modeller_n_H  = '../output_pdb/$hfile';" >> $filename
           echo "OUTPUT.model_count = ${c};" >> $filename
	   echo "doAnchoredLoc;" >> $filename
	   echo "save('../protein/$line/stage2_$line _${c}.mat')" >> $filename
	   echo "[anclocrefarr, anclocref] = doMultiAnclocRefine(wh_eq_cons_all, wh_up_bounds, wh_lo_bounds, fillQ, wh_Comp,'../output_pdb');" >> $filename
	   echo "save('../protein/$line/stage2_$line _${c}_mutiref.mat')" >> $filename
	   echo "%========================================================" >> $filename
	   echo " " >> $filename
	   echo " "
	done

	# ------------ run stage -2 script in matlab ---------------------
	   printf "\n \t\t\t anchored locatlization %s" "$filename" >> $logfile
	   matlab -nosplash -nodesktop < $filename
	# ------------- 
	upls_1=`grep upl_file $param_file | awk -F'<upl_file>' '{print $2}' | awk -F'</upl_file>' '{print $1}'`
	upls_hbo=`grep hbond_file $param_file | awk -F'<hbond_file>' '{print $2}' | awk -F'</hbond_file>' '{print $1}'`

	  if [ ! -z $upls_hbo ]; then
	    upls=$upls_1','$upls_hbo
	  else
	    upls=$upls_1
	  fi

	  IFS=', ' read -r -a array <<< "$upls"

	  new_file=../protein/$line/$line"_concatupl.upl"
	  em_file=../protein/$line/$line"_emscript.sh"
	  if [ -f $new_file ]; then
	    echo "Deleting old concatenated upl file: $new_file"
	    rm -f $new_file
	  fi
	  for files in "${array[@]}"
	  do
	     cat ../protein/$line/$files >> $new_file
	  done

	  ./rungenEMscripts.sh $line ../output_pdb $new_file | tail -1 > $em_file 
	  printf "\n \t\t\t EM script:%s " "$em_file" >> $logfile
	  em_run_cmd=`./rungenEMscripts.sh $line ../output_pdb $new_file | tail -1`
	  ./$em_run_cmd
	# -----------------------------------------------------------------------
	## -- move the folders -- ##
	  #- output_pdbs
	  mkdir ../protein/$line/output_pdb
	  #cp ../output_pdb/*$line* ../protein/$line/output_pdb/.
          mv ../output_pdb/*$line* ../protein/$line/output_pdb/.
	  #- output_plots
	  mkdir ../protein/$line/output_plots
	  #cp ../output_plots/*$line* ../protein/$line/output_plots/.
	  mv ../output_plots/*$line* ../protein/$line/output_plots/.
	  #- logs
	  mkdir ../protein/$line/log
	  #cp ../log/*$line* ../protein/$line/log/.
          mv ../log/*$line* ../protein/$line/log/.
	  #- modeller_stuff
	  mkdir ../protein/$line/modeller_stuff
	  #cp modeller_stuff/*$line* ../protein/$line/modeller_stuff/.
	  mv modeller_stuff/*$line* ../protein/$line/modeller_stuff/.
	  #-md_related
	  #mkdir ../protein/$line/md_related
	  #cp -r md_related/$line/* ../protein/$line/md_related/.
	  mkdir ../protein/$line/md_related_v2
	  #cp -r md_related_v2/$line/* ../protein/$line/md_related_v2/.
	  mv -r md_related_v2/$line/* ../protein/$line/md_related_v2/.
	  #-move temp dssp files
	  mkdir ../protein/$line/tmp_dssp_files
	  mv *$line*dssp.txt ../protein/$line/tmp_dssp_files/.

	  # - move stage2 files 
	  #cp $filename ../protein/$line/.
	  mv $filename ../protein/$line/.

	  # - move modeller files
	  mkdir ../protein/$line/modeller
	  #cp modeller_stuff/*$line* ../protein/$line/modeller/.
	  mv modeller_stuff/*$line* ../protein/$line/modeller/.
done < "$filetmp"

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
