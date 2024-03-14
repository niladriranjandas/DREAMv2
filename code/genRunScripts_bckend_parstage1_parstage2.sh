#!/bin/bash

   if [ -z "$1" ]; then
        echo "error: usage ./genRunScripts_bckend.sh <protein_name> [omit=<value(percentage)>]"
        exit 1
   fi
  
   if [ "$#" -eq 2 ]; then
	omit_percent=`echo $2 | awk -F'omit=' '{print $2}'`
	if [ -z $omit_percent ]; then
		echo "error: usage ./genRunScripts_bckend.sh <protein_name> omit=<value>"
		exit 1
	fi
   fi		

#---- check if all files are in order ---

  if [ -f ../protein/$1/$1.xml ]; then
       echo "Found PARAM file" 
  else
       echo "ERROR: PARAM file not found"
       exit 1
  fi

  forLogpath=`pwd`
  LOG_MASTER=$forLogpath"/"$1"_minimal_log.txt" 
# call a small function to check all params are present

 # -----stage-I master file 
  filename="stage1_$1.m"

  filename_preloc="stage1_$1_preloc.m"

       if [ -f $filename ];then
             echo "DELETING $filename"
             rm -f $filename
       fi
       if [ -f $filename_preloc ];then
             echo "DELETING $filename_preloc"
             rm -f $filename_preloc
       fi

  echo "cd .." > $filename_preloc
  echo "addpath(genpath(pwd));" >> $filename_preloc
  echo "cd code;" >> $filename_preloc
  echo "%------INPUTS and PARAMETERS-------" >> $filename_preloc
 

 # ----read param file and append to stage-I .m script---
 param_file=../protein/$1/$1.xml
 #./front_end/xml_xmlread/xmldom $param_file >> $filename
 #echo " " >> $filename
 ./front_end/parseXmlNew.sh $param_file >> $filename_preloc
 echo "Parsing parameter file $param_file: $?" >> $LOG_MASTER
 # ------ set paths for output files ---------
  echo " " >> $filename_preloc
  echo "%--------- OUTPUT PARAMS ----------" >> $filename_preloc
  echo "OUTPUT.local_folder_pdb           = '../output_pdb';" >> $filename_preloc
  echo "OUTPUT.local_folder_plots_postreg = '../output_plots';" >> $filename_preloc
  echo "OUTPUT.modeller_scripts           = 'modeller_stuff';" >> $filename_preloc
  echo " " >> $filename_preloc

 # ------ set up paths for logs and workspace save --------
  echo " " >> $filename_preloc
  echo "%--------- Log file paths ---------" >> $filename_preloc
  echo "LOG.div_n_conquer  = '../log';" >> $filename_preloc
  echo "LOG.considate      = '../log';" >> $filename_preloc
  LOG_modeller="../log/modeller/"
  LOG_gromacs="../log/gromacs/"
  echo "%----------workspace save ---------" >> $filename_preloc
  echo "SAVES.wrkspace = '../workspace_sace';" >> $filename_preloc
  
   # ---- wrote the omit delete ----------------
  if [ "$#" -eq 2 ]; then
	  echo " " >> $filename_preloc
	  echo "%-------- Delete upper bound parameter ------" >> $filename_preloc
	  echo "omit.percent = "$omit_percent >> $filename_preloc
  else
	  echo "omit.percent = 0" >> $filename_preloc
  fi
  

 # ----- write the preprocess module call ----------
 echo "%------preprocess and div-and-conquer -----" >> $filename_preloc
 #echo "div_conquer" >> $filename
 echo "div_conquer_preloc" >> $filename_preloc
 echo "exit" >> $filename_preloc

 # ----- do stage-1 as separate process ------------
 if [ -d "localize_tmp/$1" ]; then
    new_loc_tmp="localize_tmp/$1_bkp"
    mv "localize_tmp/$1" $new_loc_tmp
 fi
 echo  "mkdir localize_tmp/$1"
 mkdir "localize_tmp/$1"
 matlab -nodesktop -nosplash < $filename_preloc
 echo "Run divide into fragments and model $filename_preloc: $?" >> $LOG_MASTER

 filetorun_stage1="localize_tmp/$1/$1_run_to_localize.txt"

 if [ -f $filetorun_stage1 ]; then
        while IFS= read -r line || [[ -n "$line" ]]; do
            matlab -nodesktop -nosplash < "$line" &
        done < "$filetorun_stage1"
        wait
 else
        echo "ERROR: could not find file $filetorun_stage1"
        exit -1
 fi

 # ------ assimilate the files generated in stage-1 register ---
 filename_assimilate="localize_tmp/$1/assimilate_$1.m"
 file_from_preloc="localize_tmp/$1/preloc_$1.mat"
 matfilesloc="localize_tmp/$1/$1_after_localize.txt"
 if [ ! -f "$file_from_preloc" ]; then
    echo "ERROR: file from preloc $file_from_preloc not found"
    exit -1
 fi
 echo "cd .." > $filename_assimilate
 echo "addpath(genpath(pwd));" >> $filename_assimilate
 echo "cd code;" >> $filename_assimilate
 echo "load('$file_from_preloc');" >> $filename_assimilate
 echo "localize = assimilate_after_seplocalize('$matfilesloc', test.dup_flag)" >> $filename_assimilate
        # ------ global register --------------------------
 echo "curr_path = pwd;" >> $filename_assimilate
 echo "cd ../packages/cvx/" >> $filename_assimilate
 echo "cvx_startup" >> $filename_assimilate
 echo "cd (curr_path)" >> $filename_assimilate
 echo "consolidate" >> $filename_assimilate

 # ------- make scripts for stage-2 initilization ------
 echo "%-----prepare scripts for modeller --------" >> $filename_assimilate
 echo "prepForFillGap" >> $filename_assimilate

 #-------------- save the workspace ---------------
 echo "%--------- save workspace -----------------" >> $filename_assimilate
 echo "save('../protein/$1/stage1_$1.mat')" >> $filename_assimilate
 echo "exit" >> $filename_assimilate

 # ----------- run matlab for stage-I --------------------
 matlab -nodesktop -nosplash < $filename_assimilate 
 echo "Consolidating the fragmented model $filename_assimilate: $?" >> $LOG_MASTER
 # ----------- show results for stage-I ------------------
# ----------- check to see if grp register file generated --
 if [ ! -f  ../output_pdb/grp_$1.pdb ]; then
     echo "Stage-1 group registration file not generated in ../output_pdb/grp_"$1".pdb"
     exit 1
 else
     cp ../output_pdb/grp_$1.pdb modeller_stuff/.
 fi

 
 alifile=modeller_stuff/alignment_"$1".ali
 pyfile=modeller_stuff/do_loop_dntmove_"$1".py
 if [ ! -f "$alifile" ]; then
     echo "Alignment file not generated for modeller in modeller_stuff/alignment.ali"
     exit 1
 fi

 if [ ! -f "$pyfile" ]; then
      echo "Python script not generated for modeller in modeller_stuff/do_loop_dntmove.py" 
      exit 1
 fi

  
 # ----read param file and append to stage-I .m script---
 echo "==============Running Modeller==========================="
 cd modeller_stuff
# mod9.19 do_loop_dntmove.py
# mod9.23 do_loop_dntmove_"$1".py
# mod10.0 do_loop_dntmove_"$1".py
 $MOD do_loop_dntmove_"$1".py
 echo "Running gap filling do_loop_dntmove_$1: $?" >> $LOG_MASTER
 cd ..

# --------- gather the structures from modeller run reducer ----
files=`ls modeller_stuff/ | egrep ^grp_${1}_fill\.B[0-9]+\.pdb$`
     if [ `echo $files|wc -c` = 1 ]; then
             echo "ERROR: modeller files not found"
             exit 1
     fi

#filename="stage2_$1.m"
#     if [ -f $filename ];then
#           echo "DELETING $filename"
#           rm -f $filename
#     fi
#echo "%% ---------------- Do anchored localization -----------" > $filename
#c=0
#for i in `echo $files`
#do
#   c=$((c+1))
#echo $i
##wc -c $i
#  # reducer_fill_H -FLIP modeller_stuff/$i > ../output_pdb/grp_$1_${c}_H.pdb
#  #../packages/Reduce_fill_H_protein/reduce.3.23.130521.linuxi386/reducer_fill_H -FLIP  modeller_stuff/$i > ../output_pdb/grp_$1_${c}_H.pdb
#  line="$1"
#  ../packages/Reduce_fill_H_protein/reduce.3.23.130521.linuxi386/reducer_fill_H -FLIP  modeller_stuff/$i > ../output_pdb/grp_"$line"_"${c}"_H.pdb
#  echo "cd .." >> $filename
#  echo "addpath(genpath(pwd));" >> $filename
#  echo "cd code;" >> $filename
#   echo "load('../protein/$1/stage1_$1.mat')" >> $filename
#   #echo "OUTPUT.modeller_n_H  = '../output_pdb/grp_$1_${c}_H.pdb';" >> $filename
#   hfile=grp_"$line"_"${c}"_H.pdb
#   echo "OUTPUT.modeller_n_H  = '../output_pdb/$hfile';" >> $filename
#   echo "OUTPUT.model_count = ${c};" >> $filename
#   echo "doAnchoredLoc;" >> $filename
#   echo "save('../protein/$1/stage2_$1 _${c}.mat')" >> $filename
#   echo "[anclocrefarr, anclocref] = doMultiAnclocRefine(wh_eq_cons_all, wh_up_bounds, wh_lo_bounds, fillQ, wh_Comp,'../output_pdb');" >> $filename
#   echo "save('../protein/$1/stage2_$1 _${c}_mutiref.mat')" >> $filename
#   echo "%========================================================" >> $filename
#   echo " " >> $filename
#   echo " "
#done

## ------------ run stage -2 script in matlab ---------------------
#   matlab -nosplash -nodesktop < $filename
## ------------- 

allfiles_stage2=""
filename_suffix="stage2_$1"
c=0
for i in `echo $files`
do
   c=$((c+1))
   filename=$filename_suffix"_${c}.m"
   allfiles_stage2=$allfiles_stage2$filename","
echo $i
   line="$1"
   ../packages/Reduce_fill_H_protein/reduce.3.23.130521.linuxi386/reducer_fill_H -FLIP  modeller_stuff/$i > ../output_pdb/grp_"$line"_"${c}"_H.pdb
   if [ -f $filename ];then
         echo "DELETING $filename"
          rm -f $filename
   fi
   echo "%% ---------------- Do anchored localization -----------" > $filename
   echo "cd .." >> $filename
   echo "addpath(genpath(pwd));" >> $filename
   echo "cd code;" >> $filename
   echo "load('../protein/$1/stage1_$1.mat')" >> $filename
   #echo "OUTPUT.modeller_n_H  = '../output_pdb/grp_$1_${c}_H.pdb';" >> $filename
   hfile=grp_"$line"_"${c}"_H.pdb
   echo "OUTPUT.modeller_n_H  = '../output_pdb/$hfile';" >> $filename
   echo "OUTPUT.model_count = ${c};" >> $filename
   echo "doAnchoredLoc;" >> $filename
   echo "save('../protein/$1/stage2_$1 _${c}.mat')" >> $filename
   echo "[anclocrefarr, anclocref] = doMultiAnclocRefine(wh_eq_cons_all, wh_up_bounds, wh_lo_bounds, fillQ, wh_Comp,'../output_pdb');" >> $filename
   echo "save('../protein/$1/stage2_$1 _${c}_mutiref.mat')" >> $filename
   echo "%========================================================" >> $filename
   echo " " >> $filename
   echo " "
done
        
# ------------ run stage -2 script in matlab ---------------------
allfiles_stage2=${allfiles_stage2::-1}
echo "$allfiles_stage2"

IFS=, read -ra allfiles_stage2_arr <<< "$allfiles_stage2"
for file_i_stage2 in "${allfiles_stage2_arr[@]}"
do
   matlab -nosplash -nodesktop < $file_i_stage2 &  
done
wait
echo "Anchored localization $allfiles_stage2: $?" >> $LOG_MASTER
###################################################


upls_1=`grep upl_file $param_file | awk -F'<upl_file>' '{print $2}' | awk -F'</upl_file>' '{print $1}'`
upls_hbo=`grep hbond_file $param_file | awk -F'<hbond_file>' '{print $2}' | awk -F'</hbond_file>' '{print $1}'`

  if [ ! -z $upls_hbo ]; then
    upls=$upls_1','$upls_hbo
  else
    upls=$upls_1
  fi

  IFS=', ' read -r -a array <<< "$upls"

  new_file=../protein/$1/$1"_concatupl.upl"
  em_file=../protein/$1/$1"_emscript.sh"
  if [ -f $new_file ]; then
    echo "Deleting old concatenated upl file: $new_file"
    rm -f $new_file
  fi
  for files in "${array[@]}"
  do
     cat ../protein/$1/$files >> $new_file
  done

  # --- changed Jun 18 21 --- #
  #./rungenEMscripts.sh $1 ../output_pdb $new_file | tail -1 > $em_file 
  #em_run_cmd=`./rungenEMscripts.sh $1 ../output_pdb $new_file | tail -1`
  #./$em_run_cmd
  em_file1=../protein/$1/$1"_emscript_beforemd.sh"
  ./rungenEMscripts_beforemd.sh $1 ../output_pdb $new_file | tail -1 > $em_file1
  em_run_cmd1=`./rungenEMscripts_beforemd.sh $1 ../output_pdb $new_file | tail -1`
  ./$em_run_cmd1

  em_file2=../protein/$1/$1"_emscript_md.sh"
  ./rungenEMscripts_onlymd.sh $1 ../output_pdb $new_file | tail -1 > $em_file2
  em_run_cmd2=`./rungenEMscripts_onlymd.sh $1 ../output_pdb $new_file | tail -1`
  ./$em_run_cmd2

  echo "Ending anchored locatlization $em_run_cmd2: $?" >> $LOG_MASTER
  # ------------------------- #
# -----------------------------------------------------------------------
## -- move the folders -- ##
  #- output_pdbs
  mkdir ../protein/$1/output_pdb
  mv ../output_pdb/*$1* ../protein/$1/output_pdb/.
  #- output_plots   no plots are generated as of now
  #mkdir ../protein/$i/output_plots
  #mv ../output_plots/*$1* ../protein/$1/output_plots/.
  #- logs
  mkdir ../protein/$1/log
  mv ../log/*$1* ../protein/$1/log/.

  #-md_related
  #mkdir ../protein/$1/md_related
  #cp -r md_related/$1/* ../protein/$1/md_related/.
  mkdir ../protein/$1/md_related_v2
  #cp -r md_related_v2/$1/* ../protein/$1/md_related_v2/.   # DONOT change this to mv .. callBadGaps.sh will fail
  #-move temp dssp files
  mkdir ../protein/$1/tmp_dssp_files
  mv *$1*dssp.txt ../protein/$1/tmp_dssp_files/.
 
  #- matlab scripts
  #mv stage1_$1.m ../protein/$1/.
  mv "stage1_"$1"_preloc.m" ../protein/$1/.
  mv stage2_$1_*.m ../protein/$1/.

  #- modeller files and scripts
  mkdir ../protein/$1/modeller
  mv modeller_stuff/*$1* ../protein/$1/modeller/.

