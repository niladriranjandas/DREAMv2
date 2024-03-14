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
basedir=".."
  if [ -f $basedir/protein/$1/$1.xml ]; then
       echo "Found PARAM file" 
  else
       echo "ERROR: PARAM file not found"
       exit 1
  fi

 # call a small function to check all params are present

 # -----stage-I master file 
  filename="stage1_precursor_$1.m"
       if [ -f $filename ];then
             echo "DELETING $filename"
             rm -f $filename
       fi

  echo "cd .." > $filename
  echo "addpath(genpath(pwd));" >> $filename
  echo "cd code;" >> $filename
  echo "%------INPUTS and PARAMETERS-------" >> $filename
   
 # ----read param file and append to stage-I .m script---
 param_file=$basedir/protein/$1/$1.xml
 #./front_end/xml_xmlread/xmldom $param_file >> $filename
 #echo " " >> $filename
 ./front_end/parseXmlNew.sh $param_file >> $filename

 # ------ set paths for output files ---------
  echo " " >> $filename
  echo "%--------- OUTPUT PARAMS ----------" >> $filename
  echo "OUTPUT.local_folder_pdb           = '$basedir/output_pdb';" >> $filename
  echo "OUTPUT.local_folder_plots_postreg = '$basedir/output_plots';" >> $filename
  echo "OUTPUT.modeller_scripts           = 'modeller_stuff';" >> $filename
  echo " " >> $filename

 # ------ set up paths for logs and workspace save --------
  echo " " >> $filename
  echo "%--------- Log file paths ---------" >> $filename
  echo "LOG.div_n_conquer  = '$basedir/log';" >> $filename
  echo "LOG.considate      = '$basedir/log';" >> $filename
  LOG_modeller="$basedir/log/modeller/"
  LOG_gromacs="$basedir/log/gromacs/"
  echo "%----------workspace save ---------" >> $filename
  echo "SAVES.wrkspace = '$basedir/workspace_sace';" >> $filename
    
 # ---- wrote the omit delete ----------------
  if [ "$#" -eq 2 ]; then
	  echo " " >> $filename
	  echo "%-------- Delete upper bound parameter ------" >> $filename
	  echo "omit.percent = "$omit_percent >> $filename
  else
	  echo "omit.percent = 0" >> $filename
  fi

 # ----- write the preprocess module call ----------
 echo "%------preprocess and div-and-conquer -----" >> $filename
 echo "div_dont_conquer" >> $filename 
 echo "do_prediction" >> $filename
 echo "if consolidate_flag" >> $filename  # set by do_prediction
	  #echo "localize_after_div_par" >> $filename    # matlab internal exception - unable to assign threads
	  echo "localize_after_div" >> $filename
	  # ------ global register --------------------------
	  echo "curr_path = pwd;" >> $filename
	  #echo "cd /home/niladri/Documents/matlab_extras/cvx/" >> $filename
	  echo "cd $basedir/packages/cvx/" >> $filename
	  echo "cvx_setup" >> $filename     #included to run in nmr machines
	  echo "cvx_startup" >> $filename
	  echo "cd (curr_path)" >> $filename
	  echo "consolidate" >> $filename
	  echo "chk_valid_structure">> $filename    # set is_valid_structure flag
          echo "if is_valid_structure" >> $filename
	  echo "     structure_done=1;"  >> $filename
	  echo "else" >> $filename
	  echo "     structure_done=0;" >> $filename
	  echo "end" >> $filename
	 # ------- make scripts for stage-2 initilization ------
	 echo "%-----prepare scripts for modeller --------" >> $filename
	 echo "prepForFillGap" >> $filename
echo "end" >> $filename

 # if solve write in $filename_for_all
 # ---- save it ---- #
 echo "save('$basedir/protein/$1/stage1_$1.mat')" >> $filename
 echo "exit" >> $filename

