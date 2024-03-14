#!/bin/bash

#-------- coppy files from ~/tmp of httpd --------

httpd_folder="/tmp/systemd-private-305dd2c25d304e9ca343f71a63db201b-httpd.service-ZMcLBf/tmp"

   if [ -z "$1" ]; then
        echo "error: usage ./genRunScripts.sh <protein_name>"
        exit 1
   fi

   if [ ! -d "../protein/$1" ]; then
        mkdir "../protein/$1"
        sudo find $httpd_folder -name "$1*" | sudo xargs cp -t ../protein/$1/.
        sudo chmod 777 ../protein/$1/*        
   else
        echo "Backup ../protein/$1/. and re-run."
        exit 1
   fi   


#---- check if all files are in order ---

  if [ -f ../protein/$1/$1.xml ]; then
       echo "Found PARAM file" 
  else
       echo "ERROR: PARAM file not found"
       exit 1
  fi

 # call a small function to check all params are present

 # -----stage-I master file 
  filename="stage1.m"
       if [ -f $filename ];then
             echo "DELETING $filename"
             rm -f $filename
       fi

  echo "cd .." > $filename
  echo "addpath(genpath(pwd));" >> $filename
  echo "cd code;" >> $filename
  echo "%------INPUTS and PARAMETERS-------" >> $filename
   
 # ----read param file and append to stage-I .m script---
 param_file=../protein/$1/$1.xml
 ./front_end/xml_xmlread/xmldom $param_file >> $filename
 echo " " >> $filename

 # ------ set paths for output files ---------
  echo " " >> $filename
  echo "%--------- OUTPUT PARAMS ----------" >> $filename
  echo "OUTPUT.local_folder_pdb           = '../output_pdb';" >> $filename
  echo "OUTPUT.local_folder_plots_postreg = '../output_plots';" >> $filename
  echo "OUTPUT.modeller_scripts           = 'modeller_stuff';" >> $filename
  echo " " >> $filename

 # ------ set up paths for logs and workspace save --------
  echo " " >> $filename
  echo "%--------- Log file paths ---------" >> $filename
  echo "LOG.div_n_conquer  = '../log';" >> $filename
  echo "LOG.considate      = '../log';" >> $filename
  LOG_modeller="../log/modeller/"
  LOG_gromacs="../log/gromacs/"
  echo "%----------workspace save ---------" >> $filename
  echo "SAVES.wrkspace = '../workspace_sace';" >> $filename
    
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
 echo "div_conquer" >> $filename

 # ------ global register --------------------------
  echo "curr_path = pwd;" >> $filename
  echo "cd /home/niladri/Documents/matlab_extras/cvx/" >> $filename
  echo "cvx_startup" >> $filename
  echo "cd (curr_path)" >> $filename
 echo "consolidate" >> $filename

 # ------- make scripts for stage-2 initilization ------
 echo "%-----prepare scripts for modeller --------" >> $filename
 echo "prepForFillGap" >> $filename

 #-------------- save the workspace ---------------
 echo "%--------- save workspace -----------------" >> $filename
 echo "save('../protein/$1/stage1.mat')" >> $filename
 echo "exit" >> $filename

 # ----------- run matlab for stage-I --------------------
matlab -nodesktop -nosplash < $filename 

 # ----------- show results for stage-I ------------------


# ----------- check to see if grp register file generated --
 if [ ! -f  ../output_pdb/grp_$1.pdb ]; then
     echo "Stage-1 group registration file not generated in ../output_pdb/grp_"$1".pdb"
     exit 1
 else
     cp ../output_pdb/grp_$1.pdb modeller_stuff/.
 fi

 
 
 if [ ! -f modeller_stuff/alignment.ali ]; then
     echo "Alignment file not generated for modeller in modeller_stuff/alignment.ali"
     exit 1
 fi

 if [ ! -f modeller_stuff/do_loop_dntmove.py ]; then
      echo "Python script not generated for modeller in modeller_stuff/do_loop_dntmove.py" 
      exit 1
 fi

# ---------- do Modeller run ---------------------------------
 echo "==============Running Modeller==========================="
 cd modeller_stuff
 mod9.19 do_loop_dntmove.py
 cd ..

# --------- gather the structures from modeller run reducer ----
files=`ls modeller_stuff/ | egrep ^grp_${1}_fill\.B[0-9]+\.pdb$`
     if [ `echo $files|wc -c` = 1 ]; then
             echo "ERROR: modeller files not found"
             exit 1
     fi

filename="stage2.m"
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
   reducer_fill_H -FLIP modeller_stuff/$i > ../output_pdb/grp_$1_${c}_H.pdb
  echo "cd .." >> $filename
  echo "addpath(genpath(pwd));" >> $filename
  echo "cd code;" >> $filename
   echo "load('../protein/$1/stage1.mat')" >> $filename
   echo "OUTPUT.modeller_n_H  = '../output_pdb/grp_$1_${c}_H.pdb';" >> $filename
   echo "OUTPUT.model_count = ${c};" >> $filename
   echo "doAnchoredLoc;" >> $filename
   echo "save('../protein/$1/stage2_${c}.mat')" >> $filename
   echo "%========================================================" >> $filename
   echo " " >> $filename
   echo " "
done

# ------------ run stage -2 script in matlab ---------------------
   matlab -nosplash -nodesktop < $filename
# ---------- create scripts for anchored localization ----------------
#grp_5_ancloc_refine.pdb filename changed to demo_5_ancloc_refine.pdb
files=`ls ../output_pdb/ | egrep ^$1_[0-9]+_ancloc_refine.pdb`
     if [ `echo $files|wc -c` = 1 ]; then
             echo "ERROR: anchored localization files not found"
             exit 1
     fi

c=0
for i in `echo $files`
do
    c=$((c+1))
    cp ../output_pdb/$i gromacs_stuff/.
            if [ -f gromacs_stuff/gromacs_script_$c.txt ];then
                  echo " DELETING gromacs_stuff/gromacs_script_$c.txt"
                  rm -f gromacs_stuff/gromacs_script_$c.txt
            fi
    echo "./sa_md_template.sh $1_$c"
    #sed "s/grp_5/$1_$c/g" gromacs_stuff/steps_taken_gromacs_template.txt >> gromacs_stuff/gromacs_script_$c.txt
done

#-------------------- run gromacs script -------------------




