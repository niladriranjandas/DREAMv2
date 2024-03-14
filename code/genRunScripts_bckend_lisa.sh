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

 # call a small function to check all params are present

 # -----stage-I master file 
  filename="stage1_$1.m"
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
  #echo "cd /home/niladri/Documents/matlab_extras/cvx/" >> $filename
  echo "cd ../packages/cvx/" >> $filename
  echo "cvx_setup" >> $filename     #included to run in nmr machines
  echo "cvx_startup" >> $filename
  echo "cd (curr_path)" >> $filename
 echo "consolidate" >> $filename

 # ------- make scripts for stage-2 initilization ------
 echo "%-----prepare scripts for modeller --------" >> $filename
 echo "prepForFillGap" >> $filename

 #-------------- save the workspace ---------------
 echo "%--------- save workspace -----------------" >> $filename
 echo "save('../protein/$1/stage1_$1.mat')" >> $filename
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

 
 
 if [ ! -f modeller_stuff/alignment_"$1".ali ]; then
     echo "Alignment file not generated for modeller in modeller_stuff/alignment_$1.ali"
     exit 1
 fi

 if [ ! -f modeller_stuff/do_loop_dntmove_"$1".py ]; then
      echo "Python script not generated for modeller in modeller_stuff/do_loop_dntmove_$1.py" 
      exit 1
 fi

# ---------- do Modeller run ---------------------------------
 echo "==============Running Modeller==========================="
 cd modeller_stuff
 #mod9.19 do_loop_dntmove.py
 mod9.23 do_loop_dntmove_"$1".py
 cd ..

# --------- gather the structures from modeller run reducer ----
files=`ls modeller_stuff/ | egrep ^grp_${1}_fill\.B[0-9]+\.pdb$`
     if [ `echo $files|wc -c` = 1 ]; then
             echo "ERROR: modeller files not found"
             exit 1
     fi

filename="stage2_$1.m"
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
  # reducer_fill_H -FLIP modeller_stuff/$i > ../output_pdb/grp_$1_${c}_H.pdb
  ../packages/Reduce_fill_H_protein/reduce.3.23.130521.linuxi386/reducer_fill_H -FLIP  modeller_stuff/$i > ../output_pdb/grp_$1_${c}_H.pdb
  echo "cd .." >> $filename
  echo "addpath(genpath(pwd));" >> $filename
  echo "cd code;" >> $filename
   echo "load('../protein/$1/stage1_$1.mat')" >> $filename
   echo "OUTPUT.modeller_n_H  = '../output_pdb/grp_$1_${c}_H.pdb';" >> $filename
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

  ./rungenEMscripts.sh $1 ../output_pdb $new_file | tail -1 > $em_file 
  em_run_cmd=`./rungenEMscripts.sh $1 ../output_pdb $new_file | tail -1`
  ./$em_run_cmd
# -----------------------------------------------------------------------

## -- move the folders -- ##
  #- output_pdbs
  mkdir ../protein/$1/output_pdb
  cp ../output_pdb/*$1* ../protein/$1/output_pdb/.
  #- output_plots
  mkdir ../protein/$1/output_plots
  cp ../output_pdb/*$1* ../protein/$1/output_plots/.
  #- logs
  mkdir ../protein/$1/log
  cp ../log/*$1* ../protein/$1/log/.

  #-md_related
  mkdir ../protein/$1/md_related
  cp -r md_related/$1/* ../protein/$1/md_related/.
  #-move temp dssp files
  mkdir ../protein/$1/tmp_dssp_files
  mv *$1*dssp.txt ../protein/$1/tmp_dssp_files/.
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

