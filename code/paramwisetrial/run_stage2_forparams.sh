#!/bin/bash

line="$1"
filetmp="$2"
logfile="$3"

path_for_log=`pwd`
LOG_MASTER=$path_for_log/$line"_minimal_log.txt"

# ---------------------------------------------------------
#while IFS= read -r line || [[ -n "$line" ]]; do
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
         #mod10.0 do_loop_dntmove_"$line".py
 	 $MOD do_loop_dntmove_"$line".py
	 echo "running modeller $line:$?" >> $LOG_MASTER
	 cd ..

	 # --------- gather the structures from modeller run reducer ----
	files=`ls modeller_stuff/ | egrep ^grp_${line}_fill\.B[0-9]+\.pdb$`
	if [ `echo $files|wc -c` = 1 ]; then
             echo "ERROR: modeller files not found"
             exit 1
	fi
	echo $files >> $logfile  #cat $files >> $logfile

	#filename="stage2_$line.m"
        #if [ -f $filename ];then
      	#     echo "DELETING $filename"
        #     rm -f $filename
        #fi
	#echo "%% ---------------- Do anchored localization -----------" > $filename
	c=0
	stage2files=""
	for i in `echo $files`
	do
	   c=$((c+1))
	   filename="stage2_"$line"_"$c".m"
	   stage2files=$filename","$stage2files
	   if [ -f $filename ]; then
		echo "DELETING $filename"
	        rm -f $filename
	   fi
 
	   #echo $i
	   #wc -c $i
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
	   # ------------ run stage -2 script in matlab ---------------------
	   printf "\n \t\t\t anchored locatlization %s" "$filename" >> $logfile
	   matlab -nosplash -nodesktop < $filename &
	   echo "anchored location for $filename:$?" >> $LOG_MASTER
	done
	wait
      
        stage2files=${stage2files::-1}

	## ------------ run stage -2 script in matlab ---------------------
	#   printf "\n \t\t\t anchored locatlization %s" "$filename" >> $logfile
	#   matlab -nosplash -nodesktop < $filename
	## ------------- 
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
	  echo "ending anchored localization for $line:$?" >> $LOG_MASTER
	# -----------------------------------------------------------------------
	## -- move the folders -- ##
	  #- output_pdbs
	  mkdir ../protein/$line/output_pdb
	  #cp ../output_pdb/*$line* ../protein/$line/output_pdb/.
          mv ../output_pdb/*$line* ../protein/$line/output_pdb/.
	  #- output_plots  this is currently not requried
	  #mkdir ../protein/$line/output_plots
	  #cp ../output_plots/*$line* ../protein/$line/output_plots/.
	  #mv ../output_plots/*$line* ../protein/$line/output_plots/.
	  #- logs
	  mkdir ../protein/$line/log
	  #cp ../ilog/*$line* ../protein/$line/log/.
          mv ../log/*$line* ../protein/$line/log/.
	  #- modeller_stuff						copied later on
	  #mkdir ../protein/$line/modeller_stuff			copied later on
	  #cp modeller_stuff/*$line* ../protein/$line/modeller_stuff/.	copied later on
	  #mv modeller_stuff/*$line* ../protein/$line/modeller_stuff/.	copied later on
	  #-md_related
	  #mkdir ../protein/$line/md_related
	  #cp -r md_related/$line/* ../protein/$line/md_related/.
	  mkdir ../protein/$line/md_related_v2
	  #cp -r md_related_v2/$line/* ../protein/$line/md_related_v2/.
	  mv md_related_v2/$line/* ../protein/$line/md_related_v2/.
	  #-move temp dssp files
	  mkdir ../protein/$line/tmp_dssp_files
	  mv *$line*dssp.txt ../protein/$line/tmp_dssp_files/.

	  # - move stage2 files 
	  #cp $filename ../protein/$line/.
	  #mv $filename ../protein/$line/.
          IFS=',' read -r -a stage2array <<< "$stage2files"
          for mfiles in "${stage2array[@]}"
	  do
		mv $mfiles ../protein/$line/.
	  done	  

	  # - move modeller files
	  mkdir ../protein/$line/modeller
	  #cp modeller_stuff/*$line* ../protein/$line/modeller/.
	  mv modeller_stuff/*$line* ../protein/$line/modeller/.
#done < "$filetmp"

