#!/bin/bash

param_file=../protein/$1/$1.xml

#---- check if all files are in order ---

  if [ -f ../protein/$1/$1.xml ]; then
       echo "Found PARAM file" 
  else
       echo "ERROR: PARAM file not found"
       exit 1
  fi

###########################################

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

  #if [ -f $nogapfile ]; then
  if [ ! -z $nogapfile ]; then
          ./rungenEMscripts_nogap.sh $1 ../output_pdb $new_file | tail -1 > $em_file
           em_run_cmd=`./rungenEMscripts_nogap.sh $1 ../output_pdb $new_file | tail -1`
	   #echo "rungenEMscripts_nogap.sh $1 ../output_pdb $new_file | tail -1 > $em_file"
  else
          ./rungenEMscripts.sh $1 ../output_pdb $new_file | tail -1 > $em_file
          em_run_cmd=`./rungenEMscripts.sh $1 ../output_pdb $new_file | tail -1`
	  #echo "./rungenEMscripts.sh $1 ../output_pdb $new_file | tail -1 > $em_file"
  fi
  echo "$em_run_cmd"
  ./$em_run_cmd
# -----------------------------------------------------------------------

## -- move the folders -- ##
  #- output_pdbs
  mkdir ../protein/$1/output_pdb
  cp ../output_pdb/*$1* ../protein/$1/output_pdb/.
  #- output_plots
  mkdir ../protein/$1/output_plots
  cp ../output_plots/*$1* ../protein/$1/output_plots/.
  #- logs
  mkdir ../protein/$1/log
  cp ../log/*$1* ../protein/$1/log/.

  #- modeller stuff
  mkdir ../protein/$1/modeller_stuff
  cp modeller_stuff/*$1* ../protein/$1/modeller_stuff/.
  #-md_related
  mkdir ../protein/$1/md_related
  cp -r md_related/$1/* ../protein/$1/md_related/.
  #-move temp dssp files
  mkdir ../protein/$1/tmp_dssp_files
  mv *$1*dssp.txt ../protein/$1/tmp_dssp_files/.
	
