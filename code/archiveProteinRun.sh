#!/bin/bash

# -------defining small function to coppy files
coppyFiles () {
       files=`ls $1 | egrep $2`
       if [ -z "$files" ]; then
          echo "error: No files containing $2 found in $1"   
          return 1
       else          
          find $1 -name "*$2*" | xargs cp -t $3/.
          [[ $? == 0 ]] && return 0 || return 1
       fi
}

#---check if everything is proper in the command line argument supplied
if [ $# != 2 ]; then
   echo "Error: no arguments provided. Usage: ./archiveProteinRun <protein-name> <dest-folder>"
   exit 1
fi

if [ ! -d $2 ]; then
    echo "Error: destination folder for archive doesn't exist"
    exit 1
fi

if [ ! -d "../protein/$1" ]; then
     echo "Error: protein folder abscent in ../protein/$1 for archiving"
     exit 1
fi

#---------- create a folder by the name of protein ----------------------
echo "...creating $2/$1"
mkdir "$2/$1"
[[ $? == 0  ]] && echo Successfull || echo Error

#-------------archive the protein folder---------------------------------
echo "...archiving folder ../protein/$1/."
cp -r "../protein/$1/." "$2/$1/."
[[ $? == 0  ]] && echo Successfull || echo Error

#-------------archive output_pdb-----------------------------------------
echo "...archiving folder ../protein/output_pdb/."
mkdir "$2/$1/output_pdb"
[[ $? == 0  ]] && echo Successfull || echo Error
echo "Copy files from output_pdb"
coppyFiles ../output_pdb $1 $2/$1/output_pdb
[[ $?  == 0 ]] && echo successfull || echo error

# ------------archive  modeller_stuff ....................................
echo "...archiving folder modeller_stuff/."
mkdir "$2/$1/modeller_stuff"
[[ $? == 0  ]] && echo Successfull || echo Error
echo "Copy files from modeller_stuff..."
coppyFiles modeller_stuff $1 $2/$1/modeller_stuff
[[ $? == 0 ]] && echo successfull || echo error

# ------------archive gromacs_stuff ---------------------------------------
echo "...archiving folder gromacs_stuff/."
mkdir "$2/$1/gromacs_stuff"
[[ $? == 0  ]] && echo Successfull || echo Error
echo "Copy files from gromacs_stuff..."
cd gromacs_stuff
ls | egrep "^md|gromacs_script_[0-9]+.txt" | xargs cp -t ../$2/$1/gromacs_stuff
[[ $? == 0 ]] && echo successfull || echo error
cd ..

# ------------archive output_plots ----------------------------------------
echo "...archiving folder ../protein/output_plots/."
mkdir "$2/$1/output_plots"
[[ $? == 0  ]] && echo Successfull || echo Error
echo "Copy files from output_plots"
coppyFiles ../output_plots $1 $2/$1/output_plots
[[ $? == 0 ]] && echo successfull || echo error

#------------archive log files -------------------------------------------
echo "...archiving folder ../log/."
mkdir "$2/$1/log"
[[ $? == 0  ]] && echo Successfull || echo Error
echo "Copy files from log"
cp ../log/*"$1"* $2/$1/log/.

#-----------archive stage1.m and stage2.m---------------------------------
echo "...archiving stage1.m and stage2.m"
cp stage1_"$1".m $2/$1/.
[[ $? == 0  ]] && echo Successfull || echo Error
cp stage2_"$1".m $2/$1/.
[[ $? == 0  ]] && echo Successfull || echo Error

