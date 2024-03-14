#!/bin/bash

#usage: ./drawRamaFromPDB.sh ../../../archive/1d8v/gromacs_stuff/md_0_1d8v_1_nosol.pdb rama_nosol_1d8v.png

if [ $# != 2 ]; then
    echo "Usage: ./drawRamaFromPDB.sh <full_path_pdb_file> <rama_plot.png>"
    exit 1;
fi

#-------- write the torsional angles files -------------------------------

./pdbtorsions $1 > test.tsv
        if [ $? -ne 0 ]; then
                echo "Error in pdbtorsions"
                exit 1
        fi

#-------- fill in labels Pre-pro -----------------------------------------
lines_pro=`sed -n '/Proline/=' test.tsv`
pattern_sed=""
for i in `echo $lines_pro`
do
   lines_prepro=`expr $i - 1`
   if [[ ! $(echo $lines_pro | grep $lines_prepro) ]]; then
      pattern_sed="${pattern_sed};${lines_prepro}s/General/Pre-Pro/"
   fi
done
pattern_sed_=`echo $pattern_sed | cut -c 2-`  # avoid the beginning semicolon

sed -i "$pattern_sed_" test.tsv

#sed '25s/General/Pre-pro/' test.tsv
#sed '25s/General/Pre-pro/;34s/General/Pre-pro/' test.tsv

Rscript draw_rama.r $2 test.tsv
        if [ $? -eq 0 ]; then
                rm -f test.tsv
        else    
                exit 1
        fi  
