#!/bin/bash

#---- little script to list the files created -----

srchPatternDisp () {
    patter=$1
    folder=$2

    if [ ! -d $2 ]; then
         echo "NO folder found: $2"
         return 1
    fi   
    
    files=`ls $2 | egrep "$1"`
    [[ -z $files ]] && return 1
    
    printf 'Folder: %s\n Files:' $2
    for i in `echo $files`
    do
      printf '\t %s\n' $i
    done

  printf '\n'
  return 0
}

if [ $# != 1 ]; then
    echo "Usage: ./listOutputs.sh <protein_name>"
    exit 1;
fi
  
# ---- module: div_conquer.m ------
pattern="gr_$1_[0-9]+.pdb"
folder="../output_pdb"
echo "__________MODULE: div_conquer_____________________________________"
srchPatternDisp $pattern $folder
[[ $? == 1 ]] && echo "!No files found"

# ---- module: consolidate.m ------
pattern="grp_registration_$1.pdb"
folder="../output_pdb"
echo "__________MODULE: consolidate: pdb's _____________________________"
srchPatternDisp $pattern $folder
[[ $? == 1 ]] && echo "!No files found"

# ---- module: consolidate.m ------
pattern="$1_fragment_graph.eps|$1_ramachandran_run[0-9]+.fig"
folder="../output_plots"
echo "__________MODULE: consolidate: plots _____________________________"
srchPatternDisp $pattern $folder
[[ $? == 1 ]] && echo "!No files found"

# ---- module: prepForFillGap.m ------
pattern="grp_$1_[0-9]+_H.pdb"
folder="../output_pdb"
echo "__________MODULE: prepForFillGap: pdb's __________________________"
srchPatternDisp $pattern $folder
[[ $? == 1 ]] && echo "!No files found"

# ---- module: prepForFillGap.m ------
pattern=".ali$|.py$"
folder="modeller_stuff"
echo "__________MODULE: prepForFillGap: modeller scripts _______________"
srchPatternDisp $pattern $folder
[[ $? == 1 ]] && echo "!No files found"

# ---- module: doAnchoredLoc.m ------
pattern="$1_[0-9]+_ancloc.pdb|$1_[0-9]+_ancloc_refine.pdb"
folder="../output_pdb"
echo "__________MODULE: doAnchoredLoc: pdb's ___________________________"
srchPatternDisp $pattern $folder
[[ $? == 1 ]] && echo "!No files found"

# ---- module: doAnchoredLoc.m ------
pattern="$1_rama_anchor_noref_[0-9]+.fig|$1_rama_anchor_ref_[0-9]+.fig"
folder="../output_plots"
echo "__________MODULE: doAnchoredLoc: plots ___________________________"
srchPatternDisp $pattern $folder
[[ $? == 1 ]] && echo "!No files found"

# ----- module: genRunScripts.sh -----
pattern="^md|gromacs_script_[0-9]+.txt"
folder="gromacs_stuff"
echo "__________MODULE: genRunScripts.sh: gromacs scripts ______________"
srchPatternDisp $pattern $folder
[[ $? == 1 ]] && echo "!No files found"

