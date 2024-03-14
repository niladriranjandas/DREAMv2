#!/bin/bash

# Input: 1. uplfile
#	 2. suffix
#	 3. TOLSTRONG
#	 4. TOLMEDIUM
#        5. TOLLOW
#	 e.g.: ambiguous_bounds/masterPreCheckUpl.sh -suf 2m4k_ambi5r1_noe -upl 2m4k/2m4k_5perc_1_added.upl -seq 2m4k/2m4k.seq -hb '' -aco 2m4k/2m4k_concat_dihed.aco


help()
{
    echo "Usage: masterPreCheckUpl.sh [ -suf | --suffix  ]
               [ -upl | --uplfile ]
               [ -seq | --seqfile ]
               [ -hb  | --hbondfile ]
               [ -aco | --acofile ]
               [ -h | --help  ]"
    exit 2
}


# set up argparse

SHORT=suf:,upl:,seq:,hb:,aco:,str:,med:,low:,h
LONG=suffix:,uplfile:,seqfile:,hbondfile:,acofile:,help
OPTS=$(getopt -a -n masterPreCheckUpl --options $SHORT --longoptions $LONG -- "$@")

VALID_ARGUMENTS=$# # Returns the count of arguments that are in short or long options

if [ "$VALID_ARGUMENTS" -eq 0 ]; then
  help
fi

eval set -- "$OPTS"

cd ambiguous_bounds

while :
do
  case "$1" in
    -suf | --suffix )
      suf="$2"
      shift 2
      ;;  
    -upl | --uplfile )
      upl="$2"
      shift 2
      ;;
    -seq | --seqfile )
      seq="$2"
      shift 2
      ;;
    -hb | --hbondfile )
      hb="$2"
      shift 2
      ;;
    -aco | --acofile )
      aco="$2"
      shift 2
      ;;
    -h | --help)
      help
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      help
      ;;
  esac
done

# call neibhorhood upl check
  str=5
  med=5
  low=2
  suf1=$suf"_1"
  python guessAmbiByBounds.py $upl $suf1 $str $med $low
  strong_file1=$suf1"_strong_flag.csv"
  medium_file1=$suf1"_medium_flag.csv"
  low_file1=$suf1"_low_flag.csv"

  str=10
  med=10
  low=5
  suf2=$suf"_2"
  python guessAmbiByBounds.py $upl $suf2 $str $med $low
  strong_file2=$suf2"_strong_flag.csv"
  medium_file2=$suf2"_medium_flag.csv"
  low_file2=$suf2"_low_flag.csv"

  str=15
  med=15
  low=10
  suf3=$suf"_3"
  python guessAmbiByBounds.py $upl $suf3 $str $med $low
  strong_file3=$suf3"_strong_flag.csv"
  medium_file3=$suf3"_medium_flag.csv"
  low_file3=$suf3"_low_flag.csv"

cd ..
 
# call the matlab module for subgraph thing
  seq_nameonly=`echo $seq | awk -F'/' '{print $NF}'`
  upl_nameonly=`echo $upl | awk -F'/' '{print $NF}'`
  aco_nameonly=`echo $aco | awk -F'/' '{print $NF}'`
  filename_density=$suf"_density_check.m"
  echo "cd .." > $filename_density
  echo "addpath(genpath(pwd));" >> $filename_density
  echo "cd code;" >> $filename_density
  echo "%------density check-------" >> $filename_density
  echo "[upl_densities, upl_densities_strong, upl_densities_medium, upl_densities_weak, outfile] = subGraphDensityFromUpl('$suf','$seq_nameonly', {'$upl_nameonly'}, '', '$aco_nameonly');" >> $filename_density
  echo "exit" >> $filename_density
  
# ----------- run matlab for stage-I --------------------
#exit 1
 matlab -nodesktop -nosplash < $filename_density
 
 file_m_strong=$suf"_strong.csv"
 file_m_medium=$suf"_medium.csv"
 file_m_weak=$suf"_weak.csv"

# ----------- assimilate result from the above two criteria
echo " "
density_strong_lines=`wc -l ambiguous_bounds/$file_m_strong | awk '{print $1}'`
density_medium_lines=`wc -l ambiguous_bounds/$file_m_medium | awk '{print $1}'`
density_weak_lines=`wc -l ambiguous_bounds/$file_m_weak | awk '{print $1}'`

if [[ "$density_strong_lines" -ne 0 ]]; then
    strong_py_files="$strong_file1,$strong_file2,$strong_file3"
    oppfilename=$suf"_strong_overall.csv"
    cd ambiguous_bounds
    python parseAllFiles.py $file_m_strong $strong_py_files $oppfilename
    cd ..
else
    echo "EMPTY strong"
fi

if [[ "$density_medium_lines" -ne 0 ]]; then
    medium_py_files="$medium_file1,$medium_file2,$medium_file3"
    oppfilename=$suf"_medium_overall.csv"    
    cd ambiguous_bounds
    echo "python parseAllFiles.py $file_m_medium $medium_py_files $oppfilename"
    python parseAllFiles.py $file_m_medium $medium_py_files $oppfilename    
    cd ..
else
    echo "EMPTY medium"
fi


if [[ "$density_weak_lines" -ne 0 ]]; then
    low_py_files="$low_file1,$low_file2,$low_file3"
    oppfilename=$suf"_weak_overall.csv"    
    cd ambiguous_bounds
    python parseAllFiles.py $file_m_weak $low_py_files $oppfilename    
    cd ..
else
    echo "EMPTY weak"
fi

# -----------------------------------------------  #
raw_upl_=$upl".raw_upl.upl"
raw_upl=`echo $raw_upl_  | awk -F'/' '{print $NF}'`
#$suf"_strong_overall.csv"
#$suf"_medium_overall.csv"
#$suf"_weak_overall.csv"

mv $raw_upl ambiguous_bounds/.
cd ambiguous_bounds
if [[ -f $suf"_strong_overall.csv" && -f $suf"_medium_overall.csv" && -f $suf"_weak_overall.csv" ]]; then
      echo filteredUpls.py $raw_upl $suf"_strong_overall.csv" $suf"_medium_overall.csv" $suf"_weak_overall.csv" $suf"_noambi.upl"
      python filteredUpls.py $raw_upl $suf"_strong_overall.csv" $suf"_medium_overall.csv" $suf"_weak_overall.csv" $suf"_noambi.upl"
elif [[ -f $suf"_medium_overall.csv" && -f $suf"_weak_overall.csv" ]]; then
      echo filteredUpls.py $raw_upl '' $suf"_medium_overall.csv" $suf"_weak_overall.csv" $suf"_noambi.upl" 
      python filteredUpls.py $raw_upl '' $suf"_medium_overall.csv" $suf"_weak_overall.csv" $suf"_noambi.upl" 
elif [[ -f $suf"_weak_overall.csv" ]]; then
      echo filteredUpls.py $raw_upl '' '' $suf"_weak_overall.csv" $suf"_noambi.upl"
      python filteredUpls.py $raw_upl '' '' $suf"_weak_overall.csv" $suf"_noambi.upl"
else
      echo "ERROR: None of flagged upls files were generated"
fi
cd ..
# ----------- move all the files ------------------#
mkdir ../protein/$suf/ambiguous_bounds
mv ambiguous_bounds/$suf* ../protein/$suf/ambiguous_bounds/.

mv $suf* ../protein/$suf/ambiguous_bounds/.

