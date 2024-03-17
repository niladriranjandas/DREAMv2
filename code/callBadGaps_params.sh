#!/bin/bash

#
#Copyright (c) 2024 niladriranjandas
#call the gap correct for "normal mode" for DREAMv2

protname="$1"

alreadyrun="../proteinsparams/listoffiles.txt"
em_folder="md_related_v2"

protfolder="$protname"
# ------------- #
for pdbfile in `ls $em_folder/$protfolder/*/*em_nosol.pdb`
do
    proteinname=`echo $protfolder`
	
	PROTFOLDER="../protein"
	paramfile=$PROTFOLDER"/"$proteinname"/"$proteinname".xml"
	seq=`xmllint --xpath "string(//seq_file)" "$paramfile"`
	stage1data=../protein/"$protfolder"/stage1_"$protfolder".mat
	foldercurr=`dirname $pdbfile`
	uplfile=`ls $foldercurr/*concatupl.upl`
	
	# assumed that makeXplorUplAco.sh has already ran
        xplornoe_=`xmllint --xpath "string(//upl_xplor)" "$paramfile"`
	xplordihed_=`xmllint --xpath "string(//aco_xplor)" "$paramfile"`
	xplornoe=`ls ../protein/$proteinname/$xplornoe_`
	if [ ! -z "$xplordihed_" ]; then
		xplordihed=`ls ../protein/$proteinname/$xplordihed_`
	else
		xplordihed=""
	fi
       
	if [ ! -f ../protein/$proteinname/$seq ]; then
		echo "$seq file not found"
		exit -1
	fi
	if [ ! -f $stage1data ]; then
		echo "$stage1data file not found"
		exit -1
	fi
	if [ ! -f $uplfile ]; then
		echo "$uplfile file not found"
		exit -1
	fi
	if [ ! -f $xplornoe ]; then
		echo "$xplornoe file not found"
		exit -1
	fi
	if [ ! -z "$xplordihed_" ]; then
		if [ ! -f $xplordihed ]; then
			echo "$xplordihed file not found"
			exit -1
		fi
	fi
        
       name_for_gap=`basename $pdbfile`
	name_for_gap_=${name_for_gap::-4}	
	#gap_correct/correctBadGaps_old.sh "$name_for_gap_" "$pdbfile" "../protein/$proteinname/$seq" "$stage1data" "$uplfile" "$xplornoe" "$xplordihed" &
	#echo " "
	#echo    "gap_correct/correctBadGaps.sh $name_for_gap_ $pdbfile ../protein/$proteinname/$seq $stage1data $uplfile $xplornoe $xplordihed"
	gap_correct/correctBadGaps.sh "$name_for_gap_" "$pdbfile" "../protein/$proteinname/$seq" "$stage1data" "$uplfile" "$xplornoe" "$xplordihed" &
	#if [[ -z $xplordihed ]]; then
        #       gap_correct/correctBadGaps.sh "$name_for_gap_" "$pdbfile" "../protein/$proteinname/$seq" "$stage1data" "../protein/$proteinname/$uplfile" "../protein/$proteinname/$xplornoe" &
	#else
	#       echo "gap_correct/correctBadGaps.sh $name_for_gap_ $pdbfile ../protein/$proteinname/$seq $stage1data ../protein/$proteinname/$uplfile ../protein/$proteinname/$xplornoe ../protein/$proteinname/$xplordihed "	       
	       #gap_correct/correctBadGaps.sh "$name_for_gap_" "$pdbfile" "../protein/$proteinname/$seq" "$stage1data" "../protein/$proteinname/$uplfile" "../protein/$proteinname/$xplornoe" "../protein/$proteinname/$xplordihed" &
	#fi
done
wait	
# ------------- #

