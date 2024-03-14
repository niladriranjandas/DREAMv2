#!/bin/bash
protname="$1"

alreadyrun="../proteinsparams/listoffiles.txt"
em_folder="md_related_v2"

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
	protfolder=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`
	for pdbfile in `ls $em_folder/$protfolder/*/*em_nosol.pdb`
	do
	        proteinname=`echo $protfolder`
		
		PROTFOLDER="../protein"
		paramfile=$PROTFOLDER"/"$proteinname"/"$proteinname".xml"
		seq=`xmllint --xpath "string(//seq_file)" "$paramfile"`

		#stage1data=../protein/"$protfolder"/stage1_"$protfolder".mat
		foldercurr=`dirname $pdbfile`
		#uplfile=`ls $foldercurr/*concatupl.upl`
		
		# assumed that makeXplorUplAco.sh has already ran
		xplornoe=`ls ../protein/$proteinname/*noe.tbl`
		xplordihed=`ls ../protein/$proteinname/*aco.tbl`

		if [ ! -f ../protein/$proteinname/$seq ]; then
			echo "$seq file not found"
			exit -1
		fi
		#if [ ! -f $stage1data ]; then
		#	echo "$stage1data file not found"
		#	exit -1
		#fi
		#if [ ! -f $uplfile ]; then
		#	echo "$uplfile file not found"
		#	exit -1
		#fi
		#if [ ! -f $xplornoe ]; then
		#	echo "$xplornoe file not found"
		#	exit -1
		#fi
		#if [ ! -f $xplordihed ]; then
		#	echo "$xplordihed file not found"
		#	exit -1
		#fi
	        
                name_for_gap=`basename $pdbfile`
		name_for_gap_=${name_for_gap::-4}	
		#gap_correct/correctBadGaps.sh "$name_for_gap_" "$pdbfile" "../protein/$proteinname/$seq" "$stage1data" "$uplfile" "$xplornoe" "$xplordihed"
		echo " "
		echo "gap_correct/endBapResiCorrect.sh $name_for_gap_ ../protein/$proteinname/$seq $pdbfile"
	done	
else
	# -- gather the pdbs for the files -- #
	MAXPRINT=5
	RESIDUAL_path=gap_correct/RESIDUAL_FILES
        proteinname=`echo $protname`
        mkdir $RESIDUAL_path/$protname
		
	PROTFOLDER="../protein"
	paramfile=$PROTFOLDER"/"$proteinname"/"$proteinname".xml"
	
	seq=`xmllint --xpath "string(//seq_file)" "$paramfile"`

	file_to_read=$protname"_findwhichtorun.txt"
	folder_stage1=`head -1 $file_to_read`  # any one would do
	#stage1data=../protein/"$folder_stage1"/stage1_"$folder_stage1".mat
	#uplfile=`ls ../protein/$protname/*concatupl.upl`   #check ones
		
	# assumed that makeXplorUplAco.sh has already ran
	#xplornoe=`ls ../protein/$proteinname/*noe.tbl`
	#xplordihed=`ls ../protein/$proteinname/*aco.tbl`
	
	##### choose which pdb's tp run #####
	while IFS= read -r line || [[ -n "$line" ]]; do
		for pdbfile_i in `ls $em_folder/$line/*/*em_nosol.pdb`
		do	        
			allfiles=$allfiles",'"$pdbfile_i"'"
		done
	done < "$file_to_read"

	allfiles_=${allfiles:1}
	pdbfilesall="{"$allfiles_"}"

	choosefilemat=$protname"_whichwaterref_mat.m"
	oppfiletxt=$protname"_whichwaterref_mat.txt"

	echo "cd .." > $choosefilemat
	echo "addpath(genpath(pwd));" >> $choosefilemat
	echo "cd code;" >> $choosefilemat
	echo "pdbfilesall=$pdbfilesall;" >> $choosefilemat
	echo "gapstruct = findWhichToChoose(pdbfilesall,'$oppfiletxt')" >> $choosefilemat

	matlab -nodesktop -nosplash < $choosefilemat

	#####################################
	printcount=0	
	while IFS= read -r line || [[ -n "$line" ]]; do
		printcount=$((printcount+1))
		pdbfile_i=$line
		# call gap_correct/correctBadGaps with $pdbfile_i
		if [ ! -f ../protein/$proteinname/$seq ]; then
			echo "$seq file not found"
			exit -1
		fi
		#if [ ! -f $stage1data ]; then
		#	echo "$stage1data file not found"
		#	exit -1
		#fi
		#if [ ! -f $uplfile ]; then
		#	echo "$uplfile file not found"
		#	exit -1
		#fi
		#if [ ! -f $xplornoe ]; then
		#	echo "$xplornoe file not found"
		#	exit -1
		#fi
		#if [ ! -f $xplordihed ]; then
		#	echo "$xplordihed file not found"
		#	exit -1
		#fi
	        
                name_for_gap=`basename $pdbfile_i`
		name_for_gap_=${name_for_gap::-4}
		if [[ $printcount -le $MAXPRINT ]]; then 
			#gap_correct/correctBadGaps.sh "$name_for_gap_" "$pdbfile_i" "../protein/$proteinname/$seq" "$stage1data" "$uplfile" "$xplornoe" "$xplordihed"
			echo " "
			echo "gap_correct/endBapResiCorrect.sh $name_for_gap_ ../protein/$proteinname/$seq $pdbfile"
		fi
		mv $name_for_gap_"_dssp.txt" $RESIDUAL_path/$protname/.
	done < "$oppfiletxt"
	#done < "$file_to_read"
		
	mv $choosefilemat $RESIDUAL_path/$protname/.
	mv $oppfiletxt $RESIDUAL_path/$protname/.
fi
