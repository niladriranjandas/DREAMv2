#!/bin/bash

protname="$1"

RAMA_PATH="../packages/draw_rama/"
RAMA_CMD="./drawRama.sh"

if [ -z $protname ]; then
	echo "arg-1 needed"
	exit -1
fi

alreadyrun="../proteinsparams/listoffiles.txt"
protpath="../protein"
foldername="finalpdbfiles"

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
	protfolder=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`

	waterref_folder=$protpath"/"$protfolder"/"$foldername
	echo $waterref_folder
	for pdbs in `ls $waterref_folder/*pdb`
	do
		curr_path=`pwd`
		full_path=$curr_path"/"$pdbs
		
		filename=`basename $pdbs`
		foldernamelocal=`dirname $pdbs`		
		
		cd $RAMA_PATH
		echo "$RAMA_CMD $full_path"
		$RAMA_CMD $full_path
		mv $filename* $curr_path"/"$foldernamelocal/.
		cd $curr_path
	done
else
	#echo "code it"
	protfolder=$protname"_allparam"
	waterref_folder=$protpath"/"$protfolder"/"$foldername
	echo $waterref_folder
	for pdbs in `ls $waterref_folder/*pdb`
	do
		curr_path=`pwd`
		full_path=$curr_path"/"$pdbs
		
		filename=`basename $pdbs`
		foldernamelocal=`dirname $pdbs`		
		
		cd $RAMA_PATH
		echo "$RAMA_CMD $full_path"
		$RAMA_CMD $full_path
		mv $filename* $curr_path"/"$foldernamelocal/.
		cd $curr_path
	done	
fi
