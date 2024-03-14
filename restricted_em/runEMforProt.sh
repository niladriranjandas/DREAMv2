#!/bin/bash

folder="$1"
	if [ -d "$folder" ]; then
		for i in `ls -d $folder/*/`
		do
			cp getPDBFromBMRB_edit.sh $i.
			cp pseudo.csv $i.
			cp distResItpFromUpl.sh $i.
			cp nomanclature_edit.csv $i.

	                runfile=`find "$i" -iname "*$folder*sh"`
			echo "cd $i"
			echo "./$runfile"
		done
	else
		echo "ERROR: $folder not found"
	fi
