#!/bin/bash

folder="$1"
destfolder="../restricted_em"

IFS=',' read -r -a array <<< "$folder"

for ele in "${array[@]}"
do
	if [ -d "../protein/$ele" ]; then
		emfolder="../protein/$ele/md_related"
		if [ -d "$emfolder" ]; then
			cp -R $emfolder $destfolder/$ele
			curr=`pwd`
			cd "$destfolder"
			../restricted_em/runEMforProt.sh "$ele"
			cd "$curr"
		else
			echo "ERROR: $emfolder not found"
		fi
	else
		echo "ERROR: ../protein/$ele not found"
	fi
done
