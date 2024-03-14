#!/bin/bash

BASE_DIR="/home/niladri/Documents/Disco_etc_all_in_1/our_algo"

for i in `ls "$BASE_DIR/protein"`
do
	z=`grep "<size_cutoff>0.3</size_cutoff>" "$BASE_DIR/protein/$i/$i.xml"`
	if [ ! -z $z ]; then
		#--- replace the config file 
echo
echo "-----------------------------------------------------------------------------------------"
		sed -i 's/<size_cutoff>0.3/<size_cutoff>0.1/' "$BASE_DIR/protein/$i/$i.xml"

		#--- rename the archive folder
		mv "$BASE_DIR/archive/$i" "$BASE_DIR/archive/$i_oldconfig"
	fi
done


