#!/bin/bash

# run this first export LD_LIBRARY_PATH=/data2/nmr/gcc/lib:/data2/nmr/gcc/lib64:$LD_LIBRARY_PATH
protname="$1"

alreadyrun="../proteinsparams/listoffiles.txt"

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
	paramfile=`egrep ^$protname $alreadyrun | awk -F',' '{print $2}'`
        if [ -d "../protein/$paramfile" ]; then
		echo "../protein/$paramfile already exists. Backing up."
		ext=`date | tr ' ' '_' | tr ':' '_'`
		tmp=$paramfile"_"$ext
		mv ../protein/"$paramfile" ../protein/"$tmp" 
                #echo "cp -R \"../protein/$paramfile ../protein/$tmp\""
	fi
	cp -R ../proteinsparams/$paramfile ../protein/.
        #./genRunScripts_bckend.sh $paramfile
        #./genRunScripts_bckend_parstage2.sh $paramfile
        ./genRunScripts_bckend_parstage1_parstage2.sh $paramfile
        #echo "cp -R ../proteinsparams/$paramfile ../protein/."
        #echo "./genRunScripts_bckend.sh $paramfile"
        #./transferForEM.sh $paramfile
else
	paramwisetrial/makeFolderForParams.sh $protname
        #echo "paramwisetrial/makeFolderForParams.sh $protname"
        folderran="../protein/"$protname"_allparam/"$protname"_findwhichtorun.txt"
        #if [ -f $folderan ]; then
	#        argforem=`awk '{print $1}' $folderran | paste -s -d, -`
        #        ./transferForEM.sh $argforem
        #else
        #        echo "WARNING: $folderran not found"
	#fi
             
fi
