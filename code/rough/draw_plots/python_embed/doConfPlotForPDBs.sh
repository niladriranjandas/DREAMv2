#!/bin/bash

#-------------------------------------------------
# usage: ./doConfPlotForPDBs pdb1.pdb pb2.pdb etc
#-------------------------------------------------

echo $#
if [ "$#" -lt 1 ]; then
	echo "usage./doConfPlotForPDBs.sh pdb1.pdb pb2.pdb etc"
	exit 1
fi


	args=""
	c=0	
	for i in "$@"
	do
        if [ $c == 0 ]; then
            args=$i
        else
            filename="test"$c".tsv"
            printf "\n\tGenerating temporary file: %s\n" $filename
            ./pdbtorsions_mod $i > $filename
            [[ $? -eq 0 ]] && args=$args" "$filename
        fi
            c=$((c+1))
    done

	./makeConfPlotOfEnsemble $args
	if [ $? != 0 ]; then
		printf "\ngenConfPlot.cpp failed for %s\n" "$args"
	else
		printf "\nPlots displayed\n"
		printf "\nDeleting the temporary files\n"
		IFS=" " read -r -a tmp_files <<< "$args"
		for files in "${tmp_files[@]}"
		do
			sudo rm -f $files
		done
	fi
