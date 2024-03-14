#!/bin/bash

#./runConfPlotinRoot.sh /root/Downloads/compare_results/transfer_PDBS
folder="$1"

	mkdir conf_plot
	for full_path in `ls -d $folder/????/`
	do
		protein=`echo $full_path | awk -F'/' '{print $((NF-1))}'`
		mkdir -p conf_plot/$protein
		arg_pdbs_=`ls $full_path*pdb`
		arg_pdbs=`echo $arg_pdbs_ | tr '\n' ' '`	
		./doConfPlotForPDBs.sh $protein $arg_pdbs
		if [ -f "$protein.eps" ]; then
			printf "\n File generate: %s.eps Coppying to protein folder\n" $protein
			mv "$protein.eps" conf_plot/$protein/.
		else
			printf "\nFile %s.eps not found\n" $protein
		fi
	done
