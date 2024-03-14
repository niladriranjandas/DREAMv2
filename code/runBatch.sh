#!/bin/bash

#-------------------------------------#
# usage: ./runBatch.sh -file<txt file containing protein list sep by newline>
#        ./runBatch.sh -list protein1,protein2,...
#-------------------------------------#

usage () {
	echo "usage ./runBatch.sh -file <txt file having protein list one in each line>"
	printf "\tor\n"
	echo "./runBatch.sh -list protein1,protein2,..."
  exit 1
}

deleteAll() {
  #dir_list="$1"
protein="$1"
declare -a  dir_list=('../output_pdb' '../output_plots' '../log' 'modeller_stuff')

	if [ ${#dir_list[@]} -eq 0 ]; then
		printf "\nDirectory list for deleting not found"
		exit 1
	fi
	
	for dirs in ${dir_list[@]}
	do
		if [ ! -d $dirs ]; then
			printf "\n Dir %s not present " $dirs
			exit 1
		fi
		printf "\n\t deleting dir %s \n" $dirs
		files=`ls $dirs/*"$protein"*`
		if [ ! -z "$files" ]; then
			rm $dirs/*"$protein"*
			if [ $? -ne 0 ]; then
				printf "\n Error deleting dir %s " $dirs
				exit 1
			fi
		fi
	done
}

if [ $# -ne 2 ]; then
	usage
fi
	
	if [ "$1" = "-file" ]; then
		if [ ! -f $2 ]; then
			printf "\n File %s doesn't exit\n" $2
			exit 1
		fi
		c=0
		declare -a protein_list
		while IFS='' read -r line || [[ -n "$line" ]]; do
			c=$((c+1))
			protein_list[$c]=$line
		done < "$2"
	elif [ "$1" = "-list" ]; then
		IFS=',' read -r -a protein_list <<< $2	
	fi

	if [ ${#protein_list[@]} -eq 0 ]; then
		printf "\n ERROR: No protein list parsed "
		exit 1
	fi

BASE_DIR=`pwd`
declare -a MANDITORY_FILES_GMX=('minim.mdp' 'nvt.mdp' 'npt.mdp' 'ions.mdp' 'md.mdp' 'steps_taken_gromacs_template.txt')

#----------------- delete contents of arrays ------------------------#
#declare -a DEL_FILES_FULL=('../output_pdb' '../output_plots' '../log' 'modeller_stuff')

#----------loop through protein list ------------------
	for protein in ${protein_list[@]}
	do	
           printf "\n ---------------------Protein: %s ---------------------------" $protein

		printf "\n\t Cleaning up residual files from previous runs..."
		printf "\n\t Cleaning the scripts stage1_$protein.m"
		rm stage1_"$protein".m
		printf "\n\t Cleaning the scripts stage2_$protein.m"
		rm stage2_"$protein".m
		#deleteAll "${DEL_FILES_FULL[@]}"
		deleteAll $protein

	    #---------not gromacs all but the manditory files---------------
		#mkdir gromacs_stuff/tmp
		#for files in ${MANDITORY_FILES_GMX[@]}
		#do
		#	cp gromacs_stuff/$files gromacs_stuff/tmp/.
		#	if [ $? -ne 0 ]; then
		#		printf "\n Couldn't move files for protein %s to tmp folder in gromacs_stuff" $protein
		#		exit 1
		#	fi		
		#done
		#
		#find gromacs_stuff/ -maxdepth 1 -type f -delete
		#if [ $? -ne 0 ]; then
		#	printf "\n Couldn't delete gromacs folder in protein %s " $protein
		#	exit 1
		#fi
		#cp gromacs_stuff/tmp/* gromacs_stuff/.
		#if [ $? -ne 0 ]; then
		#	printf "\n Couldn't copy gromacs tmp folder in protein %s " $protein
		#	exit 1
		#fi
		
	    #--------run script ----------------------------
		printf "\n\t Call our algo begin"
		##./genRunScripts_bckend.sh "$protein" > ../log/"$protein"_run.txt
		./genRunScripts_bckend_stage1prepstage2.sh "$protein" > ../log/"$protein"_run.txt
		if [ $? -ne 0 ]; then
			printf "\n Error running genRunScripts_bckend.sh for protein %s " $protein
		fi
		printf "\n\t our algo end"
	    #----------logs --------------------------------
		./listOutputs.sh "$protein" > ../log/"$protein"_files_created.txt

  	    #-----archive ----------------------------------
		printf "\n\t archiving..."
		./archiveProteinRun.sh "$protein" ../archive
		if [ $? -ne 0 ]; then
			printf "\n Error archiving for protein %s " $protein
			exit 1
		fi
		
		cp modeller_stuff/*"$protein"* ../archive/"$protein"/modeller_stuff/.
		cp gromacs_stuff/*"$protein"* ../archive/"$protein"/gromacs_stuff/.		

	done
	





