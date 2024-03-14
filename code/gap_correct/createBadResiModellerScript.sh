#!/bin/bash

generateSeqForFill () {
        seq="$1"
        modelresi=","$modelresi_","

        amino3to1="gap_correct/aminoacid.txt"

        c=0
        while IFS= read -r line || [[ -n "$line" ]]; do
            c=$((c+1))
            resi=`echo $line | awk '{print $1}'`

            resi3to1_=`grep -i $resi $amino3to1`
            if [ ! -z "$resi3to1_" ];  then
		resi3to1=`echo $resi3to1_ | awk '{print $2}'`
                if [[ $c -gt 75 ]]; then
                        c=1
                        printf "\n%s" "$resi3to1"
                else
                                printf "%s" "$resi3to1"
                fi
            else
                echo "ERROR"
            fi
        done < "$seq"
	if [[ $c -gt 75 ]]; then
		printf "\n*"
	else
		printf "*"
	fi
        echo ""
}

generateSeqForGaps () {
	seq="$1"
	modelresi_="$2"
	modelresi=","$modelresi_","
        amino3to1="gap_correct/aminoacid.txt"

	c=0
	while IFS= read -r line || [[ -n "$line" ]]; do
	    c=$((c+1))
	    resi=`echo $line | awk '{print $1}'`
	    resi_num=`echo $line | awk '{print $2}'`
	    srchstr=","$resi_num","
    
	    resi3to1_=`grep -i $resi $amino3to1`
	    if [ ! -z "$resi3to1_" ];  then
	        resi3to1=`echo $resi3to1_ | awk '{print $2}'`
	        if [[ $c -gt 75 ]]; then
	        	c=1
	        	if [[ $modelresi == *"$srchstr"* ]]; then
	        		printf "\n-"
	        	else
        		printf "\n%s" "$resi3to1"
	        	fi
	        else
	            	if [[ $modelresi == *"$srchstr"* ]]; then	
	            		printf "-" 
	            	else
	            		printf "%s" "$resi3to1"	
	            	fi
	        fi
	    else
	    	echo "ERROR"
	    fi    
	done < "$seq"
	if [[ $c -gt 75 ]]; then
		printf "\n*"
	else
		printf "*"
	fi
	echo ""
}

deletePDBGaps () {
	pdbfile="$1"
	modelresi_="$2"
	modelresi=","$modelresi_","
	c=0
	while IFS= read -r line || [[ -n "$line" ]]; do
		beginstr=`echo $line | awk '{print $1}'`
		if [[ $beginstr == "ATOM" ]]; then
			resi=`echo $line | awk '{print $5}'`
			srchstr=","$resi","
			if [[ $modelresi == *"$srchstr"* ]]; then
				c=$((c+1))			
			else
				echo "$line"
			fi
		else
			echo "$line"
		fi
	done < "$pdbfile"
}

genAliFile () {
	proteinname="$1"
	seqfile="$2"
	pdbfile="$3"
	gaps="$4"
	beginresi="$5"

	newpdbfile_=`echo $pdbfile | awk -F'pdb' '{print $1}'`
	newpdbfile=$newpdbfile_"gap.pdb"
	newpdbfile_noext_=$newpdbfile_"gap"
	newpdbfile_noext=`echo $newpdbfile_noext_ | awk -F'/' '{print $NF}'`

	firstresi=`head -1 $seqfile | awk '{print $2}'`
	lastresi=`tail -1 $seqfile | awk '{print $2}'`
	lengthseq=$((lastresi-firstresi+1))
	# -- prepare header for ali file -- #
	printf ">P1;%s\n" "$newpdbfile_noext"
	printf "structureX:%s:%4d : :+%-5d: :::-1.00:-1.00\n" "$newpdbfile_noext"  "$beginresi" "$lengthseq"
        # -- generate seq for gaps -- #
	generateSeqForGaps "$seqfile" "$gaps"
	printf "\n>P1;%s_fill\nsequence:::::::::\n" "$proteinname"
        generateSeqForFill "$seqfile"

        # -- PDB file with gaps -- #	
        deletePDBGaps "$pdbfile" "$gaps" > "$newpdbfile"
}

breakGapsIntoSeqs () {
	gaps="$1"

	IFS=', ' read -r -a array <<< "$gaps"


	c=0
	len=${#array[@]}
	len_=$((len-1))
	
	for index in "${!array[@]}"
	do
		if [[ $index -eq $len_ ]]; then
			printf ",%s" ${array[index]}
		else
			if [[ $index -eq 0 ]]; then
				printf "%s" ${array[index]}
			fi
			curr=${array[index]}
			index_next=$((index+1))
			next=${array[index_next]}
			if [[ $((curr+1)) -ne $next ]]; then
				printf ",%s\n%s" "$curr" "$next"
			fi
		fi
		c=$((c+1))
	done
}

makeModelFile () {
	gaps="$1"

	protname="$2"
	org_num="$3"
	alignfile="$4"
	gappdb="$5"
	seqpdb="$6"
	tmpfile=$protname"_gaps"
	breakGapsIntoSeqs $gaps > $tmpfile

	gappdb_noext=${gappdb::-4}
        # --- write header ----#
        printf "\nfrom modeller import *\nfrom modeller.automodel import *    # Load the automodel class" 	
        printf "\nlog.verbose()\nenv = environ()\n"
        printf "\n# directories for input atom files\nenv.io.atom_files_directory = ['.', '../atom_files']\n"
        if [[ $org_num -eq 1 ]]; then
        		#printf "class MyModel(allhmodel):\n\tdef select_atoms(self):\n\t\treturn selection("
                        printf "class MyModel(allhmodel):\n\tdef special_patches(self, aln):\n"
                        printf "\t\tself.rename_segments(segment_ids='', renumber_residues=%d)\n" "$org_num"
                        printf "\tdef select_atoms(self):\n\t\treturn selection("
        else
                 	printf "class MyModel(allhmodel):\n\tdef special_patches(self, aln):\n"
        		printf "\t\tself.rename_segments(segment_ids='', renumber_residues=%d)\n" "$org_num"
       	        	printf "\tdef select_atoms(self):\n\t\treturn selection("
	fi

	last_count=`wc -l $tmpfile | awk '{print $1}'`
	count=0
	while IFS= read -r line || [[ -n "$line" ]]; do
		resi_i=`echo $line | awk -F',' '{print $1}'`
		resi_j=`echo $line | awk -F',' '{print $2}'`
		count=$((count+1))
		if [[ $count -eq 1 ]]; then
			printf "self.residue_range('%d','%d'),\n" "$resi_i" "$resi_j"
		elif [[ $count -gt $last_count ]]; then 
			printf "\t\t\t\t self.residue_range('%d','%d'))\n" "$resi_i" "$resi_j"
		else
			printf "\t\t\t\t self.residue_range('%d','%d'),\n" "$resi_i" "$resi_j"
		fi
	done < "$tmpfile"
	printf "\n\na = MyModel(env, alnfile = '%s',\nknowns = '%s', sequence = '%s')" "$alignfile" "$gappdb_noext" "$seqpdb"
        printf "\na.starting_model= %d\na.ending_model  = %d" "1" "5"
        printf "\n\na.make()\n" 
}

writeAliFile () {
	proteinname="$1"
	seqfile="$2"
	pdbfile="$3"
	gap="$4"
	beginresi="$5"

	if [ ! -f "$seqfile" ]; then
		echo "ERROR: $seqfile NOT found"
		exit -1
	fi
	if [ ! -f "$pdbfile"  ]; then
		echo "ERROR: $pdbfile NOT found"
		exit -1
	fi
      
       # --- define o/p file names --- #
       #alifile=$proteinname"_20.ali"
       #modelfile=$proteinname"_model20.py"
       alifile=$proteinname"_badresi.ali"
       modelfile=$proteinname"_badresi.py"

       # --- generate the ali file --- #
       genAliFile $proteinname $seqfile $pdbfile $gap $beginresi > $alifile

       # --- generate model.py file ---#
       seqbegin=`head -1 $seqfile | awk '{print $2}'`
       pdbname=`echo $pdbfile | awk -F'/' '{print $NF}' | awk -F'pdb' '{print $1}'`
       gappdb=$pdbname"gap.pdb"
       fillfile=$proteinname"_fill"
       
       makeModelFile $gap $protname $seqbegin $alifile $gappdb $fillfile > $modelfile
}

# ----- main call for 20 models ----- #
protname="$1"
seqfile="$2"
pdbfile="$3"
gaps="$4"
beginresi="$5"

writeAliFile "$protname" "$seqfile" "$pdbfile" "$gaps" "$beginresi"
