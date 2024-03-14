#!/bin/bash

generateSeqForFill () {
        local seq="$1"
        local modelresi=","$modelresi_","

        local amino3to1="gap_correct/aminoacid.txt"

        local c=0
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
        local seq="$1"
        local modelresi_="$2"
        local modelresi=","$modelresi_","

        local amino3to1="gap_correct/aminoacid.txt"

        local c=0
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
                                #printf "\n-"
				printf "\n%s" "$resi3to1"
                        else
                                #printf "\n%s" "$resi3to1"
				printf "\n-"
                        fi
                else
                        if [[ $modelresi == *"$srchstr"* ]]; then
                                #printf "-"
				printf "%s" "$resi3to1"
                        else
                                #printf "%s" "$resi3to1i"
				printf "-"
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
        local pdbfile="$1"
        local modelresi_="$2"
        local modelresi=","$modelresi_","
        local c=0
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
	local proteinname="$1"
	local seqfile="$2"
	local pdbfile="$3"
	local gaps="$4"
	local beginresi="$5"
	local newpdbfile_=`echo $pdbfile | awk -F'pdb' '{print $1}'`
	local newpdbfile=$newpdbfile"gap.pdb"

	local firstresi=`head -1 $seqfile | awk '{print $2}'`
	local lastresi=`tail -1 $seqfile | awk '{print $2}'`
	local lengthseq=$((lastresi-firstresi+1))
	# -- prepare header for ali file -- #
	printf ">P1;%s\n" "$pdbfile"
	printf "structureX:%s:%4d : :+%-5d: :::-1.00:-1.00\n" "$pdbfile"  "$beginresi" "$lengthseq"
	# -- generate seq for gaps -- #
	generateSeqForGaps "$seqfile" "$gaps"
}
includePDBresis () {
        local pdbfile="$1"
        local modelresi_="$2"
        local modelresi=","$modelresi_","
        local c=0
        while IFS= read -r line || [[ -n "$line" ]]; do
                beginstr=`echo $line | awk '{print $1}'`
                if [[ $beginstr == "ATOM" ]]; then
                        resi=`echo $line | awk '{print $5}'`
                        srchstr=","$resi","
                        if [[ $modelresi == *"$srchstr"* ]]; then
                                echo "$line"
                        else
                                c=$((c+1))
                        fi
                else
                        echo "$line"
                fi
        done < "$pdbfile"
}

writeAliFile () {
	local protname="$1"
	local seqfile="$2"
	local gapfile="$3"

	local c=0
	while IFS= read -r line || [[ -n "$line" ]]; do	
	   if [ ! -z "$line" ]; then
	        count=$((count+1))	
		pdbfilename=`echo $line | awk '{print $1}'`
		gaps=`echo $line | awk '{print $2}'`
                beginresi=`echo $line | awk '{print $3}'`
		gappdb_=`echo $pdbfilename | awk -F'/' '{print $NF}' | awk -F'pdb' '{print $1}'`
		gappdb=$gappdb_"gap_"$count

		#echo "$line"
		#echo "genAliFile $protname $seqfile $gappdb $gaps $beginresi"
                genAliFile $protname $seqfile $gappdb $gaps $beginresi
    
		if [[ "$c" -eq 1 ]]; then
			deletePDBGaps "$pdbfilename" "$gaps" > $gappdb".pdb"
                else
			includePDBresis "$pdbfilename" "$gaps" > $gappdb".pdb"
		fi
	   fi
	done < "$gapfile"
        printf "\n>P1;%s_multifill\nsequence:::::::::\n" "$protname"
        generateSeqForFill "$seqfile"
}

breakGapsIntoSeqs () {
        local gaps="$1"

        IFS=', ' read -r -a array <<< "$gaps"


        local c=0
        local len=${#array[@]}
        local len_=$((len-1))

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
        local gaps="$1"
        local protname="$2"
        local org_num="$3"
        local alignfile="$4"
        local gappdbs="$5"
        local seqpdb="$6"
        local tmpfile=$protname"_gaps_multi"
        breakGapsIntoSeqs $gaps > $tmpfile

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

        local last_count=`wc -l $tmpfile | awk '{print $1}'`
        local count=0
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
        printf "\n\na = MyModel(env, alnfile = '%s',\nknowns = (%s), sequence = '%s')" "$alignfile" "$gappdbs" "$seqpdb"
        printf "\na.starting_model= %d\na.ending_model  = %d" "1" "10"
        printf "\n\na.make()\n"
}


#### main function ####

protname="$1"
seqfile="$2"
gapfile="$3"
alifile="$4"
modfile="$5"
gapsorg="$6"

writeAliFile $protname $seqfile $gapfile > $alifile 

seqbegin=`head -1 $seqfile | awk '{print $2}'`
seqpdb=$protname"_multifill"

knowpdbs=""
pdbcount=0
while IFS= read -r line || [[ -n "$line" ]]; do
	pdbcount=$((pdbcount+1))
	curr_pdb=`echo $line | awk '{print $1}'`
	curr_pdb_=`echo $curr_pdb | awk -F'/' '{print $NF}' | awk -F'pdb' '{print $1}'`
	curr_pdb_gap=$curr_pdb_"gap_"$pdbcount
	knowpdbs=$knowpdbs",'"$curr_pdb_gap"'"
done < "$gapfile"
knowpdbs_="${knowpdbs:1}"
makeModelFile $gapsorg $protname $seqbegin $alifile $knowpdbs_ $seqpdb > $modfile
