#!/bin/bash

prot_list="$1"
upl_file="$2"
#dest_folder="$2"
prot_name="$3"

#work_folder="md_related/$prot_name"
#mdp_folder="md_related"

work_folder="md_related_v2/$prot_name"
mdp_folder="md_related_v2"

#######################
	echo "Run the following shell scripts"
	IFS=',' read -r -a array <<< "$prot_list"

	c=0
	for prot in "${array[@]}"
	do
                #echo $prot
		c=$((c+1))
		if [ -f $prot ]; then
			protein_name=`echo $prot | awk -F'/' '{print $NF}' | awk -F'.pdb' '{print $1}'`
			### file for md run ###
                        filename2=$protein_name"_"$c"_em_md.sh"
			#######################

                        if [[ $? -eq 0 ]]; then
                                curr_path_now=`pwd`
                                cd $work_folder/$c
                                ./$filename2
                                cd $curr_path_now
                        fi


                        #cat run_template.sh >> $work_folder/$filename
                       
			#gmx pdb2gmx -f "$work_folder/$protein_name".pdb -water tip3p -o "$work_folder/$protein_name"_process.gro -ignh <<< 8
                        #./distResItpFromUpl.sh $upl_file "$work_folder/$protein_name"_process.gro > "$work_folder/$protein_name"_disre.tip
                        #cp $work_folder"/topol.top" $work_folder"/topol.top_bkp"                                                
                else
                        echo "Warning: File: $prot not found"
		fi
	done
