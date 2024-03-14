#!/bin/bash

prot_list="$1"
upl_file="$2"
#dest_folder="$2"
prot_name="$3"

#work_folder="md_related/$prot_name"
#mdp_folder="md_related"

work_folder="md_related_v2/$prot_name"
mdp_folder="md_related_v2"

if [ ! -f "$upl_file" ]; then
	echo "ERROR: $upl_file not present"
	exit 1
fi

if [ ! -d $work_folder ]; then
	mkdir $work_folder
	if [[ "$?" -ne 0 ]]; then
		echo "ERROR: folder creating error $work_folder"
		exit 1
	fi	
fi

#######################
	echo "Run the following shell scripts"
	IFS=',' read -r -a array <<< "$prot_list"

	c=0
	for prot in "${array[@]}"
	do
                #echo $prot
		c=$((c+1))
		if [ -f $prot ]; then
			mkdir $work_folder/$c
			cp $prot $work_folder/$c/.
			upl_file_name=`echo $upl_file | awk -F'/' '{print $NF}'`
			cp $upl_file $work_folder/$c/.
                        cp $mdp_folder/minim_disre_posre_long.mdp $work_folder/$c/.
                        cp $mdp_folder/ions.mdp $work_folder/$c/.
                        cp $work_folder/../distResItpFromUpl.sh $work_folder/$c/.
                        cp $work_folder/../nomanclature_edit.csv $work_folder/$c/.
                        cp $work_folder/../aminoacid.txt $work_folder/$c/.
			cp $work_folder/../getPDBFromBMRB_edit.sh $work_folder/$c/.
			cp $work_folder/../pseudo.csv $work_folder/$c/.
                        
			protein_name=`echo $prot | awk -F'/' '{print $NF}' | awk -F'.pdb' '{print $1}'`
                        filename=$protein_name"_"$c"_em.sh"
                        echo "#!/bin/bash" > $work_folder/$c/$filename
                        echo "protein_name=$protein_name" >> $work_folder/$c/$filename
                        echo "work_folder=./" >> $work_folder/$c/$filename
                        echo "mdp_folder=./" >> $work_folder/$c/$filename
                        echo "upl_file=$upl_file_name" >> $work_folder/$c/$filename
                        cat $mdp_folder/run_template.sh >> $work_folder/$c/$filename
			
			echo "    $work_folder/$c/$filename"
                        chmod 755 $work_folder/$c/$filename
                        if [[ $? -eq 0 ]]; then
                                curr_path_now=`pwd`
                                cd $work_folder/$c
                                #./$filename &
				./$filename
                                cd $curr_path_now
                        fi
			#wait

                        #cat run_template.sh >> $work_folder/$filename
                       
			#gmx pdb2gmx -f "$work_folder/$protein_name".pdb -water tip3p -o "$work_folder/$protein_name"_process.gro -ignh <<< 8
                        #./distResItpFromUpl.sh $upl_file "$work_folder/$protein_name"_process.gro > "$work_folder/$protein_name"_disre.tip
                        #cp $work_folder"/topol.top" $work_folder"/topol.top_bkp"                                                
                else
                        echo "Warning: File: $prot not found"
		fi
	done
        #wait
