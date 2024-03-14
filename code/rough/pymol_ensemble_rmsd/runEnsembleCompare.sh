#!/bin/bash

basedir="proteins"
run_resi_compare=1
run_ensemble_avg_rmsd=1

for i in `ls $basedir`
do
	protein_name_=`echo $i | awk -F'/' '{print $NF}'`
        protein_name=`echo $protein_name_ | awk -F'.pdb' '{print $1}'`

	echo $protein_name

        if [ ! -d $basedir/$protein_name ]; then
  		mkdir $basedir/$protein_name
        fi
	if [ ! -f $basedir/$protein_name/$i ]; then
		cp $basedir/$i $basedir/$protein_name/.
        fi

        if [ $run_resi_compare -eq 1 ]; then
		printf "\n-----calculating the residue-wise rmsd--------\n"		
		echo "==========residue wise rmsd ============" >> log_file.txt	
		pymol -cq ensembleCompare.py -- $basedir $protein_name >> log_file.txt
		pair_rmsd_file=`ls $protein_name*pair_rmsd.dat`
		avg_rmsd_file=`ls $protein_name*avg_rmsd.dat`

		if [ ! -f $pair_rmsd_file ]; then
			printf "\n\tFile %s not created" $pair_rmsd_file
		else	
			mv $pair_rmsd_file $basedir/$protein_name/.
		fi
		if [ ! -f $avg_rmsd_file ]; then
			printf "\n\tFile %s not created" $avg_rmsd_file
		else
			mv $avg_rmsd_file $basedir/$protein_name/.
		fi
        fi

	if [ $run_ensemble_avg_rmsd -eq 1 ]; then
		printf "\n-----calculating avg rmsd for ensemble--------\n"		
                echo "========avg rmsd========" >> log_file.txt
		pymol -cq ensembleAvgRMSD.py -- $basedir $protein_name >> log_file.txt
		rmsd_file=`ls $protein_name*rmsd_file.txt`

		if [ ! -f $rmsd_file ]; then
			printf "\n\tFile %s does not exist" $rmsd_file
		else
			#add the 1ts column in the rmsd file to get the avg rmsd
			num_lines=`wc -l $rmsd_file | awk -F' ' '{print $1}'`
			sed 's/^(\|)$//g' $rmsd_file | awk -F',' -v len=$num_lines '{sum+=$i}END{avg = sum/len;print avg}' >> $rmsd_file
			mv $rmsd_file $basedir/$protein_name/.	
		fi	
	fi
done
