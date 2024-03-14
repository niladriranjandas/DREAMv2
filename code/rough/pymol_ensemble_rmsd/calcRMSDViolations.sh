basedir="/home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/????/"
result_folder="proteins_no_sa"

createFolder() {

	basedir="$1"
	result_folder="$2"

	if [ ! -d $basedir ]; then
		error "folder doesnot exist"
		return 1
	fi

	[[ ! -d $result_folder ]] && mkdir $result_folder || printf "\nFolder %s already exits" $result_folder

	for i in `ls -d $basedir/*/`
	do
		echo "--------------------------------------------------"
		protein_name=`echo $i | awk -F'/' '{print $(NF-1)}'`
		printf "\n\t%s" $protein_name

		mkdir $result_folder/$protein_name
		cp $i"gromacs_stuff/"*"pdb" $result_folder/$protein_name/.		
	done	
}

runViolationCheck() {
	basedir="$1"
	upl_lo_dir="$2"
	result_folder="$3"

	upl_ptrn="*.upl"
	lo_ptrn="*.lo"

	for i in `ls -d "$basedir/"*"/"`
	do
		protein_name=`echo $i | awk -F'/' '{print $(NF-1)}'`

		pdbs_=`ls "$i"*"pdb"`
		pdbs=""
                flag=1
		for pdbs_ in `ls "$i"*"pdb"`
		do
                     count=$((count+1))
		     pdb_name=`echo $pdbs_ | awk -F'/' '{print $NF}' | awk -F'.pdb' '{print $1}'`
		     [[ $flag -eq 1 ]] && pdbs=$pdb_name || pdbs=$pdbs","$pdb_name
		     flag=0
		done
		#pdbs=`echo $pdbs_ | sed 's/ /,/g'`

		upl_file=`find "$upl_lo_dir/$protein_name" -name "$upl_ptrn"`	
		lo_file=`find "$upl_lo_dir/$protein_name" -name "$lo_ptrn"`

		[[ -f tmp_upl.upl ]] && rm tmp_upl.upl
		[[ -f tmp_lo.lo ]] && rm tmp_lo.lo
		for tmp in `echo $upl_file`
		do
		     cat $tmp >> tmp_upl.upl
		done
		for tmp in `echo $lo_file`
		do
		    cat $tmp >> tmp_lo.lo
		done
		#if [ ! -d $result_folder/$protein_name ]; then
		#	mkdir $result_folder/$protein_name
		#fi
		
		printf "\n%s\t" $protein_name
		#pymol getViolations.py -- $i $pdbs $upl_file $lo_file | grep rmsd_up_avg
		pymol getViolations.py -- $i $pdbs tmp_upl.upl tmp_lo.lo | grep rmsd_up_avg
		rm tmp_upl.upl tmp_lo.lo
	done
}

prepEnsembleCalc() {
	pdb_name="$1"
	pdbs="$2"					 #pdb_names separated by commas
	result_folder="$3"
	IFS=',' read -r -a pdbs_arr <<< "$pdbs"
	
	op_name=$pdb_name"_ensemble.pdb"
	c=0
	for i in "${pdbs_arr[@]}"
	do
	  c=$((c+1))
	  printf "MODEL\t%s\n" $c >> $op_name
	  end_l=`grep -n END $i | awk -F':' '{print $1}'`
	  awk -v num=$end_l 'FNR>=1 && FNR<num' $i >> $op_name
	  echo "ENDMDL" >> $op_name
	done
	echo "END" >> $op_name
	mv $op_name $result_folder/.
}

ensembleResiwiseCompare() {
	basedir="$1"
	protein_name="$2"	# pdb name without the .pdb extension
	result_folder="$3"
	
	run_resi_compare=1
	run_ensemble_avg_rmsd=1

	if [ $run_resi_compare -eq 1 ]; then
		printf "\n-----calculating the residue-wise rmsd--------\n"		
		echo "==========residue wise rmsd ============" >> log_file.txt	
		pymol -cq ensembleCompare.py -- $basedir $protein_name >> log_file.txt
		pair_rmsd_file=`ls $protein_name*pair_rmsd.dat`
		avg_rmsd_file=`ls $protein_name*avg_rmsd.dat`

		if [ ! -f $pair_rmsd_file ]; then
			printf "\n\tFile %s not created" $pair_rmsd_file
		else	
			mv $pair_rmsd_file $result_folder/.
		fi
		if [ ! -f $avg_rmsd_file ]; then
			printf "\n\tFile %s not created" $avg_rmsd_file
		else
			mv $avg_rmsd_file $result_folder/.
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
			mv $rmsd_file $result_folder/.	
		fi	
	fi
}

##----------------------------------- no_sa ------------------------------------------------------------------------------------------------
#createFolder "/home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive" "exp_no_sa"
#runViolationCheck "proteins_no_sa" "/home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive" "protein_nosa"
#runViolationCheck "exp_no_sa" "/home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive" "exp_no_sa" >> dist_violation_nosa_v2.txt

##--------prepare ensemble ---------------------------------------
#mkdir exp_no_sa_ensemble
#for i in `ls -d exp_no_sa/*/`
#do
#	pdbs_=`ls "$i"*".pdb"`
#	pdbs=`echo $pdbs_ | sed 's/ /,/g'`
#
#        protein_name=`echo $i | awk -F'/' '{print $(NF-1)}'`
#	echo $protein_name
#	mkdir exp_no_sa_ensemble/$protein_name
#	prepEnsembleCalc $protein_name $pdbs exp_no_sa_ensemble/$protein_name	
#done		

##------ensemble compare --------------------------------------------
#for i in `ls -d exp_no_sa_ensemble/*/`
#do
#	pdb_ensemble=`ls "$i"*".pdb"`
#	pdb_name=`echo $pdb_ensemble | awk -F'/' '{print $(NF)}' | awk -F'.pdb' '{print $1}'`
	
#	echo $pdb_name
#	ensembleResiwiseCompare $i $pdb_name $i
#done	

#ensembleResiwiseCompare exp_no_sa_ensemble/1a7m/ 1a7m_ensemble exp_no_sa_ensemble/1a7m/


##----------------------------------- sa ------------------------------------------------------------------------------------------------
#mkdir exp_sa
#for pdbs in `ls -d /root/Downloads/compare_results/transfer_PDBS/*/`
#do 
#	pdb_name=`echo $pdbs | awk -F'/' '{print $(NF-1)}'`
#	echo $pdb_name 
#	mkdir exp_sa/$pdb_name
#	cp "$pdbs"*".pdb" exp_sa/$pdb_name/.
#done

prepEnsembleCalc 2m4k "/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/gromacs_pos_restraint/2m4k/script_1/2m4k_1_sarestraint.pdb,/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/gromacs_pos_restraint/2m4k/script_2/2m4k_2_sarestraint.pdb,/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/gromacs_pos_restraint/2m4k/script_3/2m4k_3_sarestraint.pdb,/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/gromacs_pos_restraint/2m4k/script_4/2m4k_4_sarestraint.pdb,/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/gromacs_pos_restraint/2m4k/script_5/2m4k_5_sarestraint.pdb" "/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/gromacs_pos_restraint/2m4k"

#runViolationCheck "exp_sa" "/home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive" "exp_sa" >> dist_violation_sa_v2.txt

##------prepare ensemble ---------------------
#mkdir exp_sa_ensemble

#for pdbs in `ls -d /root/Downloads/compare_results/transfer_PDBS/*/`
#do
#	pdb_name=`echo $pdbs | awk -F'/' '{print $(NF-1)}'`
#	echo $pdb_name
#	mkdir exp_sa_ensemble/$pdb_name
#	c=0
#	for pdb in `ls "$pdbs"*"pdb"`
#	do
#	    c=$((c+1))
#	    printf "MODEL%10d\n" $c >> exp_sa_ensemble/$pdb_name/$pdb_name"_ensemble.pdb"
#	    strt_l=`grep -n MODEL $pdb | awk -F':' '{print $1}'` 
#	    end_l=`grep -n ENDMDL $pdb | awk -F':' '{print $1}'` 
#	    awk -v a=$strt_l -v b=$end_l 'FNR>a && FNR<=b' $pdb >> exp_sa_ensemble/$pdb_name/$pdb_name"_ensemble.pdb"
#	done
#done

##------ensemble compare --------------------------------------------
#for i in `ls -d exp_sa_ensemble/*/`
#do
#	pdb_ensemble=`ls "$i"*".pdb"`
#	pdb_name=`echo $pdb_ensemble | awk -F'/' '{print $(NF)}' | awk -F'.pdb' '{print $1}'`

#	echo $pdb_name
#	ensembleResiwiseCompare $i $pdb_name $i
#done	


