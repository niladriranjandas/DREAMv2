#!/bin/bash

# for the distance restraints file 
#./distResItpFromUpl.sh 1xxe.upl /root/Downloads/gmx_22nov17/niladri_10_64_35_24/1xxe_anchor_pseudo_refine_2_dpos/1xxe_cluster_0001_try1_refine_2_md.gro > 1xxe_distres.tip
usage () {
	echo "./distResItpFromUpl.sh <uplfile> <grofile> <[output file]>"
}

getIndexFromGro() {
	resii="$1"
        resiname="$2"
        atomname="$3"
        groFile="$4"

        resi=`egrep [[:blank:]]$resii$resiname $groFile | awk -v atm=$atomname '{if ($2==atm) print $3}'`
        if [ -z "$resi" ]; then
		echo "XX $resii $resiname $atomname" 
        else
                echo $resi
        fi
}

####
if [[ "$#" -lt 2 ]]; then
	echo "Error: Too few arguments"
	usage
	exit -1
fi

uplfile="$1"
grofile="$2"

if [[ "$#" == 3 ]]; then
	outfile="$3"
fi

if [ ! -f $uplfile ]; then
	echo "$uplfile not found"
	exit -1
fi

if [ ! -f $grofilefile ]; then
	echo "$grolfile not found"
	exit -1
fi

diststrong=2.7
distmedium=3.7
distweak=6.5
#distweak=5.5

STRONG_WT=100
MEDIUM_WT=10
WEAK_WT=5

count=0
###
	echo "[ distance_restraints ]"
	echo "; ai aj type index type’ low up1 up2 fac"
	while IFS='' read -r line || [[ -n "$line" ]]; do
		resi_i=`echo $line | awk '{print $1}'`
                resi_name_i=`echo $line | awk '{print $2}'`
		atom_i=`echo $line | awk '{print $3}'`
                
		resi_j=`echo $line | awk '{print $4}'`
                resi_name_j=`echo $line | awk '{print $5}'`
		atom_j=`echo $line | awk '{print $6}'`
                              
                d_ij=`echo $line | awk '{print $7}'`
		
		mapped_atm_i=`./getPDBFromBMRB_edit.sh $resi_name_i $atom_i`
		mapped_atm_j=`./getPDBFromBMRB_edit.sh $resi_name_j $atom_j`

		if [ $mapped_atm_i != 'XX' ] && [ $mapped_atm_j != 'XX' ]; then
			gro_atm_i=''
			gro_atm_j=''
			gro_atm_i=`getIndexFromGro "$resi_i" "$resi_name_i" "$mapped_atm_i" "$grofile"`
			gro_atm_j=`getIndexFromGro "$resi_j" "$resi_name_j" "$mapped_atm_j" "$grofile"`	
                        #echo "$line"			
			#printf "\ngro_atm_i=%s\tgro_atm_j=%s" $gro_atm_i $gro_atm_j
			#printf "\nresi_i=%s resi_name_i=%s mapped_atm_i=%s \n" $resi_i $resi_name_i $mapped_atm_i
			#printf "\nresi_j=%s resi_name_j=%s mapped_atm_j=%s \n" $resi_i $resi_name_j $mapped_atm_j

                        if [ ! -z "$gro_atm_i" ] && [ ! -z "$gro_atm_j" ]; then	
				if [ $gro_atm_i != 'XX' ] && [ $gro_atm_j != 'XX' ]; then
					flag=`awk -v dij=$d_ij -v diststrong=$diststrong -v distmedium=$distmedium -v distweak=$distweak 'BEGIN{if (dij>0 && dij<=diststrong) 
												print 1; 
											else if (dij>diststrong && dij<=distmedium)
												print 2;
									                else if (dij>distmedium && dij<=distweak)
												print 3;}'`
									                #else if (dij>distmedium)
											#	print 3;}'
										   	#else if (dij>distmedium && dij<=distweak)
											#	print 3;}'`
					#printf "%d\t%d 10 0.0 0.25 0.3 120.0\t%d\n" $gro_atm_i $gro_atm_j $flag
					if [[ $flag == 1 ]]; then
						count=$((count+1))
						printf "%d\t%d\t1\t%d\t1\t0\t0.2\t0.27\t%d\n" $gro_atm_i $gro_atm_j $count $STRONG_WT
						#printf "%d\t%d 10 0.0 0.2 0.27 120.0\n" $gro_atm_i $gro_atm_j
						#printf "%d\t%d 10 0.0 0.25 0.3 120.0\n" $gro_atm_i $gro_atm_j
					elif [[ $flag == 2 ]]; then
						count=$((count+1))
        	                                printf "%d\t%d\t1\t%d\t1\t0\t0.28\t0.35\t%d\n" $gro_atm_i $gro_atm_j $count $MEDIUM_WT
						#printf "%d\t%d 10 0.0 0.28 0.35 120.0\n" $gro_atm_i $gro_atm_j
						#printf "%d\t%d 10 0.0 0.3 0.35 120.0\n" $gro_atm_i $gro_atm_j
					elif [[ $flag == 3 ]]; then
						count=$((count+1))
						printf "%d\t%d\t1\t%d\t1\t0\t0.36\t0.55\t%d\n" $gro_atm_i $gro_atm_j $count $WEAK_WT
						#printf "%d\t%d 10 0.0 0.36 0.55 120.0\n" $gro_atm_i $gro_atm_j
						#printf "%d\t%d 10 0.0 0.4 0.5 120.0\n" $gro_atm_i $gro_atm_j
					fi
				fi
			fi
		fi
	done < "$uplfile"
