#!/bin/bash

usage (){
	echo "usage: ./getindexFromGro.sh <uplfile> <grofile> outputfile"
}

getIndexFromGro() {
	resii="$1"
        resiname="$2"
        atomname="$3"
        groFile="$4"

        resi=`egrep [[:blank:]]$resii$resiname $groFile | awk -v atm=$atomname '{if ($2==atm) print $3}'`
        if [ -z "$resi" ]; then
		echo "XX $resii $resiname $atomname" >> $err_file
        else
                echo $resi
        fi
}

uplfile="$1"
grofile="$2"
opfile="$3"
err_file=$opfile"_errfile"

dist_tol='3.5'

        if [ -f "$opfile" ]; then
		mv $opfile $opfile"_bkp"
        fi

	while IFS='' read -r line || [[ -n "$line" ]]; do
		resi_i=`echo $line | awk '{print $1}'`
                resi_name_i=`echo $line | awk '{print $2}'`
		atom_i=`echo $line | awk '{print $3}'`
                
		resi_j=`echo $line | awk '{print $4}'`
                resi_name_j=`echo $line | awk '{print $5}'`
		atom_j=`echo $line | awk '{print $6}'`
                
              
                d_ij=`echo $line | awk '{print $7}'`
                
                flag=`awk -v dij=$d_ij -v disttol=$dist_tol 'BEGIN{if (dij<=disttol) print 1; else print 0}'`
                if [[ $flag == 1 ]]; then
			#echo "$resi_i $resi_name_i $atom_i $resi_j $resi_name_j $atom_j $d_ij" >> testfile
	                mapped_atm_i=`./getPDBFromBMRB_edit.sh $resi_name_i $atom_i`
			mapped_atm_j=`./getPDBFromBMRB_edit.sh $resi_name_j $atom_j`
	
        	        if [ $mapped_atm_i == 'XX' ]; then
        	                 echo $resi_name_i $atom_i XX >> $err_file
        	        else
  			         getIndexFromGro $resi_i $resi_name_i $mapped_atm_i $grofile >> $opfile
        	        fi
	                if [ $mapped_atm_j == 'XX'  ]; then
	                         echo $resi_name_j $atom_j XX >> $err_file
	                else
			         getIndexFromGro $resi_j $resi_name_j $mapped_atm_j $grofile >> $opfile
	                fi
		fi

	done < "$uplfile"


