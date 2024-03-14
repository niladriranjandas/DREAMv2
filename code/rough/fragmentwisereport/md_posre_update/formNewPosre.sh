#!/bin/bash

# ./formNewPosre.sh posre.itp /root/Downloads/gmx_22nov17/niladri_10_64_35_24/1xxe_anchor_pseudo_refine_2_dpos/1xxe_cluster_0001_try1_refine_2_md.gro 1xxe.upl  1xxe_newposre.itp
usage () {
	echo "./formNewPosre.sh <posre.itp file> <gro file> <upl file> <new posre.itp file>"
}

##########
if [[ "$#" -ne 4 ]]; then
	echo "too few arguments"
        usage
        exit -1
fi

posrefile="$1"
grofile="$2"
uplfile="$3"
newposrefile="$4"

if [ ! -f "$posrefile" ]; then
  	echo "$posrefile not present"
        exit -1
fi
if [ ! -f "$grofile" ]; then
  	echo "$grofile not present"
        exit -1
fi
if [ ! -f "$uplfile" ]; then
  	echo "$uplfile not present"
        exit -1
fi
if [ -f "$newposrefile" ]; then
	mv $newposrefile $newposrefile"_bkp"
fi	
###########
	if [ -f gro_index_file ]; then
		rm gro_index_file
	fi
	./getIndexFromGro.sh $uplfile $grofile gro_index_file
	if [ -f newposre_dup.itp ]; then
		rm newposre_dup.itp
	fi
	python includeInPosreitp.py $posrefile gro_index_file newposre_dup.itp
		
	lineno=`egrep -m1 -nr [[:space:]][0-9] $posrefile | cut -f1 -d:`
	cat $posrefile | awk -v lineno=$lineno '{if (NR < lineno) print $0}' > $newposrefile
	cat newposre_dup.itp | uniq >> $newposrefile

        rm newposre_dup.itp






