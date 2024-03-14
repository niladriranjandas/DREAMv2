#!/bin/bash
#usage: 

#for i in `grep  '^[[:space:]]*82[A-Z]' 2m4k_1_npt.gro | awk '{print $3}'`; do grep "^[[:blank:]]*$i" posre.itp_bkp; done

usage(){
	echo "usage: ./delLoopsPosre.sh gaps.txt 2m4k_1_npt.gro"
}

if [ "$#" -ne 2 ];then
	usage
	exit 1
fi

#file_name="gaps.txt"
#gro_file="2m4k_1_npt.gro"
file_name="$1"
gro_file="$2"

if [ ! -f $file_name ]; then
	echo "$file_name doesn't exist"
	exit 1
fi

if [ ! -f $gro_file ]; then
	echo "$gro_file doesn't exist"
	exit 1
fi

cp posre.itp posre.itp_bkp

IFS=',' read -r -a resi_arr <<< `cat $file_name`
for resi in "${resi_arr[@]}"
do
	for i in `grep ^[[:space:]]*"$resi"[A-Z] $gro_file | awk '{print $3}'`
	do
		echo $i
		#grep [[:blank:]]*"$i" posre.itp
		sed -i "/[[:space:]]*$i/d" posre.itp
	done
done
