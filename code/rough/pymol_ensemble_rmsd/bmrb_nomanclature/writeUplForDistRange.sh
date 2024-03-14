#!/bin/bash

#*** write the upl's corr to a range of distances****
# e.g., strong <= 2.7 A, medium >2.7 <= 3.3, weak >3.3 <=5
#     ./writeUplForDistRange.sh 5x1x.dist.1.upl 0,3.5   
# usage: ./writeUplForRange.sh <upl_file> <range=lo1,up1:lo2,upl2>
##########################################

usage(){
	echo "usage: /writeUplForDistRange.sh <upl_file> <lo,up>"
}

if [[ $# -ne 2 ]]; then
	usage
	exit 1
fi

upl_file="$1"
range="$2"
if [ ! -f $upl_file ]; then
	echo "upl file $upl_file not found"
	exit 1
fi 

######################################################

  lo1=`echo $range | awk -F',' '{print $1}'`
  up1=`echo $range | awk -F',' '{print $2}'`

  if [ $(bc <<< "$lo1<$up1") -eq 1 ]; then
	while IFS='' read -r line || [[ -n "$line" ]]; do
		dist_ij=`echo $line | awk '{print $7}'`
                #echo $dist_ij
		if [ $(bc <<< "$dist_ij>$lo1 && $dist_ij<=$up1") -eq 1 ]; then
			echo "$line"
		fi
	done < $upl_file
  fi
