#!/bin/bash

protname="$1"

if [ -z $protname ]; then
	echo "Need 1st argument"
	exit -1
fi

alreadyrun="../proteinsparams/listoffiles.txt"

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
       protfolder=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`
       proteinname=`echo $protfolder`

       PROTFOLDER="../protein"
       paramfile=$PROTFOLDER"/"$proteinname"/"$proteinname".xml"
       seq=`xmllint --xpath "string(//seq_file)" "$paramfile"`
       aco=`xmllint --xpath "string(//ang_file)" "$paramfile"`

       upls_1=`xmllint --xpath "string(//upl_file)" "$paramfile"`
       upls_hbo=`xmllint --xpath "string(//hbond_file)" "$paramfile"`       
 
       if [ ! -z $upls_hbo ]; then
	    upls=$upls_1','$upls_hbo
       else
	    upls=$upls_1
       fi

       IFS=', ' read -r -a array <<< "$upls"

       new_file=../protein/$proteinname/$proteinname"_concatupl.upl"
       if [ -f $new_file ]; then
	    echo "Deleting old concatenated upl file: $new_file"
	    rm -f $new_file
       fi
       for files in "${array[@]}"
       do
     	    cat ../protein/$proteinname/$files >> $new_file
       done
      
       if [ -f ../protein/$proteinname/$seq ]; then 
	    if [ -f $new_file ]; then
		if [ ! -z $aco ]; then
			./makeXplorUplAco.sh ../protein/$proteinname/$seq $new_file ../protein/$proteinname/$aco
		else
			./makeXplorUplAco.sh ../protein/$proteinname/$seq $new_file
		fi
	    else
		    echo "could not find upl file"
	    fi
       else
	    echo "could not find seq file"
       fi
       mv $proteinname"_concatupl.upl.txt" ../protein/$proteinname/. 
else
       paramfile=../protein/$protname/$protname".xml"
       seq=`xmllint --xpath "string(//seq_file)" "$paramfile"`
       aco=`xmllint --xpath "string(//ang_file)" "$paramfile"`
       upls_1=`xmllint --xpath "string(//upl_file)" "$paramfile"`
       upls_hbo=`xmllint --xpath "string(//hbond_file)" "$paramfile"`       
 
       if [ ! -z $upls_hbo ]; then
	    upls=$upls_1','$upls_hbo
       else
	    upls=$upls_1
       fi

       IFS=', ' read -r -a array <<< "$upls"

       new_file=../protein/$protname/$protname"_concatupl.upl"
       if [ -f $new_file ]; then
	    echo "Deleting old concatenated upl file: $new_file"
	    rm -f $new_file
       fi
       for files in "${array[@]}"
       do
     	    cat ../protein/$protname/$files >> $new_file
       done

       if [ -f ../protein/$protname/$seq ]; then 
	    if [ -f $new_file ]; then
		if [ ! -z $aco ]; then
			./makeXplorUplAco.sh ../protein/$protname/$seq $new_file ../protein/$protname/$aco
		else
			./makeXplorUplAco.sh ../protein/$protname/$seq $new_file
		fi
	    else
		    echo "could not find upl file"
	    fi
       else
	    echo "could not find seq file"
       fi
       mv $protname"_concatupl.upl.txt" ../protein/$protname/.
fi

