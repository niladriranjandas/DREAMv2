#!/bin/bash

 protname="$1"  # assumes folder by the name of the protein
 
 PROTFOLDER="../protein" 
 xmlprotip=$PROTFOLDER"/"$protname"/"$protname".xml"
 paramsfolder="paramwisetrial/paramfiles" 
 logfile=$protname"_log_all.log.txt"

 paramfilearg=""
 for parami in `ls $paramsfolder/`
 do
 	paramname=`echo $parami | awk -F'.xml' '{print $1}'`
        echo "$paramname"
 	name=$protname"_"$paramname
 	mkdir $PROTFOLDER/$name
 	upls=`xmllint --xpath "string(//upl_file)" "$xmlprotip"`
 	seq=`xmllint --xpath "string(//seq_file)" "$xmlprotip"` 	
 	aco=`xmllint --xpath "string(//ang_file)" "$xmlprotip"`
 	hbond=`xmllint --xpath "string(//hbond_file)" "$xmlprotip"`

        echo "<set_up>" > $PROTFOLDER/$name/$name".xml"
        cat $xmlprotip >> $PROTFOLDER/$name/$name".xml"
        cat $paramsfolder/$parami >> $PROTFOLDER/$name/$name".xml"
        echo "</set_up>" >> $PROTFOLDER/$name/$name".xml"

#echo "$name 1"
xmllint --shell $PROTFOLDER/$name/$name".xml"  << EOF
cd /set_up/inputs/protein_name
set $name
save
EOF
	
 	IFS=',' read -r -a array <<< "$upls"
 	for element in "${array[@]}"
	do
	    cp $PROTFOLDER/$protname/$element $PROTFOLDER/$name/.
	done
 	cp $PROTFOLDER/$protname/$seq $PROTFOLDER/$name/.
 	#if [[ ! -z $lo ]]; then
	# 	cp $PROTFOLDER/$protname/$lo $PROTFOLDER/$name/.
	#fi
	if [[ ! -z $aco ]]; then
	 	cp $PROTFOLDER/$protname/$aco $PROTFOLDER/$name/.
 	fi
 	if [[ ! -z $hbond ]]; then
	 	cp $PROTFOLDER/$protname/$hbond $PROTFOLDER/$name/.
	fi
 	
 	paramfilearg=$paramfilearg","$name
 done
 
 paramfilearg_="${paramfilearg:1}"
 printf "\nparamwisetrial/genRunScipts_bckend_predictDrive.sh %s %s > paramwisetrial/%s\n" "$protname" "$paramfilearg_" "$logfile"
 #paramwisetrial/genRunScipts_bckend_predictDrive.sh $protname $paramfilearg_ > paramwisetrial/$logfile
 paramwisetrial/genRunScipts_bckend_predictDrive_par.sh $protname $paramfilearg_ > paramwisetrial/$logfile
