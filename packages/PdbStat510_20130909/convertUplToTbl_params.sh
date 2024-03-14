#!/bin/bash

protname="$1"

#aco_convert_cmd="python3 ../packages/PdbStat510_20130909/acoCyanaToXplor.py"
uplTotbl_cmd="../packages/PdbStat510_20130909/pdbstat -s < "

#alreadyrun="../proteinsparams/listoffiles.txt"

#found=`egrep ^$protname, $alreadyrun`
#if [ ! -z "$found" ]; then
#        echo " "
#else
	PROTFOLDER="../protein"
	xmlprotip=$protname".xml"

        seq_cyana_file=`xmllint --xpath "string(//seq_file)" "$xmlprotip"`
 	upl_cyana_file=`xmllint --xpath "string(//upl_file)" "$xmlprotip"`
	upl_xplor_file=`xmllint --xpath "string(//upl_xplor)" "$xmlprotip"`
	
	scriptfile=$protname"_convertupl.txt"
	
	# prepare the script
	firstresi=`head -1 $seq_cyana_file | awk '{print $2}'`
	echo "rea sequence $seq_cyana_file" > $scriptfile
	echo "$firstresi" >> $scriptfile
	echo "rea cons cyana $upl_cyana_file" >> $scriptfile
	echo "write" >> $scriptfile
	echo "$upl_xplor_file" >> $scriptfile
	echo "cons" >> $scriptfile
	echo "xplor" >> $scriptfile
	echo "5" >> $scriptfile
	echo "vdw" >> $scriptfile
	echo "quit" >> $scriptfile
	
	# run pdbstat
	$uplTotbl_cmd $scriptfile
        /data2/nmr/our_algo_final/packages/PdbStat510_20130909/pdbstat -s < "$scriptfile"
	
#fi	

	
	


