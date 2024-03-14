#!/bin/bash

protname="$1"

#alreadyrun="/data2/nmr/our_algo_final/proteinsparams/listoffiles.txt"

#found=`egrep ^$protname, $alreadyrun`
#if [ ! -z "$found" ]; then
#        echo " "
#else
        POTFOLDER="../protein"
        xmlprotip=$protname".xml"

        seq_cyana_file=`xmllint --xpath "string(//seq_file)" "$xmlprotip"`
        upl_cyana_file=`xmllint --xpath "string(//upl_file)" "$xmlprotip"`
        upl_xplor_file=`xmllint --xpath "string(//upl_xplor)" "$xmlprotip"`

	python /data2/nmr/our_algo_final/packages/convert_upl/upl_tbl_conversion_script_modified.py $upl_cyana_file $upl_xplor_file

#fi       
