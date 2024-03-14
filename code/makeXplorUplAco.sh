#!/bin/bash

seq_file="$1"
upl_file="$2"
aco_file="$3"

PDBSTAT_UPL="/data2/nmr/our_algo/packages/PdbStat510_20130909/pdbstatScriptupl.sh"
PDBSTAT_cmd="/data2/nmr/our_algo/packages/PdbStat510_20130909/pdbstat"

ACOCONVERT="python /data2/nmr/our_algo/packages/PdbStat510_20130909/acoCyanaToXplor.py"

	# upl file #
	destdir=`dirname $upl_file`
	noe_tbl=`basename $upl_file`
	noe_tbl_tmp=${noe_tbl::-4}
	noe_tbl_name=$noe_tbl_tmp"_noe.tbl"
	opp_upl_file=$destdir/$noe_tbl_name
	opp_upl_script=$noe_tbl".txt"

	$PDBSTAT_UPL $seq_file $upl_file > $opp_upl_script
	$PDBSTAT_cmd -s < $opp_upl_script

	# aco file #
	if [ ! -z $aco_file ]; then
		destacodir=`dirname $aco_file`
		aco_tbl=`basename $aco_file`
		aco_tbl_tmp=${aco_tbl::-4}
		aco_tbl_name=$aco_tbl_tmp"_aco.tbl"
		opp_aco_file=$destacodir/$aco_tbl_name
		$ACOCONVERT $aco_file > $opp_aco_file
	fi

