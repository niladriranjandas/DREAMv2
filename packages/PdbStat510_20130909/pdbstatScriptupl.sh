#!/bin/bash

seqfile="$1"
uplfile="$2"

destdir=`dirname $uplfile`

firstresi=`head -1 $seqfile | awk '{print $2}'`

noe_tbl=`basename $uplfile`
noe_tbl_tmp=${noe_tbl::-4}
noe_tbl_name=$noe_tbl_tmp"_noe.tbl"

echo "rea sequence $seqfile"
echo "$firstresi"
echo "rea cons cyana $uplfile"
echo "write"
echo "$destdir/$noe_tbl_name"
echo "cons"
echo "xplor"
echo "5"
echo "vdw"


