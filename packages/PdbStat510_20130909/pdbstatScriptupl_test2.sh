#!/bin/bash

seqfile="$1"
uplfile="$2"

destdir=`dirname $uplfile`

firstresi=`head -1 $seqfile | awk '{print $2}'`

noe_tbl=`basename $uplfile`
noe_tbl_tmp=${noe_tbl::-4}
noe_tbl_name=$noe_tbl_tmp"_noe.tbl"

echo "rea sequence $seqfile" > $destdir/convertUpl.txt
echo "$firstresi" >> $destdir/convertUpl.txt
echo "rea cons cyana $uplfile" >> $destdir/convertUpl.txt
echo "write" >> $destdir/convertUpl.txt
echo "$destdir/$noe_tbl_name" >> $destdir/convertUpl.txt
echo "cons" >> $destdir/convertUpl.txt
echo "xplor" >> $destdir/convertUpl.txt
echo "5" >> $destdir/convertUpl.txt
echo "vdw" >> $destdir/convertUpl.txt


