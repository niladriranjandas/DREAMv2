#!/bin/bash

filename="$1"

if [ ! -f "$filename" ]; then
	echo "$filaname NOT found"
	exit -1
fi

pdbname=`basename "$filename"`
folder=`dirname "$filaname"`
phipsifile=$folder"/"$pdbname".tsv"
./pdbtorsions $1 > $phipsifile

if [ ! -f "$phipsifile" ]; then
	echo "$phipsifile NOT found"
	exit -1
fi

pngfile=$folder"/"$pdbname".png"
Rscript draw_rama.r "$phipsifile" "$pngfile"
