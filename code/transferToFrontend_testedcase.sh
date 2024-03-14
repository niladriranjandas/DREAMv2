#!/bin/bash

protname="$1"

uname="dream"
uip="10.24.34.4"
upass="niladri@123"
transferpath1="/var/www/html/DREAM/pdbs"
transferpath2="/var/www/html/DREAM/userpdbs"

sendfiles="sendfiles/sendfiles_list.txt"

alreadyrun="../proteinsparams/listoffiles.txt"
datafolder="../protein"

summaryfile="front_end/summary.php"

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
        protfolder=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`
        listoffiles="$datafolder/$protfolder/$sendfiles"

        if [ -f $listoffiles ]; then
        	# -- make directory remotely -- #
        	destdir="$transferpath1/$protname"
        	sshpass -p "$upass" ssh "$uname"@"$uip" "mkdir -p $destdir" 

        	# -- transfer files -- #
        	while IFS= read -r line || [[ -n "$line" ]]; do
			    sshpass -p "$upass" scp "$line" "$uname@$uip:$destdir/."
		done < "$listoffiles"
		sshpass -p "$upass" scp "$summaryfile" "$uname@$uip:$destdir/."
        else
        	echo "$listoffiles NOT found"
        	error -1
        fi
else
	protfolder=$protname"_allparam"
        listoffiles="$datafolder/$protfolder/$sendfiles"

        if [ -f $listoffiles ]; then
        	# -- make directory remotely -- #
        	destdir="$transferpath2/$protname"
        	sshpass -p "$upass" ssh "$uname"@"$uip" "mkdir -p $destdir"

        	# -- transfer files -- #
        	while IFS= read -r line || [[ -n "$line" ]]; do
			    sshpass -p "$upass" scp "$line" "$uname@$uip:$destdir/."			    
		done < "$listoffiles"
		sshpass -p "$upass" scp "$summaryfile" "$uname@$uip:$destdir/."
        else
        	echo "$listoffiles NOT found"
        	error -1
        fi
fi		
