#!/bin/bash

protname="$3"
tomail="$1"
frommail="$2"

alreadyrun="../proteinsparams/listoffiles.txt"

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
	
	protfolder=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`
	
	if [ -f ../protein/$protfolder/progress_msgs/sendfiles_list.txt ]; then	
		python3 front_end/mailsendErrorMsg_params_v3.py $tomail $frommail ../protein/$protfolder/progress_msgs/sendfiles_list.txt ../protein/$protfolder/progress_msgs/sendfiles_list 1

	        if [ -f "../protein/"$protfolder"/progress_msgs/sendfiles_list.abbr.txt" ]; then
        	        soffice --convert-to pdf "../protein/"$protfolder"/progress_msgs/sendfiles_list.abbr.txt" --outdir "../protein/"$protfolder"/progress_msgs/"
                	if [ -f "../protein/"$protfolder"/progress_msgs/sendfiles_list.png" ]; then
                        	convert "../protein/"$protfolder"/progress_msgs/sendfiles_list.png" "../protein/"$protfolder"/progress_msgs/sendfiles_list.abbr.pdf" "../protein/"$protfolder"/progress_msgs/sendfiles_list.runstatus.pdf"
	                else
        	                echo "Error: ../protein/"$protfolder"/progress_msgs/sendfiles_list.abbr.txt NOT found"
                	fi
	        else
        	        echo "../protein/"$protfolder"/progress_msgs/sendfiles_list.png NOT found"
	        fi
        else
		echo "Error: ../protein/$protfolder/progress_msgs/sendfiles_list.txt NOT found"
	fi
else
	if [ -f ../protein/$protname"_allparam/progress_msgs/sendfiles_list.txt" ]; then
        	python3 front_end/mailsendErrorMsg_params_v3.py "$tomail" "$frommail" ../protein/$protname"_allparam/progress_msgs/sendfiles_list.txt" ../protein/$protname"_allparam/progress_msgs/sendfiles_list" 2
	        if [ -f ../protein/$protname"_allparam/progress_msgs/sendfiles_list.abbr.txt" ]; then
        	        soffice --convert-to pdf ../protein/$protname"_allparam/progress_msgs/sendfiles_list.abbr.txt" --outdir ../protein/$protname"_allparam/progress_msgs/"
                	if [ -f ../protein/$protname"_allparam/progress_msgs/sendfiles_list.png" ]; then
                        	convert ../protein/$protname"_allparam/progress_msgs/sendfiles_list.png" ../protein/$protname"_allparam/progress_msgs/sendfiles_list.abbr.pdf" ../protein/$protname"_allparam/progress_msgs/sendfiles_list.runstatus.pdf"
	                else
        	                echo "Error: ../protein/"$protname"_allparam/progress_msgs/sendfiles_list.abbr.txt NOT found"
                	fi
	        else
        	        echo "../protein/"$protname"_allparam/progress_msgs/sendfiles_list.png NOT found"
	        fi
	else
        	echo "Error: ../protein/"$protname"_allparam/progress_msgs/sendfiles_list.txt not found"
	fi
fi
