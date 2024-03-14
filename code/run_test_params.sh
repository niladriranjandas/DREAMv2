#!/bin/bash
date

# - convert upl file to tbl file - #
#../packages/PdbStat510_20130909/convertUplToTbl_params.sh "$1"

# - convert aco file for new cases - #
./callConvertAco_params.sh "$1"

# - check all files #
./checkAllFiles_params.sh "$1"

# - begin the algo - #

#./driverCode.sh "$1"
./genRunScripts_bckend_parstage1_parstage2.sh "$1"

./callBadGaps_params.sh "$1"

# - move final files to protein folder - #
./gatherFiles_params.sh "$1"

# - clean up - #
./cleanUpFiles_params.sh "$1"
./makeRamaChart_params.sh "$1"

# - dist viols -#
../packages/find_viols/calcViols_params.sh "$1"
../packages/find_viols/plotCalcViols_params.sh "$1"

# - progress algo -#
./front_end/sendErrorMsg_params.sh "$1"

if [ -f ../protein/$1"/progress_msgs/sendfiles_list.txt" ]; then
	email="$2"
	if [ ! -z "$email" ]; then
	        #python front_end/mailsendErrorMsg_params_v3.py "$2" "niladrid@iisc.ac.in" "../protein/"$1"/progress_msgs/sendfiles_list.txt" "../protein/$1/progress_msgs/sendfiles_list" 3
		python3 front_end/mailsendErrorMsg_params_v3.py "$2" "niladrid@iisc.ac.in" "../protein/"$1"/progress_msgs/sendfiles_list.txt" "../protein/$1/progress_msgs/sendfiles_list" 3
	else
		#python front_end/mailsendErrorMsg_params_v3.py "abc" "niladrid@iisc.ac.in" "../protein/"$1"/progress_msgs/sendfiles_list.txt" "../protein/$1/progress_msgs/sendfiles_list" 3
		python3 front_end/mailsendErrorMsg_params_v3.py "abc" "niladrid@iisc.ac.in" "../protein/"$1"/progress_msgs/sendfiles_list.txt" "../protein/$1/progress_msgs/sendfiles_list" 3
	fi
        if [ -f "../protein/"$1"/progress_msgs/sendfiles_list.abbr.txt" ]; then
                soffice --convert-to pdf "../protein/"$1"/progress_msgs/sendfiles_list.abbr.txt" --outdir "../protein/"$1"/progress_msgs/"
                if [ -f "../protein/"$1"/progress_msgs/sendfiles_list.png" ]; then
                        convert "../protein/"$1"/progress_msgs/sendfiles_list.png" "../protein/"$1"/progress_msgs/sendfiles_list.abbr.pdf" "../protein/"$1"/progress_msgs/sendfiles_list.runstatus.pdf"
                else
                        echo "Error: ../protein/"$1"/progress_msgs/sendfiles_list.abbr.txt NOT found"
                fi
        else
                echo "../protein/"$1"/progress_msgs/sendfiles_list.png NOT found"
        fi
else
        echo "Error: ../protein/"$1"/progress_msgs/sendfiles_list.txt not found"
fi

# - send email - #
#if [ ! -z "$2" ]; then
#	#./callSendEmail.sh "$1" "$2"
#	./joinNsend_params.sh "$1" "$2"
#fi
./joinNsend_params.sh "$1" "$2"    #for creating files even if no email is there

# move the files for report
mv "$1"* "../protein/"$1"/progress_msgs/."

./transferToFrontend_params.sh "$1"

date
