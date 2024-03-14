#!/bin/bash
date

# - convert upl file to tbl file - #
#../packages/PdbStat510_20130909/convertUplToTbl.sh "$1"

# - convert aco file for new cases - #
./callConvertAco.sh "$1"

# - check if files are present - #
./checkAllFiles.sh "$1"

# - begin algo - #
./driverCode.sh "$1"
./callBadGaps.sh "$1"

# - move final files to protein folder - #
./gatherFiles.sh "$1"

# - clean up - #
./cleanUpFiles.sh "$1"
./makeRamaChart.sh "$1"

# - dist viols -#
../packages/find_viols/calcViols.sh "$1"
../packages/find_viols/plotCalcViols.sh "$1"

# - progress algo -#
./front_end/sendErrorMsg.sh "$1"
email="$2"
if [ ! -z "$email" ]; then
	./front_end/callmailsendErrorMsg_param_v2.sh "$2" "niladrid@iisc.ac.in" "$1"
else
	./front_end/callmailsendErrorMsg_param_v2.sh "abc" "niladrid@iisc.ac.in" "$1"
fi

# - send email - #
#if [ ! -z "$2" ]; then
#	#./callSendEmail.sh "$1" "$2"
#	./joinNsend.sh "$1" "$2"
#fi
./joinNsend.sh "$1" "$2"    #for creating files even if no email is there

# move the files for report
./front_end/moveReportFiles.sh "$1"

#./transferToFrontend.sh "$1"

date
