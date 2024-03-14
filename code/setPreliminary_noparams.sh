#!/bin/bash

help()
{
    echo "Usage: weather [ -c | --city1 ]
               [ -d | --city2 ]
               [ -h | --help  ]"
    exit 2
}

SHORT=n:,s:,u:,a:,h
LONG=name:,seq:,upl:,ang:,help
OPTS=$(getopt -a -n setPreliminary_noparams --options $SHORT --longoptions $LONG -- "$@")

VALID_ARGUMENTS=$# # Returns the count of arguments that are in short or long options

if [ "$VALID_ARGUMENTS" -eq 0 ]; then
  help
fi

eval set -- "$OPTS"

while :
do
  case "$1" in
    -n | --name )
      name="$2"
      shift 2
      ;;
    -s | --seq )
      seq="$2"
      shift 2
      ;;
    -u | --upl )
       upl="$2"
       shift 2
       ;;
    -a | --ang )
       ang="$2"
       shift 2
       ;;
    -h | --help)
      help
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      help
      ;;
  esac
done

######## Define #########
PROTEINFOLDER="../protein"
randno="$RANDOM"
######## check if inputs are there or not #######
	if [ -z "$name" ]; then
		echo "name of protein cannot be ommitted"
		exit -1
	fi

	if [ -z "$seq" ]; then
		echo "seq is needed"
		exit -1
	else
		if [ ! -f "$seq" ]; then
			echo "$seq NOT file"
			exit -1
		fi
	fi

	if [ -z "$upl" ]; then
		echo "upl is needed"
		exit -1
	else
		if [ ! -f "$upl" ]; then
			echo "$upl NOT file"
			exit -1
		fi
	fi
	
	if [ ! -z "$ang" ]; then
		if [ ! -f "$ang" ]; then
			echo "$ang NOT file"
			exit -1
		fi
	fi

####################################################
	
        foldername=$PROTEINFOLDER/$name"_"$randno
	if [ -d $foldername ]; then
		echo "DELETE or rename $foldername before running"
		exit -1
	fi
        mkdir $foldername

	cp $seq $foldername/.
	cp $upl $foldername/.

	seqfilename=`basename $seq` 
	uplname=`basename $upl`

	../packages/PdbStat510_20130909/pdbstatScriptupl_test2.sh "$foldername/$seqfilename"  "$foldername/$uplname"
	../packages/PdbStat510_20130909/pdbstat -s < "$foldername/convertUpl.txt"
	uplfile_=${uplfile::-4}
	xplorupl=$uplfile"_noe.tbl"
	#cp $xplorupl $foldername/.

	if [ ! -z "$ang"  ]; then
		echo "abc $ang"
		cp "$ang" "$foldername"/.
		angname=`basename $ang`
		ang_=${angname::-4}
		xplorang=$ang_".tbl"
		python3 ../packages/PdbStat510_20130909/acoCyanaToXplor.py $ang > "$foldername/$xplorang"
	else
		angname=""
		xplorang=""
	fi

	xmlname=$name"_"$randno
	xmlfile=$foldername"/"$xmlname".xml"
####################################################
	
	printf "\n\t<inputs>" > $xmlfile
	printf "\n\t\t<protein_name>%s_%s</protein_name>" "$name" "$randno" >> $xmlfile
	printf "\n\t\t<seq_file>%s</seq_file>" "$seqfilename"  >> $xmlfile
	printf "\n\t\t<upl_file>%s</upl_file>" "$uplname"  >> $xmlfile
	printf "\n\t\t<upl_xplor>%s</upl_xplor>" "$xplorupl"  >> $xmlfile
	printf "\n\t\t<hbond_file></hbond_file>"  >> $xmlfile
	printf "\n\t\t<ang_file>%s</ang_file>" "$angname"  >> $xmlfile
	printf "\n\t\t<aco_xplor>%s</aco_xplor>" "$xplorang"  >> $xmlfile
	printf "\n\t\t<protein_path>../protein/%s_%s</protein_path>" "$name" "$randno"  >> $xmlfile
	printf "\n\t\t<email></email>"  >> $xmlfile
	printf "\n\t</inputs>"  >> $xmlfile

