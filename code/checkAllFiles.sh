#!/bin/bash

#Copyright (c) 2024 niladriranjandas
# check if all the files in xml are present


isFileEmpty () {
	filename="$1"
	
	if [ ! -f $filename ]; then
		return 0
	else
		if [ ! -s $filename ]; then
			return 0
		else
			return 1
		fi
	fi
}

protname="$1"
alreadyrun="../proteinsparams/listoffiles.txt"
ERR_FLAG=1

	found=`egrep ^$protname, $alreadyrun`

	if [ ! -z "$found" ]; then
		PROTFOLDER="../protein"
                # take backup since otherwise files are not transfered at the end
	        folder_name=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`
		xml_file_name=`egrep ^$protname, $alreadyrun | awk -F',' '{print $3}'`
		tmp=`date +'%d%m%Y%H%M%S'`
		new_folder_name=$folder_name"_"$tmp

		# copy only required files
		if [ -d "$PROTFOLDER/$folder_name" ]; then
			mv "$PROTFOLDER/$folder_name" "$PROTFOLDER/$new_folder_name"
			mkdir "$PROTFOLDER/$folder_name"
	                
			xmlprotip=$PROTFOLDER"/"$new_folder_name"/"$xml_file_name
			cp $xmlprotip $PROTFOLDER/$folder_name/.
			
			if [ -f "$xmlprotip" ]; then
        		        seq_cyana_file=`xmllint --xpath "string(//seq_file)" "$xmlprotip"`
				if [ ! -z "seq_cyana_file" ]; then
					if [ -f "$PROTFOLDER/$new_folder_name/$seq_cyana_file" ]; then
						cp $PROTFOLDER/$new_folder_name/$seq_cyana_file $PROTFOLDER/$folder_name/.
					else
						echo "ERROR $seq_cyana_file NOT found"
						ERR_FLAG=0
					fi
				else
					echo "Error seq file not present in file $xmlprotip"
					ERR_FLAG=0
				fi
				
                		upl_cyana_file=`xmllint --xpath "string(//upl_file)" "$xmlprotip"`
	                	upl_xplor_file=`xmllint --xpath "string(//upl_xplor)" "$xmlprotip"`
				if [ ! -z $upl_cyana_file ]; then
        	                        if [ -f "$PROTFOLDER/$new_folder_name/$upl_cyana_file" ]; then
                	                        cp $PROTFOLDER/$new_folder_name/$upl_cyana_file $PROTFOLDER/$folder_name/.
                        	        else
                                	        echo "ERROR $upl_cyana_file file NOT found"
                                        	ERR_FLAG=0
	                                fi
				else
					echo "Error upl cyana file not found $upl_cyana_file"
					ERR_FLAG=0
				fi
				if [ ! -z "$upl_xplor_file" ]; then
        	                        if [ -f "$PROTFOLDER/$new_folder_name/$upl_xplor_file" ]; then
                	                        cp $PROTFOLDER/$new_folder_name/$upl_xplor_file $PROTFOLDER/$folder_name/.
                        	        else
                                	        echo "ERROR $upl_xplor_file NOT found"
                                        	ERR_FLAG=0
	                                fi
				else
					echo "xmllint --xpath string(//upl_xplor) $xmlprotip"
					echo "Error in upl xplor file not found $upl_xplor_file"
					ERR_FLAG=0
				fi

	        	        hbond_cyana_file=`xmllint --xpath "string(//upl_xplor)" "$xmlprotip"`
				if [ ! -z "$hbond_cyana_file" ]; then
        	                        if [ -f "$PROTFOLDER/$new_folder_name/$hbond_cyana_file" ]; then
                	                        cp $PROTFOLDER/$new_folder_name/$hbond_cyana_file $PROTFOLDER/$folder_name/.
                        	        else
                                	        echo "ERROR $hbond_cyana_file file NOT found"
                                        	ERR_FLAG=0
	                                fi
				fi

        	        	aco_cyana_file=`xmllint --xpath "string(//ang_file)" "$xmlprotip"`
				if [ ! -z "$aco_cyana_file" ]; then
        	                        if [ -f "$PROTFOLDER/$new_folder_name/$aco_cyana_file" ]; then
                	                        cp $PROTFOLDER/$new_folder_name/$aco_cyana_file $PROTFOLDER/$folder_name/.
                        	        else
                                	        echo "ERROR $aco_cyana_file file NOT found"
                                        	ERR_FLAG=0
	                                fi
				fi
	        	        aco_xplor_file=`xmllint --xpath "string(//aco_xplor)" "$xmlprotip"`
				if [ ! -z "$aco_xplor_file" ]; then
        	                        if [ -f "$PROTFOLDER/$new_folder_name/$aco_xplor_file" ]; then
                	                        cp $PROTFOLDER/$new_folder_name/$aco_xplor_file $PROTFOLDER/$folder_name/.
                        	        else
                                	        echo "ERROR $aco_xplor_file file NOT found"
                                        	ERR_FLAG=0
	                                fi
				fi
			else
				echo "Error: XML file $xmlprotip NOT found"
				ERR_FLAG=0
			fi

			if [[ $ERR_FLAG -eq 0  ]]; then
				echo "ERROR: exiting due to missing/empty file"
	                        exit -1
			fi
			echo "--BEGIN--"
		else
			echo "Error protein folder $folder_name NOT found" 
			ERR_FLAG=0
			exit -1
		fi
	else
		PROTFOLDER="../protein"
		xmlprotip=$PROTFOLDER"/"$protname"/"$protname".xml"

		seq_cyana_file=`xmllint --xpath "string(//seq_file)" "$xmlprotip"`

		upl_cyana_file=`xmllint --xpath "string(//upl_file)" "$xmlprotip"`
		upl_xplor_file=`xmllint --xpath "string(//upl_xplor)" "$xmlprotip"`

		hbond_cyana_file=`xmllint --xpath "string(//upl_xplor)" "$xmlprotip"`		

		aco_cyana_file=`xmllint --xpath "string(//ang_file)" "$xmlprotip"`
		aco_xplor_file=`xmllint --xpath "string(//aco_xplor)" "$xmlprotip"`

		 if [ ! -z $seq_cyana_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$seq_cyana_file"
        	        if [ "$?" == 0 ]; then
        	        	echo "ERROR: $seq_cyana_file NOT found"
        	        	ERR_FLAG=0
        	        	#exit
        	        fi
		 fi
		 if [ ! -z $upl_cyana_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$upl_cyana_file"
        	        if [ "$?" == 0 ]; then
	        	        echo "ERROR: $upl_cyana_file NOT found"
	        	        ERR_FLAG=0
        	        	#exit
        	        fi
		 fi
		 if [ ! -z $upl_xplor_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$upl_xplor_file"
        	        if [ "$?" == 0 ]; then
	        	        echo "ERROR: $upl_xplor_file NOT found"
	        	        ERR_FLAG=0
        	        	#exit
        	        fi
		 fi
		 if [ ! -z $hbond_cyana_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$hbond_cyana_file"
        	        if [ "$?" == 0 ]; then
	        	        echo "ERROR: $hbond_cyana_file NOT found"
	        	        ERR_FLAG=0
        	        	#exit
        	        fi
		 fi	
		 if [ ! -z $aco_cyana_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$aco_cyana_file"
        	        if [ "$?" == 0 ]; then
	        	        echo "ERROR: $aco_cyana_file NOT found"
	        	        ERR_FLAG=0
        	        	#exit
        	        fi
		 fi		 	 		 		 
		 if [ ! -z $aco_xplor_file ]; then
        	        isFileEmpty "$PROTFOLDER/$protname/$aco_xplor_file"
        	        if [ "$?" == 0 ]; then
	        	        echo "ERROR: $aco_xplor_file NOT found"
	        	        ERR_FLAG=0
        	        	#exit
        	        fi
		 fi
		  
		 if [[ $ERR_FLAG -eq 0 ]]; then
		 	echo "ERROR: exiting due to missing/empty file"
		 	exit
		 fi
		 
		 echo "--BEGIN--"
	fi




