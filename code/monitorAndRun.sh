#!/bin/bash

monitorfolder="toberun"
movefolder="ran"

monitorfolder_params="toberun_withparam"

DELAY=5 #10

while [[ true ]];
do
  sleep $DELAY
  echo "checking..."  
  folders=`ls $monitorfolder`
  if [ ! -z "$folders" ]; then
     echo "FOUND $folders"
     for folder_i in `echo $folders`
     do
	echo $folder_i
         if [ -f $monitorfolder/$folder_i ]; then
		# check if its a zipped file
		extension=${folder_i: -4}		
		if [[ "$extension" == ".zip" ]]; then
			folder_i_name=${folder_i::-4}
			unzip $monitorfolder/$folder_i -d ../protein/$folder_i_name    # after unzip update this create folder within a folder
			#unzip $monitorfolder/$folder_i -d ../protein
			if [ -d "../protein/$folder_i_name" ]; then
				xmlfile="../protein/"$folder_i_name/$folder_i_name".xml"
				email=`xmllint --xpath "string(//email)" "$xmlfile"`
                                      ######## pdb stat give error otherwise #######
	                              	  curr_dir=`pwd`
                                      cd ../protein/$folder_i_name
                                      #/data2/nmr/our_algo_final/packages/PdbStat510_20130909/convertUplToTbl.sh $folder_i_name
                                      /data2/nmr/our_algo_final/packages/convert_upl/convertUplToTbl.sh $folder_i_name
                                      cd $curr_dir 
                                      #############################################
				echo "ts ./run_test.sh $folder_i_name $email"
                                ts ./run_test.sh $folder_i_name $email
				#ts ./driverCode.sh $folder_i_name
				#ts ./callBadGaps.sh $folder_i_name
				#ts ./gatherFiles.sh $folder_i_name
				#ts ./cleanUpFiles.sh $folder_i_name
				#ts ./makeRamaChart.sh $folder_i_name
				#if [[ ! -z $email ]]; then
				#	ts ./callSendEmail.sh $folder_i_name $email
				#fi
			fi
		else
                        name=`cat $monitorfolder/$folder_i | awk -F',' '{print $1}'`
                 	email=`cat $monitorfolder/$folder_i | awk -F',' '{print $2}'`
                        # in this case running pdbstat is not required
			#../packages/PdbStat510_20130909/convertUplToTbl.sh $name
	         	echo "ts ./run_test.sh $name $email"
                	ts ./run_test.sh $name $email
		fi
         elif [ -d $monitorfolder/$folder_i ]; then
		echo "move folders and call progs"
	 fi
         mv "$monitorfolder/$folder_i" "$movefolder"/.
     done
  fi
  
  # ---------- for params mode --------- #
  folders_params=`ls $monitorfolder_params`
  if [ ! -z "$folders_params" ]; then
     echo "FOUND $folders_params"
     for folder_i_params in `echo $folders_params`
     do
	echo $folder_i_params
        if [ -f $monitorfolder_params/$folder_i_params ]; then
			# check if its a zipped file
			extension=${folder_i_params: -4}		
			if [[ "$extension" == ".zip" ]]; then
				folder_i_name=${folder_i_params::-4}
				unzip $monitorfolder_params/$folder_i_params -d ../protein/$folder_i_name   ### due to uzip update this creates folder within folder
				#unzip $monitorfolder_params/$folder_i -d ../protein
				if [ -d "../protein/$folder_i_name" ]; then
					xmlfile="../protein/"$folder_i_name/$folder_i_name".xml"
					email=`xmllint --xpath "string(//email)" "$xmlfile"`
                                      ######## pdb stat give error otherwise #######
	                              	  curr_dir=`pwd`
                                      cd ../protein/$folder_i_name
                                      #/data2/nmr/our_algo_final/packages/PdbStat510_20130909/convertUplToTbl_params.sh $folder_i_name
                                      /data2/nmr/our_algo_final/packages/convert_upl/convertUplToTbl_params.sh $folder_i_name
                                      cd $curr_dir 
                                      #############################################
					echo "ts ./run_test_params.sh $folder_i_name $email"
					ts ./run_test_params.sh $folder_i_name $email
					#ts ./driverCode.sh $folder_i_name
					#ts ./callBadGaps.sh $folder_i_name
					#ts ./gatherFiles.sh $folder_i_name
					#ts ./cleanUpFiles.sh $folder_i_name
					#ts ./makeRamaChart.sh $folder_i_name
					#if [[ ! -z $email ]]; then
					#	ts ./callSendEmail.sh $folder_i_name $email
					#fi
				fi
			fi
        elif [ -d $monitorfolder_params/$folder_i_params ]; then
			echo "move folders and call progs"
 	fi
        mv "$monitorfolder_params/$folder_i_params" "$movefolder"/.
     done
  fi
  # ------------------------------------ #

done

