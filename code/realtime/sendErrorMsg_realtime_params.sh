#!/bin/bash

protname="$1"

datafolder="../protein"
emailfileslist="sendfiles_list_realtime.txt"
errfolder="realtime"
emailfolder="$errfolder"


        #paramfile=`egrep ^$protname $alreadyrun | awk -F',' '{print $2}'`
        protfolder="$protname"

        emailfolderpath=$datafolder/$protfolder/$emailfolder
        sendfile=$emailfolderpath/$emailfileslist

        [ ! -d "$emailfolderpath" ] && mkdir "$emailfolderpath"

        #if [ ! -f "$sendfile" ] ; then  # error
		errfile=$errfolder"/"$protname"_err_msg.txt"
		
                all_param_log=`ls "$protname"_minimal_log.txt 2> /dev/null`
                #paramfiles=`ls "$protname"_param*_minimal_log.txt 2> /dev/null`
                anclocfile=`ls "$protname"_?_ancloc_refine*minimal_log.txt 2> /dev/null`

                echo "$all_param_log"
                #echo "$paramfiles"
                echo "$anclocfile"


                if [ ! -z "$all_param_log" ]; then
			            paramparse=`grep -oP "(?<=Parsing parameter file).*" "$all_param_log"`
			            div_preloc=`grep -oP "(?<=Run divide into fragments and model).*" "$all_param_log"`
			            consolidate=`grep -oP "(?<=Consolidating the fragmented model).*" "$all_param_log"`
			            anchorloc=`grep -oP "(?<=Anchored localization).*" "$all_param_log"`
			            endancloc=`grep -oP "(?<=Ending anchored locatlization).*" "$all_param_log"`
            
			            paramparse_=`echo $paramparse|awk '{print $NF}'`
			            div_preloc_=`echo $div_preloc|awk '{print $NF}'`
			            consolidate_=`echo $consolidate|awk '{print $NF}'`
			            anchorloc_=`echo $anchorloc|awk '{print $NF}'`
			            endancloc_=`echo $endancloc|awk '{print $NF}'`
            
			            if [[ $paramparse_ -ne 0 ]]; then
			            	echo "Error in parsing the input parameter files." > $sendfile
			            else
			            	echo "Success in parsing the input parameter files." > $sendfile
			            fi
			            
			            if [[ $div_preloc_ -ne 0 ]]; then
			            	echo "Error in dividing and individual modeling section." >> $sendfile
			            else
			            	echo "Successfully divided and individual modeling section." >> $sendfile
			            fi
            
			            if [[ $consolidate_ -ne 0 ]]; then
			            	echo "Error in consolidating the divided parts." >> $sendfile
			            else
			            	echo "Successfully consolidating the divided parts." >> $sendfile
			            fi
            
			            if [[ $anchorloc_ -ne 0 ]]; then
			            	echo "Error in modelling the gaps." >> $sendfile
			            else
			            	echo "Successfully modelled the gaps." >> $sendfile
			            fi
            
			            if [[ $endancloc_ -ne 0 ]]; then
			            	echo "Error in modelling the gaps post processing steps." >> $sendfile
			            else
			            	echo "Successfully modelled the gaps post processing steps." >> $sendfile
			            fi
                else
                	echo "ongoing-1"
                fi

                echo "############################" >> $sendfile
                ###-------------------------------------------------------#####
                if [ ! -z "$anclocfile" ]; then
                    for anclocfiles in `ls $anclocfile`
                    do
                    	echo " " >> $sendfile
                    	if [ ! -z "$anclocfiles" ]; then
		            	ap_correct_modeller=`grep -oP "(?<=Gap-correct modeller).*" "$anclocfiles"`
		            	Gap_correct_choosing_gaps=`grep -oP "(?<=Gap-correct choosing gaps).*" "$anclocfiles"`
		            	Gap_correct_put_together_chosen_gaps=`grep -oP "(?<=Gap-correct put together chosen gaps).*" "$anclocfiles"`
		            	Gap_correct_putting_togther_geometric_refinement=`grep -oP "(?<=Gap-correct putting togther geometric refinement).*" "$anclocfiles"`
		            	EM=`grep -oP "(?<=EM).*" "$anclocfiles"`
		            	water_refine_Convert_file=`grep -oP "(?<=water-refine Convert file).*" "$anclocfiles"`
		            	water_refine_running=`grep -oP "(?<=water-refine running).*" "$anclocfiles"`
            
		            	ap_correct_modeller_=`echo $ap_correct_modeller|awk '{print $NF}'`
		            	Gap_correct_choosing_gaps_=`echo $Gap_correct_choosing_gaps|awk '{print $NF}'`
		            	Gap_correct_put_together_chosen_gaps_=`echo $Gap_correct_put_together_chosen_gaps|awk '{print $NF}'`
		            	Gap_correct_putting_togther_geometric_refinement_=`echo $Gap_correct_putting_togther_geometric_refinement|awk '{print $NF}'`
		            	EM_=`echo $EM|awk '{print $NF}'`
		            	water_refine_Convert_file_=`echo $water_refine_Convert_file|awk '{print $NF}'`
		            	water_refine_running_=`echo $water_refine_running|awk '{print $NF}'`
            
		            	if [[ $ap_correct_modeller_ -ne 0 ]]; then
		            		echo "Error in starting structure generation for gaps." >> $sendfile
		            	else
		            		echo "Successfully generated starting structure generation for gaps." >> $sendfile
		            	fi
		            	
		            	if [[ $Gap_correct_choosing_gaps_ -ne 0 ]]; then
		            		echo "Error in choosing starting structures for gap regions." >> $sendfile
		            	else
		            		echo "Success in choosing starting structures for gap regions." >> $sendfile
		            	fi
            
		            	if [[ $Gap_correct_put_together_chosen_gaps_ -ne 0 ]]; then
		            		echo "Error in putting together choosing starting structures for gap regions." >> $sendfile
		            	else
		            		echo "Successfully put together choosing starting structures for gap regions." >> $sendfile
		            	fi
            
		            	if [[ $Gap_correct_putting_togther_geometric_refinement_ -ne 0 ]]; then
		            		echo "Error in putting together chosen structrue for the gaps." >> $sendfile
		            	else
		            		echo "Successfully put together chosen structrue for the gaps." >> $sendfile
		            	fi
            
		            	if [[ $EM_ -ne 0 ]]; then
		            		echo "Error in post processing EM step." >> $sendfile
		            	else
		            		echo "Successfully done post processing EM step." >> $sendfile
		            	fi
            
		            	if [[ $water_refine_Convert_file_ -ne 0 ]]; then
		            		echo "Error in convertion for water refinement step." >> $sendfile
		            	else
		            		echo "Successfully converted for water refinement step." >> $sendfile
		            	fi
            
            
		            	if [[ $water_refine_running_ -ne 0 ]]; then
		            		echo "Error in running water refinement step in post processing." >> $sendfile
		            	else
		            		echo "Successfully ran water refinement step in post processing." >> $sendfile
		            	fi		        	
                    	else
                    		echo "ongoing-2"
                    	fi
                    done
                fi
        
