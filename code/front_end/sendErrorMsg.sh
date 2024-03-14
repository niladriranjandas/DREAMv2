#!/bin/bash

protname_inp="$1"

alreadyrun="../proteinsparams/listoffiles.txt"
datafolder="../protein"
emailfileslist="sendfiles_list.txt"
errfolder="progress_msgs"
emailfolder="$errfolder"

found=`egrep ^$protname_inp, $alreadyrun`
if [ ! -z "$found" ]; then
        #paramfile=`egrep ^$protname $alreadyrun | awk -F',' '{print $2}'`
        protfolder=`egrep ^$protname_inp, $alreadyrun | awk -F',' '{print $2}'`
        protname=$protfolder

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
			        		echo "Error in generated starting structure generation for gaps." >> $sendfile
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
        #fi  
        

else
	protname=$protname_inp
	protfolder=$protname"_allparam"
	
	emailfolderpath=$datafolder/$protfolder/$emailfolder
	sendfile=$emailfolderpath/$emailfileslist

	[ ! -d "$emailfolderpath" ] && mkdir "$emailfolderpath"
	
	#if [ ! -f "$sendfile" ] ; then	# error
		errfile=$errfolder"/"$protname"_err_msg.txt"
		
                all_param_log=`ls "$protname"_allparam_minimal_log.txt 2> /dev/null`
                paramfiles=`ls "$protname"_param?_minimal_log.txt 2> /dev/null`
                anclocfile=`ls "$protname"_param?_?_ancloc_refine_*_minimal_log.txt 2> /dev/null`

                if [ ! -z "$all_param_log" ]; then
                	#prepforpredict_files=`grep -oP "(?<=stage1 precursor for).*" $all_param_log`
                	echo "files for prediction:" > $sendfile
                	#for prepforpredict_files_i in `echo $prepforpredict_files`
                	grep -oP "(?<=stage1 precursor for).*" $all_param_log | while read -r prepforpredict_files_i;
                	do
                		prepforpredict_files_i_nm=`echo $prepforpredict_files_i | awk '{print $1}'`
                		prepforpredict_files_i_st=`echo $prepforpredict_files_i | awk '{print $2}'`

                		if [[ $prepforpredict_files_i_st -eq 0 ]]; then
                			prepforpredict_files_i_nm=${prepforpredict_files_i_nm::-1}
                			echo "	successfuly created preprediction file $prepforpredict_files_i_nm" >> $sendfile
                		else
                			prepforpredict_files_i_nm=${prepforpredict_files_i_nm::-1}
                			echo "	error in creating preprediction file $prepforpredict_files_i_nm" >> $sendfile
                		fi
                	done

                	echo "#####################################" >> $sendfile
                	dendrogram=`grep -oP "(?<=do dengrogram).*" $all_param_log`
                	dendrogram_nm=`echo $dendrogram | awk '{print $1}'`
                	dendrogram_st=`echo $dendrogram | awk '{print $2}'`
               		if [[ $dendrogram_st -eq 0 ]]; then
               			echo "successfuly created dendrogram file $dendrogram_nm" >> $sendfile
               		else
               			echo "error in creating dendrogram file $dendrogram_nm" >> $sendfile
               		fi       

           		echo "#####################################" >> $sendfile
           		gapcorrect_files=`grep -oP "(?<=stage1 precursor for).*" $all_param_log`
                	echo "files for gap correct:" >> $sendfile
                	#for gapcorrect_files_i in `echo $gapcorrect_files`
                	grep -oP "(?<=running gap correct for).*" $all_param_log | while read -r gapcorrect_files_i
                	do
                		gapcorrect_files_i_nm=`echo $gapcorrect_files_i | awk '{print $1}'`
                		gapcorrect_files_i_st=`echo $gapcorrect_files_i | awk '{print $2}'`

                		if [[ $gapcorrect_files_i_st -eq 0 ]]; then
                			gapcorrect_files_i_nm=${gapcorrect_files_i_nm::-1}
                			echo "	successfuly created gap correct file $gapcorrect_files_i_nm" >> $sendfile
                		else
                			gapcorrect_files_i_nm=${gapcorrect_files_i_nm::-1}
                			echo "	error in creating gap correct file $gapcorrect_files_i_nm" >> $sendfile
                		fi
                	done	

		fi

		echo "#####################################" >> $sendfile
		echo "#####################################" >> $sendfile
                if [ ! -z "$paramfiles" ]; then
                	echo "files for anchored localization:" >> $sendfile
                	for paramfiles_file in `echo $paramfiles`
                	do
                		echo "#####################################" >> $sendfile
                		runnning_modeller=`grep -oP "(?<=running modeller).*" "$paramfiles_file"`
                		runnning_modeller_nm=`echo $runnning_modeller | awk -F':' '{print $1}'`
                		runnning_modeller_st=`echo $runnning_modeller | awk -F':' '{print $2}'`
				
				if [[ $runnning_modeller_st -eq 0 ]]; then
					echo "	successfully created initial file for anchored localization $runnning_modeller_nm" >> $sendfile
				else
					echo "	erorr in creating initial file for anchored localization $runnning_modeller_nm" >> $sendfile
				fi


				anchored_loc_files=`grep -oP "(?<=anchored location for).*" "$paramfiles_file"`
				if [ ! -z "$anchored_loc_files" ]; then
					#for anchored_loc_files_i in `ls anchored_loc_files`
					grep -oP "(?<=anchored location for).*" "$paramfiles_file" | while read -r anchored_loc_files_i
					do
						running_anloc_nm=`echo $anchored_loc_files_i | awk -F':' '{print $1}'`
						running_anloc_st=`echo $anchored_loc_files_i | awk -F':' '{print $2}'`

                				if [[ $running_anloc_st -eq 0 ]]; then
                					echo "	successfuly created anchored localization file $running_anloc_nm" >> $sendfile
                				else
                					echo "	error running anchored localization file $running_anloc_nm" >> $sendfile
                				fi					
					done
				else
					echo "	anchored localization file create error" >> $sendfile
				fi

                		end_anloc=`grep -oP "(?<=ending anchored localization for).*" "$paramfiles_file"`
				if [ ! -z "$end_anloc" ]; then
                			end_anloc_nm=`echo $end_anloc | awk '{print $1}'`
                			end_anloc_st=`echo $end_anloc | awk '{print $2}'`
				
					if [[ $end_anloc_st -eq 0 ]]; then
						echo "	successfully created end of anchored localization $end_anloc_nm" >> $sendfile
					else
						echo "	erorr in creating end of anchored localization $end_anloc_nm" >> $sendfile
					fi				
				else
					echo "	end anclored localization error" >> $sendfile
				fi
                	done
		fi

                ###-------------------------------------------------------#####
      		echo "#####################################" >> $sendfile
		echo "#####################################" >> $sendfile

                for anclocfiles in `ls $anclocfile`
                do
                	echo "#####################################" >> $sendfile
                	if [ ! -z "$anclocfiles" ]; then
                		#echo "		$anclocfiles" >> $sendfile

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

				_ap_correct_modeller_=`echo $ap_correct_modeller|awk '{print $(NF-1)}'`
				_Gap_correct_choosing_gaps_=`echo $Gap_correct_choosing_gaps|awk '{print $(NF-1)}'`
				_Gap_correct_put_together_chosen_gaps_=`echo $Gap_correct_put_together_chosen_gaps|awk '{print $(NF-1)}'`
				_Gap_correct_putting_togther_geometric_refinement_=`echo $Gap_correct_putting_togther_geometric_refinement|awk '{print $(NF-1)}'`
				_EM_=`echo $EM|awk '{print $(NF-1)}'`
				_water_refine_Convert_file_=`echo $water_refine_Convert_file|awk '{print $(NF-1)}'`
				_water_refine_running_=`echo $water_refine_running|awk '{print $1 $2}'`				

				_ap_correct_modeller_=${_ap_correct_modeller_::-1}
				_Gap_correct_choosing_gaps_=${_Gap_correct_choosing_gaps_::-1}
				_Gap_correct_put_together_chosen_gaps_=${_Gap_correct_put_together_chosen_gaps_::-1}
				_Gap_correct_putting_togther_geometric_refinement_=${_Gap_correct_putting_togther_geometric_refinement_::-1}
				_EM_=${_EM_::-1}
				_water_refine_Convert_file_=${_water_refine_Convert_file_::-1}
				_water_refine_running_=${_water_refine_running_::-1}
        
		        	if [[ $ap_correct_modeller_ -ne 0 ]]; then
		        		echo "	Error in starting structure generation for gaps $_ap_correct_modeller_" >> $sendfile
		        	else
		        		echo "	Successfully generated starting structure generation for gaps $_ap_correct_modeller_" >> $sendfile
		        	fi
		        	
		        	if [[ $Gap_correct_choosing_gaps_ -ne 0 ]]; then
		        		echo "	Error in choosing of starting structures for gap regions $_Gap_correct_choosing_gaps_" >> $sendfile
		        	else
		        		echo "	Success in choosing of starting structures for gap regions $_Gap_correct_choosing_gaps_" >> $sendfile
		        	fi
        
		        	if [[ $Gap_correct_put_together_chosen_gaps_ -ne 0 ]]; then
		        		echo "	Error in putting together choosing starting structures for gap regions $_Gap_correct_put_together_chosen_gaps_" >> $sendfile
		        	else
		        		echo "	Successfully put together choosing starting structures for gap regions $_Gap_correct_put_together_chosen_gaps_" >> $sendfile
		        	fi
        
		        	if [[ $Gap_correct_putting_togther_geometric_refinement_ -ne 0 ]]; then
		        		echo "	Error in putting together chosen structrue for the gaps $_Gap_correct_putting_togther_geometric_refinement_" >> $sendfile
		        	else
		        		echo "	Successfully put together chosen structrue for the gaps $_Gap_correct_putting_togther_geometric_refinement_" >> $sendfile
		        	fi
        
		        	if [[ $EM_ -ne 0 ]]; then
		        		echo "	Error in post processing EM step $_EM_" >> $sendfile
		        	else
		        		echo "	Successfully done post processing EM step $_EM_" >> $sendfile
		        	fi
        
		        	if [[ $water_refine_Convert_file_ -ne 0 ]]; then
		        		echo "	Error in convertion for water refinement step $_water_refine_Convert_file_" >> $sendfile
		        	else
		        		echo "	Successfully converted for water refinement step $_water_refine_Convert_file_" >> $sendfile
		        	fi
        
		        	if [[ $water_refine_running_ -ne 0 ]]; then
		        		echo "	Error in running water refinement step in post processing $_water_refine_running_" >> $sendfile
		        	else
		        		echo "	Successfully ran water refinement step in post processing $_water_refine_running_" >> $sendfile
		        	fi
                	else
                		echo "ongoing-3"
                	fi
                done		
	#fi	
fi

