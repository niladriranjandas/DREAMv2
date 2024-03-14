#!/bin/bash

#############################
#expected template format
#<set_up>
#	<inputs>
#		<protein_name></protein_name>
#		<seq_file></seq_file>
#		<upl_file></upl_file>
#		<hbond_file></hbond_file>
#		<ang_file></ang_file>
#		<protein_path></protein_path>
#	</inputs>
#	<params>
#		<hydrogen_omission>1</hydrogen_omission>
#		<aug_bounds>7.5</aug_bounds>
#		<break_graph>
#			<eta_lo>0.65</eta_lo>
#			<eta_hi>1</eta_hi>
#		</break_graph>
#		<include_neighbour>11</include_neighbour>
#		<multi_expand>
#			<k_2>3</k_2>
#			<size_cutoff>0.1</size_cutoff>
#			<grp_min>250</grp_min>
#			<incr_min>2</incr_min>
#		</multi_expand>
#		<grp_expand>20</grp_expand>
#	</params>
#</set_up>
#
# e.g:
# ./runFirstStage.sh /home/nmr/Documents/our_algo/code/param_files/param1.xml /home/nmr/Documents/dj_prots/2k7k param1 2k7k	
##############################

usage(){
	echo "./runFirstStage.sh <protein_folder> <param_template> <destfolder>"
}

getfiles(){
	file_ext="$1" 
	protfolder="$2"
	destfolder="$3"	
	paramext="$4"

	#cp  $protfolder/*$file_ext $destfolder/.
	new_files=""
	tmp=`ls $protfolder/*$file_ext 2>/dev/null`
	if [[ $? -eq 0 ]]; then
		for f in `ls $protfolder/*$file_ext`
		do
			name_ext=`echo $f | awk -F'/' '{print $NF}'`
			name=`echo $name_ext | awk -F".$file_ext" '{print $1}'`
			new_name="$name"_"$paramext"."$file_ext"

			if [ $file_ext=="upl" ]; then
				ind=`echo $f | grep -i hbond`
				if [ -z $ind ]; then
					cp $f $destfolder/$new_name
					new_files=$new_files","$new_name		
				fi
			else
			   cp $f $destfolder/$new_name
			   new_files=$new_files","$new_name		
			fi
		done
	
		echo ${new_files:1:${#new_files}-1}
		#echo "${new_files::-1}"
	fi
}

makexml(){
	param_file="$1"
	protfolder="$2"
	destfolder="$3"
	paramext="$4"	
	protein_name="$5"

	printf "<set_up>"
	printf "\n\t<inputs>"
	prot_name=$protein_name"_"$paramext
	printf "\n\t\t<protein_name>%s</protein_name>" $prot_name
	seq_file=$(getfiles seq $protfolder $destfolder $paramext)
	printf "\n\t\t<seq_file>%s</seq_file>" $seq_file
	upl_file=$(getfiles upl $protfolder $destfolder $paramext)
	h_file=$(getfiles "hBond*upl" $protfolder $destfolder $paramext)	
	printf "\n\t\t<upl_file>%s</upl_file>" $upl_file
	printf "\n\t\t<hbond_file>%s</hbond_file>" $h_file
	aco_file=$(getfiles "aco" $protfolder $destfolder $paramext)
	printf "\n\t\t<ang_file>%s</ang_file>" $aco_file	
	printf "\n\t\t<protein_path>../protein/%s</protein_path>" $prot_name
	printf "\n\t</inputs>"

	printf "\n\t<params>"	
	h_omm=`grep hydrogen_omission $param_file`
	printf "\n\t\t%s>" $h_omm
	aug_bounds=`grep aug_bounds $param_file`
	printf "\n\t\t%s" $aug_bounds
	printf "\n\t\t<break_graph>"
	eta_lo=`grep eta_lo $param_file`
	printf "\n\t\t\t%s" $eta_lo
	eta_hi=`grep eta_hi $param_file`
	printf "\n\t\t\t%s" $eta_hi
	printf "\n\t\t</break_graph>"
	incl_neigh=`grep include_neighbour $param_file`
	printf "\n\t\t%s" $incl_neigh
	printf "\n\t\t<multi_expand>"
	k2=`grep k_2 $param_file`
	printf "\n\t\t\t%s" $k2
	size_cutoff=`grep size_cutoff $param_file`
	printf "\n\t\t\t%s" $size_cutoff
	grp_min=`grep grp_min $param_file`
	printf "\n\t\t\t%s" $grp_min
	incr_min=`grep incr_min $param_file`
	printf "\n\t\t\t%s" $incr_min
	printf "\n\t\t</multi_expand>"
	expnd_grp=`grep grp_expand $param_file`
	printf "\n\t\t%s" $expnd_grp
	printf "\n\t</params>"	
	printf "\n</set_up>"
}

runStage1prepStage2(){
	#----------- i/p ---------#
	param_file="$1"
	protfolder="$2"
	destfolder="$3"
	paramext="$4"
	protein_name="$5"

	# ---------- xml -------- #
	xml_file=$protein_name"_"$paramext".xml"
	makexml $param_file $protfolder $destfolder $paramext $protein_name  > $destfolder/$xml_file
	#makexml $param_file $protfolder $destfolder $paramext $protein_name > $destfolder/$xml_file
}

runStage1genStage2(){
	param_file="$1"
	protfolder="$2"
	paramext="$3"
	protein_name="$4"
	newproteinfolder="../protein"		
	
	destfolder=$newproteinfolder"/"$protein_name"_"$paramext
	mkdir $destfolder
	runStage1prepStage2 $param_file	$protfolder $destfolder	$paramext $protein_name
	xmlfile=`find $destfolder -iname "*xml"`
	
	prot_xml=`echo $xmlfile | awk -F'/' '{print $NF}' | awk -F'.xml' '{print $0}'`
	prot=`echo $prot_xml | awk -F'.xml' '{print $1}'`
	echo "./runBatch.sh -list $prot > /dev/null &"
	echo "# run archive after 2nd stage run"
	echo "./archiveProteinRun.sh $prot ../archive/"
	# change runBatch to call genRunScripts_bckend_new.sh
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    runStage1genStage2 "$@"  # Call the main function if the script is in shell
fi


