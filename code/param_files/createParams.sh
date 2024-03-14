#!/bin/bash

printxmls() {

h_omi="$1"
aug_bounds="$2"
nu_lo="$3"
include_neigh="$4"
grp_exp="$5"

printf "<set_up>\n"
printf "\t<inputs>\n"
printf "\t\t<protein_name></protein_name>\n"
printf "\t\t<seq_file></seq_file>\n"
printf "\t\t<upl_file></upl_file>\n"
printf "\t\t<hbond_file></hbond_file>\n"
printf "\t\t<ang_file></ang_file>\n"
printf "\t\t<protein_path></protein_path>\n"
printf "\t</inputs>\n"
printf "\t<params>\n"
printf "\t\t<hydrogen_omission>%s</hydrogen_omission>\n" "$h_omi"
printf "\t\t<aug_bounds>%s</aug_bounds>\n" "$aug_bounds"
printf "\t\t<break_graph>\n"
printf "\t\t\t<eta_lo>%s</eta_lo>\n" "$nu_lo"
printf "\t\t\t<eta_hi>1</eta_hi>\n"
printf "\t\t</break_graph>\n"
printf "\t\t<include_neighbour>%s</include_neighbour>\n" "$include_neigh"
printf "\t\t<multi_expand>\n"
printf "\t\t\t<k_2>3</k_2>\n"
printf "\t\t\t<size_cutoff>0.1</size_cutoff>\n"
printf "\t\t\t<grp_min>250</grp_min>\n"
printf "\t\t\t<incr_min>2</incr_min>\n"
printf "\t\t</multi_expand>\n"
printf "\t\t<grp_expand>%s</grp_expand>\n" "$grp_exp"
printf "\t</params>\n"
printf "</set_up>\n"

}

param_file_dest="$1"

h_omi="1"
aug_bounds="0,6,6.5,7,7.5,8,8.5"
nu_lo="0.45,0.5,0.6,0.65,0.7,0.75,0.8"
include_neigh="4,5,6,7,8,10,11"
grp_exp="30,35,45,50,70,80,90"


IFS="," read -r -a h_omi_arr <<< "$h_omi"
IFS="," read -r -a aug_bounds_arr <<< "$aug_bounds"
IFS="," read -r -a nu_lo_arr <<< "$nu_lo"
IFS="," read -r -a include_neigh_arr <<< "$include_neigh"
IFS="," read -r -a grp_exp_arr <<< "$grp_exp"

count=0
for h_omi_ in "${h_omi_arr[@]}"
do
	for aug_bounds_ in "${aug_bounds_arr[@]}"
	do
		for nu_lo_ in "${nu_lo_arr[@]}"
		do
			for include_neigh_ in "${include_neigh_arr[@]}"
			do
				for grp_exp_ in "${grp_exp_arr[@]}"
				do
					count=$((count+1))			
					filename="param"$count".xml"
					printxmls $h_omi_ $aug_bounds_ $nu_lo_ $include_neigh_ $grp_exp_ > $param_file_dest/$filename
				done
			done
		done
	done
done


