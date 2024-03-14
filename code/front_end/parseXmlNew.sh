#!/bin/bash

xmlprotip="$1"

if [ -f $xmlprotip ]; then

	# ---
        proteinname=`xmllint --xpath "string(//protein_name)" "$xmlprotip"`
        seq_file=`xmllint --xpath "string(//seq_file)" "$xmlprotip"`
        upls_file=`xmllint --xpath "string(//upl_file)" "$xmlprotip"`
        hbond_file=`xmllint --xpath "string(//hbond_file)" "$xmlprotip"`
        ang_file=`xmllint --xpath "string(//ang_file)" "$xmlprotip"`
        protein_path=`xmllint --xpath "string(//protein_path)" "$xmlprotip"`
        
        hydrogen_omission=`xmllint --xpath "string(//hydrogen_omission)" "$xmlprotip"`
        aug_bounds=`xmllint --xpath "string(//aug_bounds)" "$xmlprotip"`
        eta_lo=`xmllint --xpath "string(//eta_lo)" "$xmlprotip"`
        eta_hi=`xmllint --xpath "string(//eta_hi)" "$xmlprotip"`
        include_neighbour=`xmllint --xpath "string(//include_neighbour)" "$xmlprotip"`
	multi_expand_k2=`xmllint --xpath "string(//k_2)" "$xmlprotip"`
	multi_expand_size_cutoff=`xmllint --xpath "string(//size_cutoff)" "$xmlprotip"`
	multi_expand_grp_min=`xmllint --xpath "string(//grp_min)" "$xmlprotip"`
	multi_incr_min=`xmllint --xpath "string(//incr_min)" "$xmlprotip"`
	grp_expand=`xmllint --xpath "string(//grp_expand)" "$xmlprotip"`
	
	# ---
	if [ -z "$proteinname" ]; then
		echo "$proteinname NOT found in $xmlprotip"
		exit -1
	fi
	if [ -z "$seq_file" ]; then
		echo "seq_file NOT found in $xmlprotip"
		exit -1
	fi
	if [ -z "$upls_file" ]; then
		echo "upls_file NOT found in $xmlprotip"
		exit -1
	fi
	if [ -z "$hydrogen_omission" ]; then
		echo "hydrogen_omission NOT found in $xmlprotip"
		exit -1
	fi
	if [ -z "$aug_bounds" ]; then
		echo "aug_bounds NOT found in $xmlprotip"
		exit -1
	fi
	if [ -z "$eta_lo" ]; then
		echo "eta_lo NOT found in $xmlprotip"
		exit -1
	fi
	if [ -z "$eta_hi" ]; then
		echo "eta_hi NOT found in $xmlprotip"
		exit -1
	fi
	if [ -z "$include_neighbour" ]; then
		echo "include_neighbour NOT found in $xmlprotip"
		exit -1
	fi			
	if [ -z "$multi_expand_k2" ]; then
		echo "multi_expand_k2 NOT found in $xmlprotip"
		exit -1
	fi
	if [ -z "$multi_expand_size_cutoff" ]; then
		echo "multi_expand_size_cutoff NOT found in $xmlprotip"
		exit -1
	fi
	if [ -z "$multi_expand_grp_min" ]; then
		echo "multi_expand_grp_min NOT found in $xmlprotip"
		exit -1
	fi	
	if [ -z "$multi_incr_min" ]; then
		echo "multi_incr_min NOT found in $xmlprotip"
		exit -1
	fi
	if [ -z "$grp_expand" ]; then
		echo "grp_expand NOT found in $xmlprotip"
		exit -1
	fi				
        # --- 	
        echo "INPUTS.protein_name = '$proteinname';"
        echo "INPUTS.seq_file = '$seq_file';"
        echo "INPUTS.upl_file = {'$upls_file'};"
        echo "INPUTS.hbond_file = '$hbond_file';"
        echo "INPUTS.ang_file = '$ang_file';"
        echo "INPUTS.protein_path = '$protein_path';"
        
        echo "PARAMS.hydrogen_omission = $hydrogen_omission;"
        echo "PARAMS.aug_bounds = $aug_bounds;"
        echo "PARAMS.eta_lo = $eta_lo;"
        echo "PARAMS.eta_hi = $eta_hi;"
        echo "PARAMS.include_neighbour = $include_neighbour;"
	echo "PARAMS.multi_expand_k2 = $multi_expand_k2;"
	echo "PARAMS.multi_expand_size_cutoff = $multi_expand_size_cutoff;"
	echo "PARAMS.multi_expand_grp_min = $multi_expand_grp_min;"
	echo "PARAMS.multi_incr_min = $multi_incr_min;"
	echo "PARAMS.grp_expand = $grp_expand;"
else
	echo "$xmlprotip NOT found"
	exit -1
fi
