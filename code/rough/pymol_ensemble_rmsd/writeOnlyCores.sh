# ./writeOnlyCores.sh  /home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/pymol_ensemble_rmsd/2m4k_cluster_0002.pdb /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/2m4k/output_pdb/grp_registration_2m4k.pdb

genPDB="$1"
grpregPDB="$2"

resi=`cat $grpregPDB | awk '{print $5}' | uniq`

while IFS='' read -r line || [[ -n "$line" ]]; do
	resi_i=`echo $line | awk '{print $5}'`
	if [[ ! -z $resi_i ]]; then
		abc=`echo $resi | grep $resi_i`
		if [[ ! -z $abc ]]; then
			echo "$line"
		fi
	fi
done < "$genPDB"
