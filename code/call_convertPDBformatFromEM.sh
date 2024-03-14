#!/bin/bash

if [ $# -ne 2 ]; then
	echo "call_convertPDBformatFromEM 2 args needed"
	exit -1
fi

ip_pdb_file="$1"
#map_file_em="map_names_em.csv"
map_file_em="map_names_em_gro21.csv"
op_pdb_file="$2"

#python convertPDBformatFromEM.py "$ip_pdb_file" "$map_file_em" "$op_pdb_file"
python3 convertPDBformatFromEM.py "$ip_pdb_file" "$map_file_em" "$op_pdb_file"

sed -i 's/CD1  ILE/CD1 ILE/g' $op_pdb_file

sed -i 's/ HD11 ILE/HD11 ILE/g' $op_pdb_file
sed -i 's/ HD12 ILE/HD12 ILE/g' $op_pdb_file
sed -i 's/ HD13 ILE/HD13 ILE/g' $op_pdb_file

sed -i 's/HG SER/HG  SER/g' $op_pdb_file
sed -i 's/HG CYS/HG  CYS/g' $op_pdb_file

tmp_file=$op_pdb_file"_tmp"
# sed -i "/)/d" $op_pdb_file #is slower than grep
grep -v ")" $op_pdb_file > $tmp_file; mv $tmp_file $op_pdb_file
