protein_old="1bf8,1d8v,1de3,1eza,1f16,1fpw,1p4s,1pbu,1tfb,1zxf,2k4t,2kok,2kty,2la3,2lfc,2lhf,2lvv,2m4k,2n8x,2o4e,2xfm,2ysz,demo"

protein_old="1bf8,1d8v,1de3,1eza,1f16,1fpw,1p4s,1pbu,1tfb,1zxf,2k4t,2kok,2kty,2la3,2lfc,2lhf,2lvv,2m4k,2n8x,2o4e,2xfm,2ysz,demo"

IFS=',' read -r -a protein_i <<< $protein_old

for i in ${protein_i[@]}
do
	mv /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/$i /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/$i_old

	if [ $? -ne 0 ]; then
		echo
		echo $i
		echo
 	fi
done
