#!/bin/bash

proteinname="$1"
seqfile="$2"
pdbfile="$3"

if [ ! -f $seqfile ]; then
	echo "$seqfile NOT found"
	exit -1
fi

if [ ! -f $pdbfile ]; then
	echo "$pdbfile NOT found"
	exit -1
fi

BAD_RESI_PATH=gap_correct/BAD_RESI

run_matlab_file=$proteinname"_outlier.m"

oppfileend=$proteinname"_outliers.txt"
echo "cd .." > $run_matlab_file
echo "addpath(genpath(pwd));" >> $run_matlab_file
echo "cd code;" >> $run_matlab_file
echo "outliers = getOutliers('$pdbfile')" >> $run_matlab_file
echo "filename = '$oppfileend'" >> $run_matlab_file
echo "fid = fopen(filename,'w');" >> $run_matlab_file
echo "for i=1:length(outliers)" >> $run_matlab_file
echo "	if i ~= length(outliers)" >> $run_matlab_file
echo "          fprintf(fid,'%d,',outliers(i))" >> $run_matlab_file
echo "	else" >> $run_matlab_file
echo "          fprintf(fid,'%d',outliers(i))" >> $run_matlab_file
echo "	end" >> $run_matlab_file
echo "end" >> $run_matlab_file

matlab -nodesktop -nosplash < $run_matlab_file
#########
outliers=`cat $proteinname"_outliers.txt"`
beginresi=`head -1 $seqfile | awk '{print $2}'`

#gap_correct/createModellerScript.sh proteinname seqfile pdbfile outliers beginresi
gap_correct/createBadResiModellerScript.sh $proteinname $seqfile $pdbfile $outliers $beginresi
#############################

modellergappdb_=`echo $pdbfile | awk -F'pdb' '{print $1}'`
modellergappdb=$modellergappdb_"gap.pdb"
alifile=$proteinname"_badresi.ali"
modelfile=$proteinname"_badresi.py"

mkdir "$BAD_RESI_PATH/$proteinname"
#cp $modellergappdb $MOD20_path/$proteinname/.
#cp $alifile $MOD20_path/$proteinname/.
#cp $modelfile $MOD20_path/$proteinname/.

mv $modellergappdb $BAD_RESI_PATH/$proteinname/.
mv $alifile $BAD_RESI_PATH/$proteinname/.
mv $modelfile $BAD_RESI_PATH/$proteinname/.

curr_path=`pwd`
cd $BAD_RESI_PATH/$proteinname
#mod9.23 $modelfile
mod10.0 $modelfile
cd $curr_path

##########################
mv $oppfileend $BAD_RESI_PATH/$proteinname/.
mv $run_matlab_file $BAD_RESI_PATH/$proteinname/.
mv $proteinname"_gaps" $BAD_RESI_PATH/$proteinname/.
