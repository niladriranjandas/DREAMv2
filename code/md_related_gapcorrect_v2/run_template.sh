


$GMX pdb2gmx -f "$work_folder/$protein_name".pdb -water tip3p -o "$work_folder/$protein_name"_process.gro -ignh <<< 8

$GMX editconf -f "$work_folder/$protein_name"_process.gro -o "$work_folder/$protein_name"_newbox.gro -c -d 1.0 -bt cubic
$GMX solvate -cp "$work_folder/$protein_name"_newbox.gro -cs spc216.gro -o "$work_folder/$protein_name"_solv.gro -p $work_folder/topol.top
$GMX grompp -f $mdp_folder/ions.mdp -c "$work_folder/$protein_name"_solv.gro -p $work_folder/topol.top -o $work_folder/ions.tpr -maxwarn 1
$GMX genion -s $work_folder/ions.tpr -o "$work_folder/$protein_name"_solv_ions.gro -p $work_folder/topol.top -pname NA -nname CL -neutral <<< 13

$GMX grompp -f $mdp_folder/minim.mdp -c "$work_folder/$protein_name"_solv_ions.gro -p $work_folder/topol.top -o "$work_folder/$protein_name"_em.tpr -r "$work_folder/$protein_name"_solv_ions.gro
$GMX mdrun -deffnm "$work_folder/$protein_name"_em
#$GMX trjconv -s "$work_folder/$protein_name"_em.tpr -f "$work_folder/$protein_name"_em.gro -o "$work_folder/$protein_name"_em.pdb <<< 0
#egrep -v 'SOL|NA' "$work_folder/$protein_name"_em.pdb > "$work_folder/$protein_name"_em_nosol.pdb
$GMX trjconv -f "$work_folder/$protein_name"_em.gro -s "$work_folder/$protein_name"_em.tpr -pbc mol -ur compact -o "$work_folder/$protein_name"_em_nosol.pdb <<< 1

