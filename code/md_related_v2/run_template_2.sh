$GMX mdrun -deffnm "$work_folder/$protein_name"_em
#$GMX trjconv -s "$work_folder/$protein_name"_em.tpr -f "$work_folder/$protein_name"_em.gro -o "$work_folder/$protein_name"_em.pdb <<< 0
#egrep -v 'SOL|NA' "$work_folder/$protein_name"_em.pdb > "$work_folder/$protein_name"_em_nosol.pdb
$GMX trjconv -f "$work_folder/$protein_name"_em.gro -s "$work_folder/$protein_name"_em.tpr -pbc mol -ur compact -o "$work_folder/$protein_name"_em_nosol.pdb <<< 1
