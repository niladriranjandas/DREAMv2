function resi_full_sidec = getResiFullSidec(grpreg_atom_index, Comp, resi_noH_details)
%%
%
%
%
%%
   if nargin ~= 3
      error('Module: getResiFullSidec: incorrect no. of inputs.');
   end
%% 
  grpreg_Comp_residue = Comp.residue(grpreg_atom_index);
  grpreg_resi         = unique(grpreg_Comp_residue);
  
    resi_full_sidec = nan(1,length(grpreg_resi));
    for i=1:length(grpreg_resi)
        Comp_num_seq_indx =    find(Comp.num_seq == grpreg_resi(i));
        resi_no = Comp.seq(Comp_num_seq_indx);
        grpreg_index_resii = find(grpreg_Comp_residue == grpreg_resi(i));
           tmp_index = grpreg_atom_index(grpreg_index_resii);
       resi_full_sidec(i) = (i_has_full_sidec(Comp.atom_names(tmp_index),resi_noH_details{resi_no,2})) * grpreg_resi(i);
    end
end

function one_zero = i_has_full_sidec(grpreg_atom_names, resi_i_noH)
 
   one_zero = 0;
   resi_i_noH = resi_i_noH(3:end-1);  %starts with C,O,N....,C,O,N
   
   if length(grpreg_atom_names) < length(resi_i_noH)
       return
   end
   
   for i=1:length(resi_i_noH)
      if ~ismember(resi_i_noH(i), grpreg_atom_names)
          return
      end
   end
   one_zero =1 ;
end